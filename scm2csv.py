import os
import sys
import pandas as pd
import opencog.bioscience
from opencog.atomspace import AtomSpace, types
from opencog.scheme import scheme_eval
from config import RESULT_DIR
import time
import json

def find_txt_name(node, atomspace):
    # Finds the text name of the bio-entity
    name = \
    '''
    (let ([name (cog-outgoing-set (cog-execute! (GetLink
        (Evaluation (PredicateNode "has_name") (ListLink {} (Variable "$name"))))))])
        (if (null? name) "" (cog-name (car name))))
    '''.format(node)
    name = scheme_eval(atomspace, name).decode("utf-8")
    return name if len(name) > 1 else None

def find_inherits_or_members(node, atomspace, linktype="Inheritance"):
    # Returns concepts where the given bio-entity inherits or is a member from
    # and their count by type
    if linktype == "Inheritance":
        atoms = atomspace.get_atoms_by_type(types.InheritanceLink)
    else:
        atoms = atomspace.get_atoms_by_type(types.MemberLink)
    _from = []
    has_ = []
    count_by_type = {}
    for i in atoms:
        element1 = i.out[0]
        element2 = i.out[1]
        if element1 == node:
            txt_name = find_txt_name(element2, atomspace)
            node_name = element2.name
            _from.append("{} ({})".format(node_name, txt_name))
            # Count member_from by type for summary
            elm2_type = element2.type_name
            if elm2_type in count_by_type.keys():
                count_by_type[elm2_type] = count_by_type[elm2_type] + 1
            else:
                count_by_type[elm2_type] = 1
        elif element2 == node:
            txt_name = find_txt_name(element1, atomspace)
            node_name = element1.name
            has_.append("{} ({})".format(node_name, txt_name))
    _from = "\n".join(_from) if _from else None
    has_ = "\n".join(has_) if has_ else None
    return _from, has_, count_by_type


def find_interactions(node, atomspace):
    # Finds and returns interactions of the given bio-entity with others, the type of interaction
    # and their count
    interaction_result = \
    '''
    (let ([interactions (list "expresses" "interacts_with" "binding" "reaction" "inhibition" "activation" "expression" "catalysis" "ptmod" 
           "GO_regulates" "GO_positively_regulates" "GO_negatively_regulates" "has_part" "has_role" "translated_to" "transcribed_to" )]
          [node {}])
        (string-join (map (lambda (i)
          (let ([result   
            (cog-outgoing-set (cog-execute! (GetLink (ChoiceLink 
            (Evaluation (PredicateNode i) (ListLink node (Variable "$n")))
            (Evaluation (PredicateNode i) (ListLink (Variable "$n") node))
            (Evaluation (PredicateNode i) (SetLink node (Variable "$n")))))))])
          (if (null? result) ""
          (string-join (list i (string-join (map (lambda (x) (cog-name x))result) ",")) ":")))
        )interactions) " "))
    '''.format(node)
    interaction_result = scheme_eval(atomspace, interaction_result).decode("utf-8")
    interaction_result = "\n".join(interaction_result.split())
    count = {}
    if interaction_result:
        for i in interaction_result.split("\n"):
            interaction = i.split(":")
            count[interaction[0]] = len(interaction[1].split(","))
    return interaction_result, count

def find_locations(node, atomspace):
    # Finds cellular location of the given bio-entity
    locations = \
    '''
    (let ([loc (cog-outgoing-set (cog-execute! (GetLink
        (Evaluation (PredicateNode "has_location") (ListLink {} (Variable "$loc"))))))])
        (if (null? loc) "" (string-join (map (lambda (l) (cog-name l)) loc) ",")))
    '''.format(node)
    locations = scheme_eval(atomspace,locations).decode("utf-8")
    return locations if locations else None

def filter_df(df):
    # Drop empty columns from the dataframe
    empty_cols = [col for col in df.columns if df[col].isnull().all()]
    df.drop(empty_cols, axis=1, inplace=True)
    # Fill null values with NA
    filtered_df = df.fillna("N/A")
    filtered_df = filtered_df.drop_duplicates()
    return filtered_df

def generate_url(node, node_type):
    # Generate the source URL for the bio-entity
    if "Uniprot" in node_type:
        return "https://www.uniprot.org/uniprot/{}".format(node)
    elif "Gene" in node_type:
        return "https://www.ncbi.nlm.nih.gov/gene/?term={}".format(node)
    elif node_type == "Reactome":
        return "http://www.reactome.org/content/detail/{}".format(node)
    elif node_type in ["CellularComponent", "BiologicalProcess","MolecularFunction"]:
        return "http://amigo.geneontology.org/amigo/term/{}".format(node)
    elif node_type == "Chebi":
        return "https://www.ebi.ac.uk/chebi/searchId.do?chebiId={}".format(node)
    elif node_type == "Smp":
        return "http://smpdb.ca/view/{}".format(node)
    else:
        return node_type.replace("Node"," database")

def to_csv(file_dir, main_nodes=False):
    start = time.time()
    if main_nodes:
        main_nodes = [i["geneName"] for i in main_nodes]
    path = os.path.join(RESULT_DIR, file_dir)
    input_file = os.path.join(path, "result.scm")
    atomspace = AtomSpace()
    modules = "(use-modules (opencog) (opencog bioscience) (opencog persist-file) (opencog exec))"
    scheme_eval(atomspace, modules)
    scheme_eval(atomspace, "(load-file \"{}\")".format(input_file))
    df_columns = ["Bio-entity","Type","Name","Source db","Member_of","Has_members","Inherits_from","Has_inherits","Interaction/type","Cell Location"]
    df1 = pd.DataFrame([], columns=df_columns)
    df2 = pd.DataFrame([], columns=df_columns)

    summary = dict()
    cross_an = dict()
    main_input = dict()

    # Find all Biological entities
    bio_types = scheme_eval(atomspace, "(cog-get-types)").decode("utf-8")
    bio_types = bio_types.replace("(","").replace(")","").split(" ")
    bio_types = [i for i in bio_types if "Node" in i]
    molecules = ["Uniprot", "Gene", "Chebi", "Pubchem","Drubank","Refseq","Enst"]
    atoms = []
    atoms_count = [{}] 
    for t in bio_types:
        try:
            atom = atomspace.get_atoms_by_type(getattr(types, t))
            if len(atom) > 0 and not t in ["Node","PredicateNode","ConceptNode"]:
                atoms_count[0].update({t.replace("Node",""):len(atom)})
                atoms = atoms + atom
        except:
            continue

    for i,val in enumerate(atoms):
        txt_name = find_txt_name(val, atomspace)
        if txt_name:
           node = val.name
           node_type = val.type_name.replace("Node","")
           source = generate_url(node, node_type)
           count = {}
           member_of, has_members, count_mem = find_inherits_or_members(val, atomspace,linktype="member")
           count.update(count_mem)
           inherits_from, has_inherits, count_inh  = find_inherits_or_members(val, atomspace)
           count.update(count_inh)
           location = find_locations(val, atomspace)
           if node_type in molecules:
               interacts_with, count_int = find_interactions(val, atomspace)
               count.update(count_int)
               if main_nodes and node in main_nodes:
                  main_input[node] = [count]
               else:  
                  cross_an[node] = [count] 
               df1.loc[len(df1)] = [node, node_type, txt_name, source, member_of,has_members, inherits_from,has_inherits,interacts_with,location]
           else:
               interacts_with = None 
               df2.loc[len(df2)] = [node, node_type, txt_name, source, member_of,has_members, inherits_from,has_inherits,interacts_with,location]
    
    summary["A Reference Databases"] = "mozi.ai/datasets/"
    summary["Cross Annotations"] = cross_an
    summary["Input Genes"] = main_input
    summary["Total count"] = {"Count": atoms_count}
    with open(os.path.join(path, "summary.json"), "w") as j:
        json.dump(summary, j)
    
    df1 = filter_df(df1)
    if len(df1) > 0:
        df1.to_csv(os.path.join(path, "result1.csv"), index=False)
    df2 = filter_df(df2)
    if len(df2) > 0:
        df2.to_csv(os.path.join(path, "result2.csv"), index=False)
    end = time.time()
    
    return "Time to parse atomese to csv and generate summary: {}".format(end-start)

if __name__ == "__main__":
    print(to_csv(sys.argv[1]))