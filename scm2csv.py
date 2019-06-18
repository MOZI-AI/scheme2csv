import sys
import re
import pandas as pd
import uuid
import os

CSV_FOLDER = "/root/csv_result"

def find_name(str):
	name = re.findall('"([^"]*)"', str)
	if len(name) > 0:
		return re.findall('"([^"]*)"', str)[0]
	return ""
def checkdic(dic, str):
	try:
		return dic[str]
	except KeyError:
		return ""
def uni(stri):
	return ','.join(set(stri.split(','))).replace('Uniprot:', '').replace('pubmed:', '')

def find_codingGene(prot, dec):
	for k in dec.keys():
		if prot in dec[k].split(","):
			return k

def to_csv(file):
	member = []
	evalun = []
	interaction = []
	go_annotation = 0
	gene_pathway = 0
	biogrid = 0
	userid = str(uuid.uuid4())
	result = []

	f = open(file, "r")
	lines=open(file, "r").readlines()

	for num, line in enumerate(f , 0):
		if "InheritanceLink" in line or "MemberLink" in line:
			member.append(num)
		elif "EvaluationLink" in line:
			evalun.append(num)
		elif "gene-go-annotation" in line:
			go_annotation = num
		elif "gene-pathway-annotation" in line:
			gene_pathway = num
		elif "biogrid-interaction-annotation" in line:
			biogrid = num

	GO_ns = {}
	node_name = {}
	node_defn = {}
	pubmed ={}
	location = {}
	express_pw = {}
	expresses_bg = {}

	for i in evalun:
		if "GO_namespace" in lines[i+1]:
			GO_ns.update({find_name(lines[i+3]): find_name(lines[i+4])})
		elif "has_name" in lines[i+1]:
			node_name.update({find_name(lines[i+3]): find_name(lines[i+4])})
		elif "has_definition" in lines[i+1]:
			node_defn.update({find_name(lines[i+3]): find_name(lines[i+4])})
		elif "interacts_with" in lines[i+1]:
			interaction.append(i+1)
		elif "expresses" in lines[i+1]:
			if biogrid != 0:
				if i+1 < biogrid:
					express_pw.update({find_name(lines[i+3]): checkdic(express_pw, find_name(lines[i+3])) + find_name(lines[i+4]) + ','})
				else:
					expresses_bg.update({find_name(lines[i+3]): checkdic(expresses_bg, find_name(lines[i+3])) + find_name(lines[i+4]) + ','})
			else:
				express_pw.update({find_name(lines[i+3]): checkdic(express_pw, find_name(lines[i+3])) + find_name(lines[i+4]) + ','})
		elif "has_location" in lines[i+1]:
			location.update({find_name(lines[i+3]): find_name(lines[i+4])})
		elif "has_pubmed" in lines[i+1]:
			pubmedline = i+11
			while "pubmed" in lines[pubmedline]:
				pubmed.update({find_name(lines[i+6]) + find_name(lines[i+7]): checkdic(pubmed, find_name(lines[i+6]) + find_name(lines[i+7])) + find_name(lines[pubmedline]) + ','})
				pubmedline +=1
	
	col = ["Gene_ID","GO_cellular_componenet", "GO_Molecular_function", "GO_Biological_process", "pathway", "proteins", "small_mol"]
	gene_go = pd.DataFrame([], columns= col)
	for i in member:
		gene = find_name(lines[i+1])
		go = find_name(lines[i+2])
		if "GeneNode" in lines[i+1] and "GO:" in lines[i+2] and GO_ns != {}:
			if GO_ns[go] == "cellular_component":
				gene_go = gene_go.append(pd.DataFrame([[gene, go +" (" +node_name[go] +")", "", "", "", "", ""]], columns=col))
			elif GO_ns[go] == "biological_process":
				gene_go = gene_go.append(pd.DataFrame([[gene, "", "", go +" (" +node_name[go] +")", "", "", ""]], columns=col))	
			elif GO_ns[go] == "molecular_function":
				gene_go = gene_go.append(pd.DataFrame([[gene, "", go +" (" +node_name[go] +")", "", "", "", ""]], columns=col))	
		elif "GeneNode" in lines[i+1] and ("R-HSA" in lines[i+2] or 'ConceptNode "SMP' in lines[i+2]):
			prot = []
			sm = []
			pathway = find_name(lines[i+2])
			for j in member:
				if 'ChEBI:' in lines[j+1] and find_name(lines[j+2]) == pathway:
					sm.append(find_name(lines[j+1]) + ' (' + node_defn[find_name(lines[j+1])] + ")" )
				elif 'Uniprot:' in lines[j+1] and find_name(lines[j+2]) == pathway:
					prot.append(find_name(lines[j+1]) + ' ('+ node_defn[find_name(lines[j+1])] + ")" )
			gene_go = gene_go.append(pd.DataFrame([[gene, "","", "", pathway, str(prot), str(sm)]], columns=col))
		elif "Uniprot:" in lines[i+1] and ("R-HSA" in lines[i+2] or 'ConceptNode "SMP' in lines[i+2]):
			gene = find_codingGene(find_name(lines[i+1]), express_pw)
			prot = []
			sm = []
			pathway = find_name(lines[i+2])
			for j in member:
				if 'ChEBI:' in lines[j+1] and find_name(lines[j+2]) == pathway:
					sm.append(find_name(lines[j+1]) + ' (' + node_defn[find_name(lines[j+1])] + ")" )
				elif 'Uniprot:' in lines[j+1] and find_name(lines[j+2]) == pathway:
					prot.append(find_name(lines[j+1]) + ' ('+ node_defn[find_name(lines[j+1])] + ")" )
			gene_go = gene_go.append(pd.DataFrame([[gene, "","", "", pathway, str(prot), str(sm)]], columns=col))
	
	# Gene Go annotation
	if go_annotation != 0:
		col_go = ["Gene_ID[location]", "Gene_Name", "Gene_definition", "GO_cellular_componenet", "GO_Molecular_function", "GO_Biological_process"]
		go_df = pd.concat([pd.DataFrame([[g + "[" +checkdic(location, g) + "]", node_name[g], node_defn[g], "\n".join(filter(None, gene_go[gene_go['Gene_ID'] == g]['GO_cellular_componenet'].get_values())), 
		"\n".join(filter(None, gene_go[gene_go['Gene_ID'] == g]['GO_Molecular_function'].get_values())),
		"\n".join(filter(None, gene_go[gene_go['Gene_ID'] == g]['GO_Biological_process'].get_values()))]], columns= col_go) for g in set(gene_go['Gene_ID'])], ignore_index=True)
		go_df = go_df.dropna(axis=1, how='all')
		go_df.to_csv(os.path.join(CSV_FOLDER,userid+"-Gene_GO_annotation.csv"))
		result.append({"displayName":"GO" ,"fileName":  userid+"-Gene_GO_annotation.csv"})
	
	# Gene Pathway annotation
	if gene_pathway != 0:
		col_pw = ["Pathway", "Pathway_detail", "Gene[Uniprot]", "Proteins", "Small Molecules"]
		pw_df = pd.concat([pd.DataFrame([[p, checkdic(node_name, p) + '[' + checkdic(node_defn, p) + ']',",".join(set(gene_go[gene_go['pathway'] == p]['Gene_ID'].get_values())), ",".join(set(gene_go[gene_go['pathway'] == p]['proteins'].get_values())), 
		"\n".join(set(filter(None, gene_go[gene_go['pathway'] == p]['small_mol'].get_values())))]], 
		columns= col_pw) for p in set(filter(None, gene_go['pathway']))], ignore_index=True)
		pw_df.to_csv(os.path.join(CSV_FOLDER,userid+"-Gene_pathway_annotations.csv"))
		result.append({"displayName":"PATHWAY" ,"fileName": userid+"-Gene_pathway_annotations.csv"})
	
	# Biogrid annotation
	if biogrid != 0:
		col_int = ["Interactor-1","Uniprot", "Interactor-1_Name", "Interactor-1_definition", "Interaction", "Interactor-2", "Uniprot", "Interactor-2_Name", "Interactor-2_definition", "Pubmed ID"]
		bg_df = pd.concat([pd.DataFrame([[find_name(lines[i+2]), uni(checkdic(expresses_bg, find_name(lines[i+2]))), checkdic(node_name, find_name(lines[i+2])), checkdic(node_defn, find_name(lines[i+2])), "Interacts_with", 
		find_name(lines[i+3]), uni(checkdic(expresses_bg, find_name(lines[i+3]))), checkdic(node_name, find_name(lines[i+3])), checkdic(node_defn, find_name(lines[i+3])), uni(checkdic(pubmed, find_name(lines[i+2])+find_name(lines[i+3])))]], columns= col_int) for i in interaction], ignore_index=True)
		bg_df.to_csv(os.path.join(CSV_FOLDER,userid+"-biogrid_annotation.csv"))
		result.append({"displayName":"BIOGRID" ,"fileName": userid+"-biogrid_annotation.csv"})

	return result

if __name__ == "__main__":
	to_csv(sys.argv[1])
