import sys
from collections import defaultdict
import re
import pandas as pd
import uuid
import os

CSV_FOLDER = "/root/csv_result"

def find_name(str):
	name = re.findall('"([^"]*)"', str)
	if len(name) > 0 and "VariableNode" not in str:
		return re.findall('"([^"]*)"', str)[0]
	return ""
def checkdic(dic, str):
	try:
		return dic[str]
	except KeyError:
		return ""
def find_codingGene(prot, dec):
	for k in dec.keys():
		if prot in dec[k].split(",") and prot != "":
			return k
	return ""
def flatten_list(lists):
	flattened_list = [y for x in lists for y in x]
	return flattened_list
def find_go(lst, identifier='ID'):
	result = [i.split(" (") for i in lst]
	if identifier is 'name':
		return [n[1].replace(')', '') for n in result]
	else:
		return [n[0].replace(')', '') for n in result]

def to_csv(file):
	member = []
	evalun = []
	interaction = []
	go_annotation = 0
	gene_pathway = 0
	main_nodes = 0
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
		elif 'ConceptNode "main"' in line:
			main_nodes = num
		elif "gene-go-annotation" in line:
			go_annotation = num
		elif "gene-pathway-annotation" in line:
			gene_pathway = num
		elif "biogrid-interaction-annotation" in line:
			biogrid = num
	main_genes = [find_name(i) for i in lines[main_nodes:min([x for x in [go_annotation, gene_pathway,biogrid] if x!=0])] if "GeneNode" in i]
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
				if i+1 < biogrid and gene_pathway < biogrid:
					express_pw.update({find_name(lines[i+3]): checkdic(express_pw, find_name(lines[i+3])) + find_name(lines[i+4]) + ','})
				elif i+1 > biogrid and gene_pathway > biogrid:
					express_pw.update({find_name(lines[i+3]): checkdic(express_pw, find_name(lines[i+3])) + find_name(lines[i+4]) + ','})
				else:
					expresses_bg.update({find_name(lines[i+3]): checkdic(expresses_bg, find_name(lines[i+3])) + find_name(lines[i+4]) + ','})
			else:
				express_pw.update({find_name(lines[i+3]): checkdic(express_pw, find_name(lines[i+3])) + find_name(lines[i+4]) + ','})
		elif "has_location" in lines[i+1]:
			if find_name(lines[i+3]) in location.keys():
				location.update({find_name(lines[i+3]): location[find_name(lines[i+3])]+[find_name(lines[i+4])]})
			location.update({find_name(lines[i+3]): [find_name(lines[i+4])]})
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
					sm.append(find_name(lines[j+1])  )
				elif 'Uniprot:' in lines[j+1] and find_name(lines[j+2]) == pathway:
					prot.append(find_name(lines[j+1])  )
			gene_go = gene_go.append(pd.DataFrame([[gene, "","", "", pathway, tuple(prot), tuple(sm)]], columns=col))
		elif "Uniprot:" in lines[i+1] and ("R-HSA" in lines[i+2] or 'ConceptNode "SMP' in lines[i+2]):
			gene = find_codingGene(find_name(lines[i+1]), express_pw)
			prot = []
			sm = []
			pathway = find_name(lines[i+2])
			for j in member:
				if 'ChEBI:' in lines[j+1] and find_name(lines[j+2]) == pathway:
					sm.append(find_name(lines[j+1]) )
				elif 'Uniprot:' in lines[j+1] and find_name(lines[j+2]) == pathway:
					prot.append(find_name(lines[j+1]) )
			gene_go = gene_go.append(pd.DataFrame([[gene, "","", "", pathway, tuple(prot), tuple(sm)]], columns=col))
	# Gene GO annotation
	if go_annotation != 0:
		go_genes_list = set(main_genes)
		gene_description = [checkdic(node_name, g) for g in go_genes_list]
		namespaces = ['GO_Molecular_function', 'GO_Biological_process', 'GO_cellular_componenet' ]
		namespace_details = ['Name', 'ID']
		go_data = []
		col_length = len(namespaces)*len(namespace_details)
		go_column_arrays = [flatten_list([[i]*col_length for i in go_genes_list]),
						flatten_list([[i]*col_length for i in gene_description]),
						flatten_list([[n]*len(namespace_details) for n in namespaces]*len(go_genes_list)),
						flatten_list([namespace_details]*len(namespaces)*len(go_genes_list))]
		for g in go_genes_list:
			go_data.append(find_go(set(filter(None,gene_go[gene_go['Gene_ID'] == g]['GO_Molecular_function'].get_values())), 'name'))
			go_data.append(find_go(set(filter(None,gene_go[gene_go['Gene_ID'] == g]['GO_Molecular_function'].get_values()))))
			go_data.append(find_go(set(filter(None,gene_go[gene_go['Gene_ID'] == g]['GO_Biological_process'].get_values())), 'name'))
			go_data.append(find_go(set(filter(None,gene_go[gene_go['Gene_ID'] == g]['GO_Biological_process'].get_values()))))
			go_data.append(find_go(set(filter(None,gene_go[gene_go['Gene_ID'] == g]['GO_cellular_componenet'].get_values())), 'name'))
			go_data.append(find_go(set(filter(None,gene_go[gene_go['Gene_ID'] == g]['GO_cellular_componenet'].get_values()))))
		if go_data:
			index_length = max([len(i) for i in go_data])
			go_data = [i+['']*(index_length- len(i)) for i in go_data] 
			go_data_lst = []
			for d in range(index_length):
				try:
					go_data_lst.append([i[d] for i in go_data])
				except IndexError:
					continue
			cols = pd.MultiIndex.from_arrays(go_column_arrays, names=('Gene', 'Description', 'Features', 'Details'))
			go_df = pd.DataFrame(go_data_lst, columns=cols)
			go_df.to_csv(os.path.join(CSV_FOLDER,userid+"-Gene_GO_annotation.csv"))
			go_df.to_excel(os.path.join(CSV_FOLDER,userid+"-Gene_GO_annotation.csv.xlsx"), sheet_name='GO-annotation')
			result.append({"displayName":"GO" ,"fileName":  userid+"-Gene_GO_annotation.csv"})
	# Gene Pathway annotation
	if gene_pathway != 0:
		pathways_list = set(filter(None, gene_go['pathway']))
		pathways_description = [checkdic(node_name, p) for p in pathways_list]
		features = ['Gene', 'Protein', 'Small molecule']
		pathway_data = []
		column_arrays = [flatten_list([[i]*len(features) for i in pathways_list]),
						flatten_list([[i]*len(features) for i in pathways_description]),
						flatten_list([features]*len(pathways_list))]
		for p in pathways_list:
			pathway_data.append(list(filter(None, set(gene_go[gene_go['pathway'] == p]['Gene_ID'].get_values()))))			
			pathway_data.append([p + "  " +str(find_codingGene(p, express_pw)) for p in list(set(gene_go[gene_go['pathway'] == p]['proteins'].get_values()[0]))])
			pathway_data.append(list(set(gene_go[gene_go['pathway'] == p]['small_mol'].get_values()[0])))
		if pathway_data:
			index_length = max([len(i) for i in pathway_data])
			pathway_data = [i+['']*(index_length- len(i)) for i in pathway_data] 
			pathway_data_lst = []
			for d in range(index_length):
				try:
					pathway_data_lst.append([i[d] for i in pathway_data])
				except IndexError:
					continue

			cols = pd.MultiIndex.from_arrays(column_arrays, names=('Pathway', 'Name', 'Feature'))
			pathway_df = pd.DataFrame(pathway_data_lst, columns=cols)
			pathway_df.to_csv(os.path.join(CSV_FOLDER,userid+"-Gene_pathway_annotations.csv"))
			pathway_df.to_excel(os.path.join(CSV_FOLDER,userid+"-Gene_pathway_annotations.xlsx"), sheet_name='pathway-annotation')
			result.append({"displayName":"PATHWAY" ,"fileName": userid+"-Gene_pathway_annotations.csv"})
	# Biogrid annotation
	if biogrid != 0:
		gene_interactions = defaultdict(list)
		interaction = [(find_name(lines[i+2]),find_name(lines[i+3])) for i in interaction]
		for key, val in interaction:
			if not checkdic(pubmed,key+val) == "":
				gene_interactions[key].append(val)
		gene_list = list(gene_interactions.keys())
		gene_description = [checkdic(node_name, g) for g in gene_list]
		features = ['Location', 'Proteins', 'Interacting_genes', 'PMID' ]
		biogrid_data = []
		col_length = 0
		biogrid_column_arrays = [flatten_list([[i]*len(features) for i in gene_list]),
						flatten_list([[i]*len(features) for i in gene_description]),
						flatten_list([features]*len(gene_list))]
		for g in gene_list:
			if checkdic(location, g) == "":
				biogrid_data.append([])
			else:
				biogrid_data.append(checkdic(location, g))			
			biogrid_data.append([p for p in (checkdic(expresses_bg, g).split(','))])
			biogrid_data.append(gene_interactions[g])
			biogrid_data.append([",\n".join(checkdic(pubmed,g+i).split(",")) if checkdic(pubmed,g+i) != "" else ",\n".join(checkdic(pubmed,g+i).split(","))  
			for i in gene_interactions[g]])
		if biogrid_data:
			index_length = max([len(i) for i in biogrid_data])
			biogrid_data = [i+['']*(index_length- len(i)) for i in biogrid_data] 
			biogrid_data_lst = []
			for d in range(index_length):
				try:
					biogrid_data_lst.append([i[d] for i in biogrid_data])
				except IndexError:
					continue

			cols = pd.MultiIndex.from_arrays(biogrid_column_arrays, names=('Gene', 'Description', 'Feature'))
			biogrid_df = pd.DataFrame(biogrid_data_lst, columns=cols)
			biogrid_df.to_csv(os.path.join(CSV_FOLDER,userid+"-biogrid_annotation.csv"))
			biogrid_df.to_excel(os.path.join(CSV_FOLDER,userid+"-biogrid_annotation.xlsx"), sheet_name='Biogrid-annotation')
			result.append({"displayName":"BIOGRID" ,"fileName": userid+"-biogrid_annotation.csv"})
	return result

if __name__ == "__main__":
	to_csv(sys.argv[1])
