import sys
import re
import os


def build_summary(go_features=None,interactions=None,pathways=None,main_dict=None, cross_dict=None, main_genes=None):
    if go_features:
        gene = go_features["gene"]
        go_namespace = go_features["ns"]
        goterms = go_features["go"]
        if gene in main_genes: 
            try:
               main_dict[gene][0][go_namespace]=len(set(goterms))
            except KeyError:
               main_dict[gene] =[{}]
               main_dict[gene][0][go_namespace]=len(set(goterms))
        else:
            try:
               cross_dict[gene][0][go_namespace]=len(set(goterms))
            except KeyError:
               cross_dict[gene] =[{}]
               cross_dict[gene][0][go_namespace]=len(set(goterms))
    if interactions:
        for node in interactions.keys():
            if node in main_genes:
                try:
                    main_dict[node][0]["Biogrid Interactors"]= len(set(interactions[node]))
                    main_dict[node][0]["Biogrid Interactions"]= len(interactions[node])
                except KeyError:
                    main_dict[node] =[{}]
                    main_dict[node][0]["Biogrid Interactors"]= len(set(interactions[node]))
                    main_dict[node][0]["Biogrid Interactions"]= len(interactions[node])
            else:
                try:
                    cross_dict[node][0]["Biogrid Interactors"]= len(set(interactions[node]))
                    cross_dict[node][0]["Biogrid Interactions"]= len(interactions[node])
                except KeyError:
                    cross_dict[node] =[{}]
                    cross_dict[node][0]["Biogrid Interactors"]= len(set(interactions[node]))
                    cross_dict[node][0]["Biogrid Interactions"]= len(interactions[node])
    if pathways:
        pw_features_lst = [pathways['gene'],pathways["small_molecule"], pathways["protein"]]
        pathway = pathways["pathway"]
        flattened_list = [y for x in pw_features_lst for y in x]      
        for node in flattened_list:
            if node in main_genes:
                try:
                    main_dict[node][0]["Pathways"]= len(set(pathway))
                except KeyError:
                    main_dict[node] =[{}]
                    main_dict[node][0]["Pathways"]= len(set(pathway))
            else:
                try:
                    cross_dict[node][0]["Pathways"]= len(set(pathway))
                except KeyError:
                    cross_dict[node] =[{}]
                    cross_dict[node][0]["Pathways"]= len(set(pathway))

    return main_dict, cross_dict
