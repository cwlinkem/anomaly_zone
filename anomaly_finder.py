#! /usr/bin/ python

import math
from math import exp
import dendropy


def anomaly_calc(nodes):
	lambda1=nodes[0]
	lambda2=nodes[1]
	if lambda1 is not None:
		Z1=math.log(2.0/3+((3*exp(2*lambda1)-2)/(18*(exp(3*lambda1)-exp(2*lambda1))))) #function a(x) from Degnan and Rosenberg (2006)
		if lambda2 <= Z1:
			return 1
		else:
			return 0
	else:
		return 0

def split_mapper(tree, taxon_set):
	split_dict={}
	for node in tree.postorder_node_iter():
		if node.parent_node is None:
			node.value = 1.0
		else:
			taxlist=[]
			tax = taxon_set.split_taxa_list(node.edge.split_bitmask)
			for i in tax:
				taxlist.append(i.label)
			split_dict[dendropy.treesplit.split_as_string(node.edge.split_bitmask, width=len(taxon_set))] = taxlist
	return split_dict

def split_freq(tree,taxon_set, split_occ_dict):
	for node in tree.postorder_node_iter():
		if node.parent_node is None:
			node.value = 1.0
		else:
			split=dendropy.treesplit.split_as_string(node.edge.split_bitmask, width=len(taxon_set))
			if split not in split_occ_dict.keys():
				split_occ_dict[split] = [1]
			else:
				split_occ_dict[split].append(1)
	return split_occ_dict
			
def get_nodes(tree, master_dict):
	pair_list=[]
	for node in tree.postorder_node_iter():
		if node.parent_node is None:
			node.value = 1.0
		else:
			if node.taxon is None:
				node_pair = (node.parent_node.edge_length, node.edge_length)
				anomalous=anomaly_calc(node_pair)
				edge_pair = (node.parent_node.edge.split_bitmask, node.edge.split_bitmask)
				if edge_pair not in master_dict.keys():
					master_dict[edge_pair] = [anomalous]
				else:
					master_dict[edge_pair].append(anomalous)	
	return master_dict		


taxon_set = dendropy.TaxonSet()
trees = dendropy.TreeList.get_from_path('mpest_boots.tre', "nexus", taxon_set = taxon_set)
mrc = dendropy.Tree.get_from_path('mpest_eMRC.tre', "nexus", taxon_set = taxon_set)

master_dict={}
mrc_dict={}
split_occ_dict={}
for sp_tree in trees:
	dendropy.treesplit.encode_splits(sp_tree)
	master_dict = get_nodes(sp_tree, master_dict)
	split_occ_dict=split_freq(sp_tree, taxon_set, split_occ_dict)	

dendropy.treesplit.encode_splits(mrc)
split_dict=split_mapper(mrc, taxon_set) #match a list of taxa to the split pattern from the bitmask for the whole tree
mrc_dict=get_nodes(mrc, mrc_dict)


for k,v in master_dict.iteritems():
	taxlab1=[]
	taxlab2=[]
	if k in mrc_dict.keys():
		tax1 = dendropy.treesplit.split_as_string(k[0], width = len(taxon_set))
		tax11 = taxon_set.split_taxa_list(k[0])
		tax2 = dendropy.treesplit.split_as_string(k[1], width = len(taxon_set))
		tax22 = taxon_set.split_taxa_list(k[1])
		for i in tax11:
			taxlab1.append(i.label)
		for j in tax22:
			taxlab2.append(j.label)
			
		print tax1, tax2, len(v), sum(v)
		
