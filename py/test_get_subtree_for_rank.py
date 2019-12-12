import get_subtree_for_rank
# you cand ownload the taxonomy file from https://tree.opentreeoflife.org/about/taxonomy-version/ott3.1
# OR download the taxonomy file with this function, it will be saved as ott3.1.tsv in your current dir
get_subtree_for_rank.download_taxonomy_file()

all_gen = get_subtree_for_rank.get_ott_ids_for_rank(rank = "genus", taxonomy_file = 'taxonomy_clean.tsv', clean = False)
all_subgen = get_subtree_for_rank.get_ott_ids_for_rank(rank = "subgenus", taxonomy_file = 'taxonomy_clean.tsv', clean = False)

get_subtree_for_rank.get_tree(rank = "family", taxonomy_file = 'replace by path and name of tsv file')
get_subtree_for_rank.get_tree(rank = "family", taxonomy_file = '../ott3.1.tsv')
tree = get_subtree_for_rank.get_tree(rank = "family", taxonomy_file = 'taxonomy_clean.tsv', clean = False)
tree = get_subtree_for_rank.get_tree(rank = "genus", taxonomy_file = 'taxonomy_clean.tsv', clean = False, label_format = "name_and_id")


get_subtree_for_rank.get_tree_from_synth(ott_ids = ['544595'])
group = ['544595']
subtree = get_subtree_for_rank.get_tree_from_synth(ott_ids = group, label_format="name_and_id", citation="cites.txt", tree_type = "subtree")

from Bio import Phylo
Phylo.read(tree)
import pylab # not working in physcraper virtual environment, need to import it
Phylo.draw(tree)
