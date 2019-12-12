import get_subtree_for_rank
# you cand ownload the taxonomy file from https://tree.opentreeoflife.org/about/taxonomy-version/ott3.1
# OR download the taxonomy file with this function, it will be saved as ott3.1.tsv in your current dir
get_subtree_for_rank.download_taxonomy_file()

all_gen = get_ott_ids_for_rank(rank = "genus", taxonomy_file = 'taxonomy_clean.tsv', clean = False)

get_subtree_for_rank.get_tree(rank = "family", taxonomy_file = 'replace by path and name of tsv file')
get_subtree_for_rank.get_tree(rank = "family", taxonomy_file = '../ott3.1.tsv')
tree = get_subtree_for_rank.get_tree(rank = "family", taxonomy_file = 'taxonomy_clean.tsv', clean = False)

from Bio import Phylo
Phylo.draw(tree)
