import get_subtree_for_rank
# you cand ownload the taxonomy file from https://tree.opentreeoflife.org/about/taxonomy-version/ott3.1
# OR download the taxonomy file with this function, it will be saved as ott3.1.tsv in your current dir
get_subtree_for_rank.download_taxonomy_file()

get_subtree_for_rank.clean_taxonomy_file(taxonomy_file = 'ott3.1/taxonomy.tsv')

all_gen = get_subtree_for_rank.get_ott_ids_for_rank(rank = "genus", taxonomy_file = 'taxonomy_clean.tsv', clean = False)
all_subgen = get_subtree_for_rank.get_ott_ids_for_rank(rank = "subgenus", taxonomy_file = 'taxonomy_clean.tsv', clean = False)

get_subtree_for_rank.get_tree(rank = "family", taxonomy_file = 'replace by path and name of tsv file')
get_subtree_for_rank.get_tree(rank = "family", taxonomy_file = '../ott3.1.tsv')
tree = get_subtree_for_rank.get_tree(rank = "family", taxonomy_file = 'taxonomy_clean.tsv', clean = False)
tree = get_subtree_for_rank.get_tree(rank = "genus", taxonomy_file = 'taxonomy_clean.tsv', clean = False, label_format = "name_and_id")


get_subtree_for_rank.get_tree_from_synth(ott_ids = ['544595'])
group = ['544595'] # amphibia
subtree = get_subtree_for_rank.get_tree_from_synth(ott_ids = group, label_format="name_and_id", citation="cites.txt", tree_type = "subtree")

get_subtree_for_rank.get_ott_ids_for_rank_and_group(rank = "family", group_ott_id = ['544595'], taxonomy_file = 'taxonomy_clean.tsv', clean = False,

children_ott_ids = get_subtree_for_rank.get_ott_ids_for_group(group_ott_id = ['544595'], write_file = "children_amphibia.txt")
len(children_ott_ids) # 10318
rank_ott_id = get_subtree_for_rank.get_ott_ids_for_rank(rank = "family", taxonomy_file = 'taxonomy_clean.tsv', clean = False)
len(rank_ott_id) # 9616
amph_fam = get_subtree_for_rank.get_ott_ids_group_and_rank(group_ott_id = ['544595'], rank = "family", taxonomy_file = 'taxonomy_clean.tsv', clean = False)
len(amph_fam) # 47
amph_gen = get_subtree_for_rank.get_ott_ids_group_and_rank(group_ott_id = ['544595'], rank = "genus", taxonomy_file = 'taxonomy_clean.tsv', clean = False)
len(amph_gen) # 386
amph_fam = get_subtree_for_rank.get_ott_ids_group_and_rank(group_ott_id = ['544595'], group_ott_ids_file = 'children_amphibia.txt', rank = "family", taxonomy_file = 'taxonomy_clean.tsv', clean = False)

amph_fam = get_subtree_for_rank.get_ott_ids_X(group_ott_ids_file = 'children_amphibia.txt', rank = ["superfamily", "family", "subfamily"], taxonomy_file = 'taxonomy_clean.tsv', clean = False)
amph_fam = get_subtree_for_rank.get_ott_ids_X(group_ott_id = ['544595'], rank = ["superfamily", "family", "subfamily"], taxonomy_file = 'taxonomy_clean.tsv', clean = False)

tree = get_subtree_for_rank.get_tree(group_ott_id = ['544595'], rank = "family", taxonomy_file = 'taxonomy_clean.tsv', clean = False)
tree = get_subtree_for_rank.get_tree(group_ott_ids_file = 'children_amphibia.txt', rank = "family", taxonomy_file = 'taxonomy_clean.tsv', clean = False)


from Bio import Phylo
Phylo.read(tree)
import pylab # not working in physcraper virtual environment, need to import it
Phylo.draw(tree)
