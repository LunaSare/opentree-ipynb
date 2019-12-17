ott_ids = ['1098176',
 '578817',
 '441950',
 '767327',
 '538123',
 '17900',
 '135338',
 '328373',
 '1039354',
 '1039356']

tre = opentree_helpers.get_tree_from_synth(ott_ids, citation = "cites_test1.txt")


children_ott_ids = ['1098176',
 '578817',
 '441950',
 '767327',
 '538123',
 '17900',
 '135338',
 '328373',
 '1039354',
 '1039356']

rank_ott_ids = ['1098176',
 '441950',
 '538123',
 '135338',
 '1039354']

ii = [rank_ott_ids.index(i) for i in children_ott_ids]

res = [key for key, val in enumerate(children_ott_ids)
    if val in set(rank_ott_ids)]

[children_ott_ids[i] for i in res]

children_ott_ids = get_subtree_for_rank.get_ott_ids_for_group(group_ott_id)
[item for sublist in children_ott_ids for item in sublist]
rank_ott_ids = get_subtree_for_rank.get_ott_ids_for_rank(rank, taxonomy_file, clean)

indices = [rank_ott_ids.index(i) for i in children_ott_ids]
indices = [children_ott_ids.index(i) for i in rank_ott_ids]
ott_ids = [children_ott_ids[i] for i in indices]
