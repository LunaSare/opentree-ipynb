# Get an opentree synthetic tree for a rank X, from a clade Y,
# distributed in geographic region Z.

# arguments
# ott_ids: Numeric. If given only one, it will get a subtree for the corresponding node
         # If given 2 or more, it will give back an induced subtree
# rank: Character string. Indicates the taxonomic rank to search for.
# clade: Character string or number (?). Indicates a phylogenetic context to search for.
# range: Character string. Indicates the geogaphic region to search for.
# ott_version: Character string
_DEBUG = 1
def debug(msg):
    """short debugging command
    """
    if _DEBUG == 1:
        print(msg)

import os
from physcraper import opentree_helpers
from physcraper.treetaxon import TreeTax

from dendropy import Tree, DnaCharacterMatrix, DataSet, datamodel

# import get_citations_from_json # should not need this anymore because
# it is already incorporated in module opentree_helpers

import gzip
import shutil

def download_taxonomy_file(version = '3.1'):
    # download ott taxonomy
    os.system('wget "http://files.opentreeoflife.org/ott/ott' + version + '/ott' + version + '.tgz"')
    # unzip it
    with gzip.open('ott' + version + '.tgz', 'rb') as fin:
        with open('ott' + version + '.tsv', 'wb') as fout:
            shutil.copyfileobj(fin, fout)


def get_ott_ids_for_rank(rank, taxonomy_file = 'ott3.1/taxonomy.tsv'):
    sys.stdout.write("Gathering ott ids from {}\n".format(taxonomy_file))
    all_ranks = ['species', 'genus', 'family', 'order', "class"]
    # clean taxonomy file
    os.system('grep -a "' + rank + '" ' + taxonomy_file + ' | egrep -v "Incertae" | egrep -v "no rank" | egrep -v "major_rank_conflict" | egrep -v "uncultured" | egrep -v "barren" | egrep -v "extinct" | egrep -v "incertae" | egrep -v "unplaced" | egrep -v "hidden" | egrep -v "inconsistent"  | egrep -v "synonym" | egrep -v "in ' + rank + '" | egrep -v "species" | egrep -v "genus" | egrep -v "super' + rank + '" | egrep -v "sub' + rank + '" > taxonomy_red.tsv')
    # extract ott ids from taxonomy reduced file
    taxonomy_tsv = 'taxonomy_red.tsv'
    fi = open(taxonomy_tsv).readlines()
    ott_ids = []
    for lin in fi:
        lii = lin.split('\t')
        ott_ids.append(lii[0])
    ott_ids = list(ott_ids)
    return ott_ids

def get_tree(rank, taxonomy_file = 'ott3.1/taxonomy.tsv'):
    ott_ids = get_ott_ids_for_rank(rank, taxonomy_file)
    # tre = get_citations_from_json.get_tree_from_synth(ott_ids, citation = 'citations_' + rank + '.txt')
    tre = opentree_helpers.get_tree_from_synth(ott_ids, citation = "citations_" + rank + ".txt")
    return tre
