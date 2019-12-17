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
import re
import sys
from physcraper import opentree_helpers
from physcraper.treetaxon import TreeTax
from dendropy import Tree, DnaCharacterMatrix, DataSet, datamodel
import gzip
import shutil
from Bio import Phylo
import get_subtree_for_rank_DEV

def download_taxonomy_file(version = '3.1'):
    # download ott taxonomy
    os.system('wget "http://files.opentreeoflife.org/ott/ott' + version + '/ott' + version + '.tgz"')
    # unzip it
    with gzip.open('ott' + version + '.tgz', 'rb') as fin:
        with open('ott' + version + '.tsv', 'wb') as fout:
            shutil.copyfileobj(fin, fout)

# the function cleans up the word 'species' and the flag 'no rank - terminal', which is not associated to higher taxonomic ranks
def clean_taxonomy_file(taxonomy_file = 'ott3.1/taxonomy.tsv'):
    sys.stdout.write('Cleaning {} file... '.format(taxonomy_file))
    # clean taxonomy file, writes cleaned file to taxonomy_clean.tsv
    os.system('grep -a -v "major_rank_conflict" ' + taxonomy_file + ' | egrep -a -v "species" | egrep -a -v "varietas" | egrep -a -v "no rank" | egrep -a -v "Incertae" | egrep -a -v "incertae" | egrep -a -v "uncultured" | egrep -a -v "barren" | egrep -a -v "extinct" | egrep -a -v "unplaced" | egrep -a -v "hidden" | egrep -a -v "inconsistent" | egrep -a -v "synonym" > taxonomy_clean.tsv')
    sys.stdout.write("Done.\n")

#accepts several ranks at the same time
def get_ott_ids_for_rank(rank, taxonomy_file = 'ott3.1/taxonomy.tsv', clean = True, write_file = 'rank_ott_ids.txt'):
    taxonomy_tsv = taxonomy_file
    # clean taxonomy file
    if clean:
        clean_taxonomy_file(taxonomy_file)
        taxonomy_tsv = 'taxonomy_clean.tsv'
    # extract ott ids from taxonomy reduced file
    sys.stdout.write('Gathering ott ids from {}...\n'.format(rank))
    fi = open(taxonomy_tsv).readlines()
    ott_ids = []
    for lin in fi:
        lii = re.split('\t*\|\t*', lin)
        if re.match('[0-9]', lii[0]):
            if len(lii) > 2:
                if lii[3] in rank:
                    ott_ids.append(lii[0])
    if isinstance(write_file, str):
        with open(write_file, 'w') as f:
            for item in ott_ids:
                print >> f, item
    ott_ids = list(ott_ids)
    return ott_ids

def get_ott_ids_for_group(group_ott_id, write_file = 'children_ott_ids.txt'):
    sys.stdout.write('Gathering ott ids from group with ott id {}...\n'.format(group_ott_id[0]))
    debug(group_ott_id)
    subtree = get_subtree_for_rank_DEV.get_tree_from_synth(ott_ids = group_ott_id, label_format='name_and_id', citation='cites.txt', tree_type = 'subtree')
    newick = str(subtree)
    # get ott ids from newick
    ott_ids0 = re.findall(r'ott\d+', newick)
    ott_ids = []
    for item in ott_ids0:
        ott_ids.append(re.findall(r'\d+', item))
    ott_ids = [item for sublist in ott_ids for item in sublist] # flattens the list
    if isinstance(write_file, str):
        with open(write_file, 'w') as f:
            for item in ott_ids:
                print >> f, item
    return ott_ids



def get_ott_ids_X(group_ott_id = None, group_ott_ids_file = None, rank = "family", taxonomy_file = 'ott3.1/taxonomy.tsv', clean = True):
    taxonomy_tsv = taxonomy_file
    # clean taxonomy file
    if clean:
        clean_taxonomy_file(taxonomy_file)
        taxonomy_tsv = 'taxonomy_clean.tsv'
    # extract ott ids from taxonomy reduced file
    if isinstance(group_ott_ids_file, str):
        sys.stdout.write('Getting ott ids from file {}...\n'.format(group_ott_ids_file))
        children_ott_ids = [line.rstrip('\n') for line in open(group_ott_ids_file)]
    else:
        children_ott_ids = get_ott_ids_for_group(group_ott_id)
    sys.stdout.write('Gathering ott ids from {} in group...\n'.format(rank))
    fi = open(taxonomy_tsv).readlines()
    ott_ids = []
    # debug(len(children_ott_ids))
    for line in fi:
        lii = re.split('\t*\|\t*', line)
        if re.match('[0-9]', lii[0]): # skips the headers and other weird lines
            if len(lii) > 2:
                if lii[0] in children_ott_ids:
                    # debug(lii[3])
                    # if re.match(rank, lii[3]):
                    if lii[3] in rank:
                        ott_ids.append(lii[0])
    ott_ids = list(ott_ids)
    return ott_ids


def get_tree(group_ott_id = None, group_ott_ids_file = None, rank = 'family', taxonomy_file = 'ott3.1/taxonomy.tsv', clean = True, label_format='name'):
    # assert group_ott_id is not None and group_ott_ids_file is not None and rank is not None
    if group_ott_id is not None and rank is None:
        ott_ids = get_ott_ids_for_group(group_ott_id)
    if group_ott_id is None and group_ott_ids_file is None and rank is not None:
        ott_ids = get_ott_ids_for_rank(rank, taxonomy_file, clean = clean)
    if(group_ott_id is not None or group_ott_ids_file is not None and rank is not None):
        ott_ids = get_ott_ids_X(group_ott_id, group_ott_ids_file, rank, taxonomy_file, clean)
    # tre = get_citations_from_json.get_tree_from_synth(ott_ids, citation = 'citations_' + rank + '.txt')
    tre = opentree_helpers.get_tree_from_synth(ott_ids, label_format = label_format, citation = 'citations_' + rank + '.txt')
    return tre
