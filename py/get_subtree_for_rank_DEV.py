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

# the function cleans up the word 'species' and the flag 'no rank - terminal', which is not associated to higher taxonomic ranks
def clean_taxonomy_file(taxonomy_file = 'ott3.1/taxonomy.tsv'):
    sys.stdout.write('Cleaning {} file... '.format(taxonomy_file))
    # clean taxonomy file, writes cleaned file to taxonomy_clean.tsv
    os.system('grep -a -v "major_rank_conflict" ' + taxonomy_file + ' | egrep -a -v "species" | egrep -a -v "varietas" | egrep -a -v "no rank" | egrep -a -v "Incertae" | egrep -a -v "incertae" | egrep -a -v "uncultured" | egrep -a -v "barren" | egrep -a -v "extinct" | egrep -a -v "unplaced" | egrep -a -v "hidden" | egrep -a -v "inconsistent" | egrep -a -v "synonym" > taxonomy_clean.tsv')
    sys.stdout.write("Done.\n")


def get_ott_ids_for_rank(rank, taxonomy_file = 'ott3.1/taxonomy.tsv', clean = True):
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
        # lii = re.split('\t*', lin)
        lii = re.split('\t*\|\t*', lin)
        if re.match('[0-9]', lii[0]):
            if len(lii) > 2:
                if re.match(rank, lii[3]):
                    ott_ids.append(lii[0])
                    sys.stdout.write(".")
    ott_ids = list(ott_ids)
    return ott_ids


def get_tree(rank, taxonomy_file = 'ott3.1/taxonomy.tsv', clean = True, label_format='name'):
    ott_ids = get_ott_ids_for_rank(rank, taxonomy_file, clean = clean)
    # tre = get_citations_from_json.get_tree_from_synth(ott_ids, citation = 'citations_' + rank + '.txt')
    tre = opentree_helpers.get_tree_from_synth(ott_ids, label_format = label_format, citation = 'citations_' + rank + '.txt')
    return tre


import json
import requests
def get_tree_from_synth(ott_ids, label_format='name', citation='cites.txt', tree_type = 'subtree'):
    assert label_format in ['id', 'name', 'name_and_id']
    assert tree_type in ['subtree', 'induced_subtree']
    assert isinstance(ott_ids, list) # if it is a string it is very hard to know how to delimit ott ids
    if 'subtree' in tree_type:
            ott_ids = ott_ids[0] # take only first element of ott_ids if it is a list of multiple items
    url = 'https://api.opentreeoflife.org/v3/tree_of_life/' + tree_type
    headers = {'content-type':'application/json'}
    pass_number = 0
    while pass_number <= 1:
        if 'subtree' in tree_type:
            payload = json.dumps(dict(ott_id=ott_ids, label_format = label_format))
        else:
            payload = json.dumps(dict(ott_ids=ott_ids, label_format = label_format))
        res = requests.post(url, data=payload, headers=headers)
        if res.status_code == 200:
            pass_number += 2
            break
        else:
            pass_number += 1
            if 'unknown' in res.json():
                bad_ids = res.json()['unknown'].keys()
                ott_ids = set(ott_ids)
                for bad_ott_id in bad_ids:
                    num = bad_ott_id.strip('ott')
                    ott_ids.remove(num)
                ott_ids = list(ott_ids)
        if pass_number == 2:
            sys.stderr.write('Error getting synth tree, {}, {}, {}, (full error ott ids hidden)\n'.format(res.status_code, res.reason, res.json().get('message'), res.json()))
            return None
    synth_json = res.json()
    tre = Tree.get(data=synth_json['newick'],
                   schema='newick',
                   suppress_internal_node_taxa=True)
    assert 'supporting_studies' in synth_json.keys(), synth_json.keys()
    opentree_helpers.get_citations_from_json(synth_json, citation)
    tre.suppress_unifurcations()
    return tre

from Bio import Phylo
def get_ott_ids_for_group(group_ott_id, write_file = 'children_ott_ids.txt'):
    sys.stdout.write('Gathering ott ids from group with ott id {}...\n'.format(group_ott_id[0]))
    debug(group_ott_id)
    subtree = get_tree_from_synth(ott_ids = group_ott_id, label_format='name_and_id', citation='cites.txt', tree_type = 'subtree')
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

def get_ott_ids_group_and_rank(group_ott_id = None, group_ott_ids_file = None, rank = None, rank_ott_ids_file = None, taxonomy_file = 'ott3.1/taxonomy.tsv', clean = True):
    # clean taxonomy file
    taxonomy_tsv = taxonomy_file
    if clean:
        clean_taxonomy_file(taxonomy_file)
        taxonomy_tsv = 'taxonomy_clean.tsv'
    # get group ott ids
    if isinstance(group_ott_ids_file, str):
        sys.stdout.write('Getting ott ids from file {}...\n'.format(group_ott_ids_file))
        children_ott_ids = [line.rstrip('\n') for line in open(group_ott_ids_file)]
    else:
        children_ott_ids = get_ott_ids_for_group(group_ott_id)
    # get rank ott ids
    rank_ott_ids = get_ott_ids_for_rank(rank, taxonomy_tsv, clean)
    # get rank ott ids that are in children ott ids:
    ott_ids = []
    for item in children_ott_ids:
        if item in rank_ott_ids:
            ott_ids.append(item)
    # alternatively: subset taxonomy file by group_ott_ids and write to taxonomy_group.tsv:
    # tried it, but it is REALLYYYY slow (it seems following line is not working as expected)
    # os.system('grep -Fwf children_ott_ids.txt ' + taxonomy_tsv + ' > taxonomy_clean_group.tsv')
    # not working yet...
    # os.system('xargs -I {} grep "^{}" ' + taxonomy_tsv + ' < children_ott_ids.txt > taxonomy_clean_group3.tsv')
    # fi = open(taxonomy_tsv).readlines()
    # ott_ids = []
    # for lin in fi:
    #     lii = re.split('\t*\|\t*', lin)
    #     if re.match('[0-9]', lii[0]):
    #         if len(lii) > 2:
    #             if re.match(lii[0], ): # filter children ott ids
    #                 if re.match(rank, lii[3]): # filter rank
    #                     ott_ids.append(lii[0])
    #                     sys.stdout.write(".")
    # ott_ids = list(ott_ids)
    return ott_ids

#accepts several ranks at the same time
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
    sys.stdout.write('Gathering ott ids from {}...\n'.format(rank))
    fi = open(taxonomy_tsv).readlines()
    ott_ids = []
    # debug(len(children_ott_ids))
    for line in fi:
        lii = re.split('\t*\|\t*', line)
        if re.match('[0-9]', lii[0]): # skips the headers and other weird lines
            if len(lii) > 2:
                if lii[0] in children_ott_ids:
                    debug(lii[3])
                    # if re.match(rank, lii[3]):
                    if lii[3] in rank:
                        ott_ids.append(lii[0])
    ott_ids = list(ott_ids)
    return ott_ids
