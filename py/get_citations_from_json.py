
_DEBUG = 0
def debug(msg):
    """short debugging command
    """
    if _DEBUG == 1:
        print(msg)

import sys
from physcraper.helpers import cd, standardize_label, to_string
import json
import requests

def get_citations_from_json(study_id, citations_file):
    assert isinstance(citations_file, str)
    f = open(citations_file,"w+")
    sys.stdout.write("Gathering citations ...")
    for study in study_id.json()['supporting_studies']:
        study = study.split('@')[0]
        index_url = 'https://api.opentreeoflife.org/v3/studies/find_studies'
        payload = json.dumps({"property":"ot:studyId","value":study,"verbose":"true"})
        res_cites = requests.post(index_url, data=payload, headers={'content-type':'application/json'})
        new_cite = res_cites.json()['matched_studies']
        debug(new_cite)
        f.write(to_string(new_cite[0]['ot:studyPublicationReference']) + '\n' + new_cite[0]['ot:studyPublication'] + '\n')
    f.close()
    sys.stdout.write("Citations printed to {}\n".format(citations_file))

# this "new" get_tree_from_synth function returns citations from supporting studies
def get_tree_from_synth(ott_ids, label_format="name", citation="cites.txt"):
    assert label_format in ['id', 'name', 'name_and_id']
    url = 'https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree'
    headers = {'content-type':'application/json'}
    payload = json.dumps(dict(ott_ids=ott_ids, label_format = label_format))
    res = requests.post(url, data=payload, headers=headers)
    if res.status_code == 200:
        pass
    else:
        sys.stderr.write("error getting synth tree, {}, {}, {}\n".format(res.status_code, res.reason, res.json()['message']))
        return None
    get_citations_from_json(res, citation) # returns file with citations
    tre = Tree.get(data=res.json()['newick'],
                   schema="newick",
                   suppress_internal_node_taxa=True)
    tre.suppress_unifurcations()
    return tre
