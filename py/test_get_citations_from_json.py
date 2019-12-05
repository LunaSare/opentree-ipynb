### EXAMPLE OF SUBTREE TO TEST FUNCTION get_citations_from_json

import os
from physcraper import opentree_helpers
from physcraper.treetaxon import TreeTax
json_file = "../../../OpenTree_SSB2020/tutorial/main.json"
assert os.path.isfile(json_file) #check the file exists and the path is correct
otu_dict = opentree_helpers.bulk_tnrs_load(json_file)
ott_ids =set()
for otu in otu_dict:
   ott_ids.add(otu_dict[otu].get("^ot:ottId"))
#turn it back into a list
ott_ids = list(ott_ids)
import json
import requests
label_format="name"
url = 'https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree'
payload = json.dumps(dict(ott_ids=ott_ids, label_format = label_format))
res = requests.post(url, data=payload, headers={'content-type':'application/json'})

import get_citations_from_json as jsoncites # should not need this anymore because
# it is already incorporated in module opentree_helpers

jsoncites.get_citations_from_json(res, "citations_test_2019.12.04.txt")

# opentree_helpers.get_citations_from_json(res, "citations_test_2019.12.04.txt") # throwing error:
#   File "test_get_citations_from_json.py", line 26, in <module>
#     opentree_helpers.get_citations_from_json(res, "citations_test_2019.12.04.txt")
#   File "/Users/luna/physcraper/physcraper/opentree_helpers.py", line 104, in get_citations_from_json
#     assert 'supporting_studies' in synth_response.keys(), synth_response.keys()
# AttributeError: 'Response' object has no attribute 'keys'
