#BEGIN_HEADER
# The header block is where all import statements should live
import sys
import traceback
import uuid
import json
import requests
import random
from pprint import pprint, pformat
from biokbase.workspace.client import Workspace as workspaceService
#END_HEADER


class HomologySearch:
    '''
    Module Name:
    HomologySearch

    Module Description:
    A KBase module: HomologySearch
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    # Class variables and functions can be defined in this block
    workspaceURL = None
    def formatHspList(self, list):
        name_map = {"Hsp_align-len": "align_len", "Hsp_bit-score": "bit_score", "Hsp_evalue": "evalue",
                    'Hsp_hit-frame': "hit_frame", "Hsp_hit-from": "hit_from", "Hsp_hit-to": "hit_to",
                    "Hsp_hit-strand": "hit_strand",
                    "Hsp_hseq": "hseq", "Hsp_identity": "identity", "Hsp_midline": "midline", "Hsp_num": "num",
                    "Hsp_positive": "positive", "Hsp_qseq": "qseq",
                    "Hsp_query-frame": "query_frame", "Hsp_query-from": "query_from", "Hsp_query-to": "query_to",
                    "Hsp_query-strand": "query_strand", "Hsp_score": "score", "Hsp_gaps": "gaps"}
        returnVal = []
        for item in list:
            returnVal.append(dict(map(lambda (k, v): (k, str(item.get(v, ''))), name_map.iteritems())))
        return returnVal
    def formatHitList(self, list, metadata):
        name_map = {"Hit_accession": "", "Hit_def": "", "Hit_id": "", "Hit_len": "len",
        "Hit_num": "num", "Hit_hsps": ""}
        returnVal = []
        for item in list:
            newItem = {}
            for name, val in name_map.items():
                if name is "Hit_accession":
                    newItem[name] = item['description'][0]['accession']
                elif name is "Hit_id":
                    newItem[name] = item['description'][0]['id']
                elif name is "Hit_def":
                    newItem[name] = metadata[item['description'][0]['id']].get('function', '')
                elif name is "Hit_hsps":
                    newItem[name] = {'Hsp': self.formatHspList(item['hsps'])}
                elif val in item:
                    newItem[name] = str(item[val])
            returnVal.append(newItem)
        return returnVal
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.homologyServiceURL = config['homologyservice-url']
        #END_CONSTRUCTOR
        pass

    def log(self, target, message):
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    def run_blast_fasta(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN run_blast_fasta
        console = []
        self.log(console, 'Starting run_blast_fasta with params=')
        self.log(console, pformat(params))

        # Step 1 - parse parameters
        if 'sequence' not in params:
            raise ValueError('Parameter sequence is not set in input arguments')
        if 'database' not in params:
            raise ValueError('Parameter database is not set in input arguments')
        if 'program' not in params:
            raise ValueError('Parameter program is not set in input arguments')
        if 'evalue_cutoff' not in params:
            raise ValueError('Parameter evalue_cutoff is not set in input arguments')
        if 'max_hit' not in params:
            raise ValueError('Parameter max_hit is not set in input arguments')

        # Step 2 - query
        if params['database'] is not "":
            req_method = "HomologyService.blast_fasta_to_database"
            req_params = [params['sequence'], params['program'], params['database'], params['evalue_cutoff'], params['max_hit'], 0]
        else:
            req_method = "HomologyService.blast_fasta_to_genomes"
            req_params = [params['sequence'], params['program'], params['genome_ids'], params['search_type'], params['evalue_cutoff'], params['max_hit'], 0]

        rpc = {"version": "1.1", "params": req_params, "method": req_method, "id": str(random.random())[2:]}

        self.log(console, 'sending query')
        self.log(console, pformat(rpc))

        req = requests.post(self.homologyServiceURL, json=rpc)

        # Step 3 - wrap results
        jsonQueryReturn = req.json()
        report = jsonQueryReturn['result'][0][0]['report']
        query = report['results']['search']
        hitsList = report['results']['search']['hits']
        metadata = jsonQueryReturn['result'][1]
        # identical = jsonQueryReturn['result'][2];

        # hits = self.formatHitList(hitsList)
        # print hits
        # hspsList = report["results"]["search"]["hits"][0]["hsps"]
        # hsps = self.formatHspList(hspsList)
        # print hsps
        self.log(console, 'hits:')
        self.log(console, pformat(hitsList))
        self.log(console, 'metadata:')
        self.log(console, pformat(metadata))

        returnVal = {
            "BlastOutput_db": "",
            "BlastOutput_program": report['program'],
            "BlastOutput_query-ID": query['query_id'],
            "BlastOutput_query-def": query.get('query_title', ''),
            "BlastOutput_query-len": str(query['query_len']),
            "BlastOutput_reference": report['reference'],
            "BlastOutput_version": report['version'],
            "BlastOutput_param": {
                "Parameters": {
                    "Parameters_expect": str(report['params'].get('expect', '')),
                    "Parameters_filter": report['params'].get('filter', ''),
                    "Parameters_gap-extend": str(report['params'].get('gap_extend', '')),
                    "Parameters_gap-open": str(report['params'].get('gap_open', '')),
                    "Parameters_matrix": report['params'].get('matrix', ''),
                    "Parameters_sc-match": "",
                    "Parameters_sc-mismatch": ""
                }
            },
            "BlastOutput_iterations": {
                "Iteration": [
                    {
                        "Iteration_hits": {
                            "Hit": self.formatHitList(hitsList, metadata)
                        },
                        "Iteration_iter-num": "",
                        "Iteration_query-ID": query['query_id'],
                        "Iteration_query-def": query.get('query_title', ''),
                        "Iteration_query-len": str(query['query_len']),
                        "Iteration_stat": {
                            "Statistics": {
                                "Statistics_db-len": str(query['stat']['db_len']),
                                "Statistics_db-num": str(query['stat']['db_num']),
                                "Statistics_eff-space": str(query['stat']['eff_space']),
                                "Statistics_entropy": str(query['stat']['entropy']),
                                "Statistics_hsp-len": str(query['stat']['hsp_len']),
                                "Statistics_kappa": str(query['stat']['kappa']),
                                "Statistics_lambda": str(query['stat']['lambda'])
                            }
                        }
                    }
                ]
            }
        }
        # print returnVal
        # Step 4 - save in workspace
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']

        token = ctx['token']
        ws = workspaceService(self.workspaceURL, token=token)
        res = ws.save_objects(
            {
                "workspace": params["workspace_name"],
                "objects": [{
                    "type": "GenomeUtil.BlastOutput",
                    "data": returnVal,
                    "provenance": provenance,
                    "name": params["output_name"]
                }]
            }
        )
        self.log(console, 'saved to workspace')
        returnVal = {'blast_output_name': params["output_name"], 'workspaceName': params['workspace_name']}
        self.log(console, pformat(returnVal))
        #END run_blast_fasta

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method run_blast_fasta return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
