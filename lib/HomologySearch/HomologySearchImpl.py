#BEGIN_HEADER
# The header block is where all import statments should live
import sys
import traceback
import uuid
import json
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
        name_map = {"Hsp_align-len": "align_len", "Hsp_bit-score": "bit_score",
        "Hsp_evalue": "evalue","Hsp_hit-from": "hit_from",
        "Hsp_hit-to": "hit_to","Hsp_hseq": "hseq", "Hsp_identity": "identity",
        "Hsp_midline": "midline","Hsp_num": "num", "Hsp_positive": "positive",
        "Hsp_qseq": "qseq","Hsp_query-from": "query_from", "Hsp_query-to": "query_to",
        "Hsp_score": "score"}
        returnVal = []
        for item in list:
            returnVal.append(dict(map(lambda (k,v): (k, str(item[v])), name_map.iteritems())))
        return returnVal
    def formatHitList(self, list):
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
                    newItem[name] = item['description'][0]['title']
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
        #END_CONSTRUCTOR
        pass

    def blast_fasta(self, ctx, params):
        # ctx is the context object
        # return variables are: output
        #BEGIN blast_fasta

        print('Starting blast_fasta')

        # Step 1 - parse parameters
        # if 'sequence' not in params:
        #     raise ValueError('Parameter sequence is not set in input arguments')
        # if 'search_type' not in params:
        #     raise ValueError('Parameter search type is not set in input arguments')
        # if 'program' not in params:
        #     raise ValueError('Parameter program is not set in input arguments')

        # Step 2 - query

        # Step 3 - wrap results
        # TODO: don't forget to add extra \ for \"
        queryReturn = '{"result":[[{"report":{"search_target":{"db":"/vol/patric3/kbase/blastdb/genomes/kb|g.0.faa"},"params":{"gap_open":11,"gap_extend":1,"cbs":2,"expect":1e-05,"filter":"F","matrix":"BLOSUM62"},"results":{"search":{"hits":[{"len":428,"description":[{"accession":"3656","title":"kb|g.0.peg.4288   Threonine synthase (EC 4.2.3.1)   [Escherichia coli K12]","id":"kb|g.0.peg.4288"}],"num":1,"hsps":[{"identity":428,"qseq":"MKLYNLKDHNEQVSFAQAVTQGLGKNQGLFFPHDLPEFSLTEIDEMLKLDFVTRSAKILSAFIGDEIPQEILEERVRAAFAFPAPVANVESDVGCLELFHGPTLAFKDFGGRFMAQMLTHIAGDKPVTILTATSGDTGAAVAHAFYGLPNVKVVILYPRGKISPLQEKLFCTLGGNIETVAIDGDFDACQALVKQAFDDEELKVALGLNSANSINISRLLAQICYYFEAVAQLPQETRNQLVVSVPSGNFGDLTAGLLAKSLGLPVKRFIAATNVNDTVPRFLHDGQWSPKATQATLSNAMDVSQPNNWPRVEELFRRKIWQLKELGYAAVDDETTQQTMRELKELGYTSEPHAAVAYRALRDQLNPGEYGLFLGTAHPAKFKESVEAILGETLDLPKELAERADLPLLSHNLPADFAALRKLMMNHQ","hit_to":428,"hit_from":1,"positive":428,"bit_score":882.093,"score":2278,"midline":"MKLYNLKDHNEQVSFAQAVTQGLGKNQGLFFPHDLPEFSLTEIDEMLKLDFVTRSAKILSAFIGDEIPQEILEERVRAAFAFPAPVANVESDVGCLELFHGPTLAFKDFGGRFMAQMLTHIAGDKPVTILTATSGDTGAAVAHAFYGLPNVKVVILYPRGKISPLQEKLFCTLGGNIETVAIDGDFDACQALVKQAFDDEELKVALGLNSANSINISRLLAQICYYFEAVAQLPQETRNQLVVSVPSGNFGDLTAGLLAKSLGLPVKRFIAATNVNDTVPRFLHDGQWSPKATQATLSNAMDVSQPNNWPRVEELFRRKIWQLKELGYAAVDDETTQQTMRELKELGYTSEPHAAVAYRALRDQLNPGEYGLFLGTAHPAKFKESVEAILGETLDLPKELAERADLPLLSHNLPADFAALRKLMMNHQ","gaps":0,"evalue":0,"query_to":428,"query_from":1,"hseq":"MKLYNLKDHNEQVSFAQAVTQGLGKNQGLFFPHDLPEFSLTEIDEMLKLDFVTRSAKILSAFIGDEIPQEILEERVRAAFAFPAPVANVESDVGCLELFHGPTLAFKDFGGRFMAQMLTHIAGDKPVTILTATSGDTGAAVAHAFYGLPNVKVVILYPRGKISPLQEKLFCTLGGNIETVAIDGDFDACQALVKQAFDDEELKVALGLNSANSINISRLLAQICYYFEAVAQLPQETRNQLVVSVPSGNFGDLTAGLLAKSLGLPVKRFIAATNVNDTVPRFLHDGQWSPKATQATLSNAMDVSQPNNWPRVEELFRRKIWQLKELGYAAVDDETTQQTMRELKELGYTSEPHAAVAYRALRDQLNPGEYGLFLGTAHPAKFKESVEAILGETLDLPKELAERADLPLLSHNLPADFAALRKLMMNHQ","num":1,"align_len":428}]}],"query_title":"fig|83333.1.peg.4","query_len":428,"query_id":"fig|83333.1.peg.4","stat":{"lambda":0.267,"hsp_len":86,"kappa":0.041,"entropy":0.14,"db_num":4309,"eff_space":335033460,"db_len":1350204}}},"version":"BLASTP 2.3.0+","program":"blastp","reference":"Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), \\"Gapped BLAST and PSI-BLAST: a new generation of protein database search programs\\", Nucleic Acids Res. 25:3389-3402."}}],{"kb|g.0.peg.4288":{"genome_name":"Escherichia coli K12","function":"Threonine synthase (EC 4.2.3.1)","genome_id":"kb|g.0"}}],"id":"8527913568541408","version":"1.1"}';
        # output = {'json_output': 'this is a sample output result'}
        jsonQueryReturn = json.loads(queryReturn)
        report = jsonQueryReturn['result'][0][0]['report']
        hitsList = report["results"]["search"]["hits"]
        # hits = self.formatHitList(hitsList)
        # print hits
        # hspsList = report["results"]["search"]["hits"][0]["hsps"]
        # hsps = self.formatHspList(hspsList)
        # print hsps

        returnVal = {
            "BlastOutput_db": "",
            "BlastOutput_program": "",
            "BlastOutput_query-ID": "",
            "BlastOutput_query-def": "",
            "BlastOutput_query-len": "",
            "BlastOutput_reference": "",
            "BlastOutput_version": "",
            "BlastOutput_param": {
                "Parameters": {
                    "Parameters_expect": "",
                    "Parameters_filter": "",
                    "Parameters_gap-extend": "",
                    "Parameters_gap-open": "",
                    "Parameters_matrix": "",
                    "Parameters_sc-match": "",
                    "Parameters_sc-mismatch": ""
                }
            },
            "BlastOutput_iterations": {
                "Iteration": [
                    {
                        "Iteration_hits": {
                            "Hit": self.formatHitList(hitsList)
                        },
                        "Iteration_iter-num": "",
                        "Iteration_query-ID": "",
                        "Iteration_query-def": "",
                        "Iteration_query-len": ""
                    }
                ]
            }
        }
        # print returnVal
        # Step 4 - save in workspace
        token = ctx['token']
        ws = workspaceService(self.workspaceURL, token=token)
        res = ws.save_objects(
            {
                "workspace": params["workspace_name"],
                "objects": [{
                    "type": "GenomeUtil.BlastOutput",
                    "data": returnVal,
                    "name": params["output_name"]
                }]
            }
        )
        #END blast_fasta

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method blast_fasta return value ' +
                             'output is not type dict as required.')
        # return the results
        return [returnVal]
