import unittest
import os
import json
import time

from os import environ
from ConfigParser import ConfigParser
from pprint import pprint

from biokbase.workspace.client import Workspace as workspaceService
from HomologySearch.HomologySearchImpl import HomologySearch


class HomologySearchTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        cls.ctx = {'token': token, 'provenance': [{'service': 'HomologySearch',
            'method': 'please_never_use_it_in_production', 'method_params': []}],
            'authenticated': 1}
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('HomologySearch'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = HomologySearch(cls.cfg)

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_HomologySearch_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_blast_fasta_to_database(self):
        params = {
            'workspace_name': self.getWsName(),
            'sequence': ">fig|83333.1.peg.4\nMKLYNLKDHNEQVSFAQAVTQGLGKNQGLFFPHDLPEFSLTEIDEMLKLDFVTRSAKILS\nAFIGDEIPQEILEERVRAAFAFPAPVANVESDVGCLELFHGPTLAFKDFGGRFMAQMLTH\nIAGDKPVTILTATSGDTGAAVAHAFYGLPNVKVVILYPRGKISPLQEKLFCTLGGNIETV\nAIDGDFDACQALVKQAFDDEELKVALGLNSANSINISRLLAQICYYFEAVAQLPQETRNQ\nLVVSVPSGNFGDLTAGLLAKSLGLPVKRFIAATNVNDTVPRFLHDGQWSPKATQATLSNA\nMDVSQPNNWPRVEELFRRKIWQLKELGYAAVDDETTQQTMRELKELGYTSEPHAAVAYRA\nLRDQLNPGEYGLFLGTAHPAKFKESVEAILGETLDLPKELAERADLPLLSHNLPADFAAL\nRKLMMNHQ\n",
            'program': 'blastp',
            'database': 'kbase_nr.faa',
            'evalue_cutoff': '1e-5',
            'max_hit': 10,
            'output_name': 'blast_output_0'
        }

        result = self.getImpl().blast_fasta(self.getContext(), params)
        ws = result[0]
        blast_outputs = self.getWsClient().get_objects([{'workspace': ws['workspaceName'], 'name': ws['blast_output_name']}])
        hits = blast_outputs[0]['data']['BlastOutput_iterations']['Iteration'][0]['Iteration_hits']['Hit']
        # print "hits", hits
        self.assertTrue(len(hits) >= 1)

    def test_blast_fasta_to_genomes(self):
        params = {
            'workspace_name': self.getWsName(),
            'sequence': ">fig|83333.1.peg.4\nMKLYNLKDHNEQVSFAQAVTQGLGKNQGLFFPHDLPEFSLTEIDEMLKLDFVTRSAKILS\nAFIGDEIPQEILEERVRAAFAFPAPVANVESDVGCLELFHGPTLAFKDFGGRFMAQMLTH\nIAGDKPVTILTATSGDTGAAVAHAFYGLPNVKVVILYPRGKISPLQEKLFCTLGGNIETV\nAIDGDFDACQALVKQAFDDEELKVALGLNSANSINISRLLAQICYYFEAVAQLPQETRNQ\nLVVSVPSGNFGDLTAGLLAKSLGLPVKRFIAATNVNDTVPRFLHDGQWSPKATQATLSNA\nMDVSQPNNWPRVEELFRRKIWQLKELGYAAVDDETTQQTMRELKELGYTSEPHAAVAYRA\nLRDQLNPGEYGLFLGTAHPAKFKESVEAILGETLDLPKELAERADLPLLSHNLPADFAAL\nRKLMMNHQ\n",
            'program': 'blastp',
            'database': '',
            'genome_ids': ['kb|g.0'],
            'search_type': 'features',
            'evalue_cutoff': '1e-5',
            'max_hit': 10,
            'output_name': 'blast_output_1'
        }

        result = self.getImpl().blast_fasta(self.getContext(), params)
        ws = result[0]
        blast_outputs = self.getWsClient().get_objects([{'workspace': ws['workspaceName'], 'name': ws['blast_output_name']}])
        hits = blast_outputs[0]['data']['BlastOutput_iterations']['Iteration'][0]['Iteration_hits']['Hit']
        # print "hits", hits
        self.assertTrue(len(hits) >= 1)
