# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from poolcount.poolcountImpl import poolcount
from poolcount.poolcountServer import MethodContext
from poolcount.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class poolcountTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('poolcount'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'poolcount',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = poolcount(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_your_method(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods
       
        pool_ref = "58816/36/1"


        set_output_name = "Burk376_Test1"

        genome_ref = "58816/34/1"

        fastq_ref_1 = '58816/39/1' 
        fastq_ref_2 = '58816/37/1'

        fastq_refs = [fastq_ref_1, fastq_ref_2]

        test_local_bool = False 

        # added
        save_ignore_bool = True

        # added
        max_Reads = None

        # added
        minQuality = 0

        debug = False

        # 'custom' (requires preseq/postseq) or 'dntag' (requires index_file)
        # or 'base' or 'bs3' or 'n25' or 'Unknown'
        protocol_type = "n25"
        # Below only if protocol_type == 'custom'
        preseq = None
        postseq = None

        doOff1 = False

        ret = self.serviceImpl.run_poolcount(self.ctx, 
                {
                'workspace_name': self.wsName,
                "poolfile_ref": pool_ref,
                "fastq_files": fastq_refs,
                "genome_ref": genome_ref, 
                "KB_PoolCount_Bool": "yes",
                "poolcount_description": "Testing",
                "output_name": set_output_name,
                "test_local_bool": test_local_bool,
                "save_ignore_bool": save_ignore_bool,
                "max_Reads": max_Reads,
                "minQuality": minQuality,
                "debug": debug,
                "protocol_type": protocol_type,
                "doOff1": doOff1,
                "preseq": preseq,
                "postseq": postseq
                    })
