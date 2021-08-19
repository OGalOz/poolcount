# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os, shutil, json

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from poct.parse_and_test_params import parse_and_check_params, get_FullRun_d
from poct.downloader import download_fastq_and_prepare_mc, download_mutantpool
from poct.FullProgram import PC_RunAll  
from poct.pool_util import clean_output_dir
from poct.upload_poolcount import upload_poolcount_to_KBase
#END_HEADER


class poolcount:
    '''
    Module Name:
    poolcount

    Module Description:
    A KBase module: poolcount
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.ws_url = config['workspace-url']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_poolcount(self, ctx, params):
        """
        This example function accepts any number of parameters and returns 
        results in a KBaseReport
        :param params: instance of mapping from String to unspecified object

        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        
        Args:
            params:
                ['workspace_name']: self.wsName,
                "mutantpool_ref": pool_ref (str),
                "fastq_files": list<fastq_refs (str)>,
                "genome_ref": genome_ref (str), 
                "KB_BarcodeCount_Bool": "yes"/"no" - create a poolcount file?
                "poolcount_description": (str) A text description of the pool file,
                "output_name": (str),
                ## "test_local_bool": test_local_bool Deprecated
                ## "save_ignore_bool": bool, Deprecated
                "maxReads": int or None,
                "minQuality": int,
                "debug": bool,
                "protocol_type": str,
                "doOff1": bool 

        """
        # ctx is the context object
        # return variables are: output
        logging.basicConfig(level=logging.DEBUG)
        logging.warning("INPUT PARAMS:")
        logging.warning(params)



        #BEGIN run_poolcount
        report = KBaseReport(self.callback_url)

        myToken = os.environ.get('KB_AUTH_TOKEN', None)
        ws = Workspace(self.ws_url, token=myToken)
        

        ws_id = ws.get_workspace_info({
                        'workspace': params['workspace_name']})[0]

        #Creating Data File Util Object
        dfu = DataFileUtil(self.callback_url)

        #We make the BarcodeCount output directory in scratch:
        outputs_dir = os.path.join(self.shared_folder, "BarcodeCount_Outputs")
        HTML_dir = os.path.join(self.shared_folder, "HTML_OP")
        MC_dir = os.path.join(self.shared_folder, "MC_op_dir")
        for x_dir in [outputs_dir, HTML_dir, MC_dir]:
            if os.path.isdir(x_dir):
                logging.info(f"{x_dir} contents: " + ",\n".join(
                              os.listdir(x_dir)))
            else:
                os.mkdir(x_dir)
        main_HTML_fp = os.path.join(HTML_dir, "index.html")


        # parsed_params_dict contains keys: 
        # mutantpool_ref, fastq_files_refs_list, genes_table_ref, output_name,
        # KB_BarcodeCount_Bool, poolcount_description
        parsed_params_dict = parse_and_check_params(params)
        # We get the username for later
        parsed_params_dict['username'] = ctx['user_id']


        # We set the mutantpool's path
        mutantpool_path = os.path.join(self.shared_folder, "kb_pool.pool")

        download_mutantpool(
            parsed_params_dict['mutantpool_ref'], mutantpool_path, dfu,
            genome_ref = parsed_params_dict['genome_ref'])
        
        poolcount_prefix = os.path.join(outputs_dir,
                                        parsed_params_dict["output_name"]
                                        )


        fastq_dicts_list = download_fastq_and_prepare_mc(parsed_params_dict, 
                dfu, self.shared_folder, outputs_dir)


        logging.info(fastq_dicts_list)

        #First we download all the Fastq Files and the Pool File

        report_dict = {} 
        report_str = ""

        mc_run_list = fastq_dicts_list
        mc_run_num = len(mc_run_list)
        logging.info("Total MultiCodes Runs: {}".format(mc_run_num))

        FullRun_d = get_FullRun_d(parsed_params_dict, MC_dir, fastq_dicts_list,
                                  poolcount_prefix, mutantpool_path,
                                  main_HTML_fp)

        # Running all programs:
        PC_RunAll(FullRun_d)


        # Now we upload the poolcount file to KBase to make a BarcodeCount Object
        if parsed_params_dict['KB_BarcodeCount_Bool']:
            upload_params = {
                    'username': parsed_params_dict['username'],
                    'fastq_refs': parsed_params_dict['fastq_files'],
                    'genome_ref': parsed_params_dict['genome_ref'],
                    'poolcount_description': parsed_params_dict[
                        'poolcount_description'] ,
                    'workspace_id': ws_id,
                    'ws_obj': ws,
                    'protocol_type': parsed_params_dict["protocol_type"],
                    'poolcount_fp': poolcount_prefix + ".poolcount",
                    'poolcount_name': parsed_params_dict['output_name'],
                    'dfu': dfu,
                    "scratch_dir": self.shared_folder,
                    "set_name": parsed_params_dict['output_name']
            }
            logging.info("UPLOADING BarcodeCount FILE to KBASE through DFU")
            upload_mutantpool_results = upload_poolcount_to_KBase(upload_params)
            logging.info("Upload BarcodeCount File Results:")
            logging.info(upload_mutantpool_results)


        # DEBUGGING WHAT REPORT DICT LOOKS LIKE
        with open(os.path.join(self.shared_folder, "Report.JSON"), 'w') as g:
            g.write(json.dumps(report_dict, indent=2))


        #Cleaning outputs dir (removing .codes, .close, .counts files)
        clean_output_dir(outputs_dir, self.shared_folder)


        #Writing report string:
        report_fp = os.path.join(outputs_dir, "Run_Report.txt")
        with open(report_fp, "w") as f:
            f.write(report_str)

        #Return Outputs Dir to user:
        dir_zip_shock_id = dfu.file_to_shock({"file_path": outputs_dir,
            'pack': 'zip'})['shock_id']
        dir_link_dict = {
            'shock_id': dir_zip_shock_id,
            'name': parsed_params_dict['output_name'] + ".zip",
            'label': 'RBTnSeqBarcodeCount_dir',
            'description': 'The folder containing outputs from this app'
        }




        HTML_report_shock_id = dfu.file_to_shock({
                "file_path": HTML_dir,
                "pack": "zip"
                })['shock_id']

        HTML_report_d_l = [{"shock_id": HTML_report_shock_id,
                            "name": os.path.basename(os.path.join(HTML_dir,"index.html")),
                            "label": "BarcodeCount Report",
                            "description": "HTML Summary Report for MultiCodes and Combine BarSeq"
                            }]




        extended_report_params = {
                'workspace_name': params['workspace_name'],
                "html_links": HTML_report_d_l,
                "direct_html_link_index": 0,
                "html_window_height": 333,
                "report_object_name": "KB_BarcodeCount_Report",
                'file_links' : [dir_link_dict],
                'message': ""
                }
        report_info = report.create_extended_report(extended_report_params)

        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
        }
        #END run_poolcount

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_poolcount return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
