# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os, shutil

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from pool_main.parse_and_test_params import parse_and_check_params
from pool_main.downloader import download_fastq_and_prepare_mc
from pool_main.run_multi_codes import run_multi_codes_from_dict
from pool_main.run_combine_barseq import run_combine_barseq_from_dict
from pool_main.pool_util import clean_output_dir
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
        """
        # ctx is the context object
        # return variables are: output

        #BEGIN run_poolcount
        report = KBaseReport(self.callback_url)

        #We make an output directory in scratch:
        outputs_dir = os.path.join(self.shared_folder, "PoolCount_Outputs")
        os.mkdir(outputs_dir)


        parsed_params_dict = parse_and_check_params(params)

        dfu = DataFileUtil(self.callback_url) 

        # Download pool file using dfu:
        GetObjectsParams = {
                'object_refs': [parsed_params_dict['poolfile_ref']]
                }
        # We get the handle id
        PoolFileObjectData = dfu.get_objects(GetObjectsParams)['data'][0]['data']
        logging.info(PoolFileObjectData)
        poolfile_handle = PoolFileObjectData['poolfile']
        # We set the download path:
        poolfile_path = os.path.join(self.shared_folder, 
                "kb_poolcount_pool.n10")
        ShockToFileParams = {
                "handle_id": poolfile_handle,
                "file_path": poolfile_path,
                "unpack": "uncompress"
                }
        ShockToFileOutput = dfu.shock_to_file(ShockToFileParams)
        logging.info(ShockToFileOutput)



        fastq_dicts_list = download_fastq_and_prepare_mc(parsed_params_dict, 
                dfu, self.shared_folder, outputs_dir)


        logging.info(fastq_dicts_list)

        #First we download all the Fastq Files and the Pool File
        my_mod_dir = '/kb/module/lib/pool_main'
        #pool_file_path = os.path.join(my_mod_dir, "feba_148_pool.n10")

        
        """
        multi_codes_config_fp = os.path.join(my_mod_dir, "sample_run_dict.json")
        with open(multi_codes_config_fp, "r") as f:
            mc_total_info = json.loads(f.read())
        """

        report_dict = {} 
        report_str = ""

        mc_run_list = fastq_dicts_list
        mc_run_num = len(mc_run_list)
        logging.info("Total MultiCodes Runs: {}".format(mc_run_num))

        codes_fp_list = []
        for i in range(mc_run_num):
            mc_run_dict = mc_run_list[i]
            mc_out_dict = run_multi_codes_from_dict(mc_run_dict)
            logging.info("Completed {}/{} Multi Codes Runs".format(
            i+1, mc_run_num))
            report_dict["multi_codes_" + str(i+1)] = mc_out_dict["mc_report_dict"]
            report_str += "---Multi Codes {} Report----\n{}".format(i+1,mc_out_dict[
                "mc_report_str"])
            report_str += "---Multi Codes {} Warnings---\n{}".format(i+1,
                    "\n".join(mc_out_dict["mc_report_dict"]["warnings"]))
            codes_fp_list.append(mc_out_dict['codes_file'])

        cmbarseq_out_prefix = os.path.join(outputs_dir, 
                parsed_params_dict['output_name'])
        #combine barseq dict: 
        cmb_bs_dict = {"pool": poolfile_path, 
                "codes_filepaths_list": codes_fp_list,
                "main_out_prefix": cmbarseq_out_prefix }
        logging.info("Running Combine Bar Seq")
        cmb_bs_out = run_combine_barseq_from_dict(cmb_bs_dict)
        report_dict["combine_bar_seq_report_dict"] = cmb_bs_out["cbs_report_dict"] 
        report_str += "---Combine BarSeq Report ---\n{}".format(cmb_bs_out[
            "cbs_report_str"])
        report_str += "---Combine BarSeq Warnings---\n{}".format(
                    "\n".join(cmb_bs_out["cbs_report_dict"]["warnings"]))

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
            'label': 'RbTnSeqPoolCount_dir',
            'description': 'The folder containing outputs from this app'
                }

        extended_report_params = {
                'workspace_name': params['workspace_name'],
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
