# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from installed_clients.KBaseReportClient import KBaseReport
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
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_poolcount
        report = KBaseReport(self.callback_url)

        #First we download all the Fastq Files and the Pool File
        multi_codes_config_fp = "sample_run_dict.json"
        pool_file_path = "Test_Files/feba_148_pool.n10"
        main_out_prefix = "tmp/cbs_out_Friday"

        with open(multi_codes_config_fp, "r") as f:
            mc_total_info = json.loads(f.read())

        report_dict = {} 
        report_str = ""

        mc_run_list = mc_total_info["multi_codes_run_list"]
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

        
        cmb_bs_dict = {"pool": pool_file_path, 
                "codes_filepaths_list": codes_fp_list,
                "main_out_prefix": main_out_prefix }
        logging.info("Running Combine Bar Seq")
        cmb_bs_out = run_combine_barseq_from_dict(cmb_bs_dict)
        report_dict["combine_bar_seq_report_dict"] = cmb_bs_out["cbs_report_dict"] 
        report_str += "---Combine BarSeq Report ---\n{}".format(cmb_bs_out[
            "cbs_report_str"])
        report_str += "---Combine BarSeq Warnings---\n{}".format(
                    "\n".join(cmb_bs_out["cbs_report_dict"]["warnings"]))

        logging.debug(report_str)
        poolcount_fp = cmb_bs_out["poolcount"]
        colsum_fp = cmb_bs_out["colsum"]

        return report_dict








        report_info = report.create({'report': {'objects_created':[],
                                                'text_message': params['parameter_1']},
                                                'workspace_name': params['workspace_name']})
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
