#python3



import subprocess
import gzip
import shutil
import sys
import os
import json
import logging
from poct.MultiCodes import RunMultiCodes
from poct.CombineBarSeq import RunCombineBarSeq
from poct.HTMLReport import CreateHTMLString



def PC_RunAll(inp_d):
    """
    inp_d: (d) Must contain

        MC_config_d: (d) MultiCodes config json file path. MC_cgf_d must contain:
            out_prefix: str,
            maxReads: int or None, 
            index_name: str,
            minQuality: int, 
            debug: bool,
            protocol_type: str,
            bs3_fp: File path to barseq3.index2 file
            doOff1: bool, 
            MC_seqs:
                    dnt_pre: GTCTCGTAG,
                    dnt_post: CGATGAATT,
                    bs_pre: CAGCGTACG,
                    bs_post: AGAGACCTC
            fastq_fp: str 
        fq_index_list: (l) list<fq_ind_d> where
            fq_ind_d: (d) FASTQ INDEX DICT
                fq_fp: (str) Fastq file path
                [index_name] OR (str)
                [indexfile_fp] (str) 
                debug:
        CBS_config_d: (d) CombineBarSeq config json file path. CBS_cfg_d contains:
            out_prefix_fp: (s) Output PoolCount/Colsum/Ignore File to write to
            pool_fp: (s) Input pool file to parse
            save_ignore: (b) If True, we save ignored lines to out_prefix_fp.ignored 
            <<codes_fp_l>>: list<code_fp> List of all codes filepaths GENERATED LATER
                code_fp: (str) Path to codes file
        HTML_op_fp: (s) Path to write HTML file to

    Description:
        In this function we run multiple scripts.
    """

    logging.info("Beginning PC Run All")
   
    '''
    with open(inp_d["MC_config_fp"], "r") as f:
        MC_cfg_d = json.loads(f.read())
    '''
    MC_cfg_d = inp_d["MC_config_d"]

    pre_HTML_d = {}
    MultiCodesHTML_list = []
    codes_fp_list = []
    logging.info("Starting to Run MultiCodes.")
    for fq_ind_d in inp_d["fq_index_list"]:
        # Every time we update the entire Multi Codes run dict
        MC_cfg_d["fastq_fp"] = fq_ind_d["fq_fp"]
        op_prefix = os.path.dirname(MC_cfg_d["out_prefix"]) \
                + os.path.basename(fq_ind_d["fq_fp"]).split(".")[0]
        MC_cfg_d["out_prefix"] = op_prefix

        logging.info("Starting to Run MultiCodes with codes destination: " + \
                     op_prefix + ".codes")
        codes_fp_list.append(op_prefix + ".codes")
        if "index_name" in fq_ind_d:
            # We remove indexfile_fp key if it formerly existed
            MC_cfg_d.pop("indexfile_fp", None)
            MC_cfg_d["index_name"] = fq_ind_d["index_name"]
             
        elif "indexfile_fp" in fq_ind_d:
            # We remove index_name key if it formerly existed exists
            MC_cfg_d.pop("index_name", None)
            MC_cfg_d["indexfile_fp"] = fq_ind_d["indexfile_fp"]
        else:
            raise Exception("Neither index_name nor indexfile_fp found in info dict")

        MC_report_d = RunMultiCodes(MC_cfg_d)
        MultiCodesHTML_list.append(MC_report_d)


    '''
    # For debugging and writing out to one HTML Report
    with open("tmp/MC_HTML.json", "w") as f:
        f.write(json.dumps(MultiCodesHTML_list, indent=4))
    logging.info("Wrote MC HTML to tmp/MC_HTML.json")
    '''

    '''
    with open(inp_d["CBS_config_fp"], "r") as f:
        CBS_d = json.loads(f.read())
    '''
    CBS_d = inp_d['CBS_config_d']

    CBS_d["codes_fp_l"] = codes_fp_list

    logging.info("Running CombineBarSeq")
    ctg_d = RunCombineBarSeq(CBS_d)


    '''
    with open("tmp/CBS_HTML.json", "w") as f:
        f.write(json.dumps(ctg_d, indent=4))
    logging.info("Wrote CBS HTML to tmp/CBS_HTML.json")
    '''


    inp_HTML_d = {
            "MultiCodes_reports_list": MultiCodesHTML_list,
            "CombineBarSeqReport_d": ctg_d
    }

    logging.info("Writing HTML.")
    HTML_str = CreateHTMLString(inp_HTML_d)

    with open(inp_d['HTML_op_fp'], "w") as f:
        f.write(HTML_str)

    logging.info("Finished running PoolCount's RunAll")

    #End this part of program
    return None
    




def test_run_mult_cod():
    mult_cod_cfg_dict = create_test_vars()
    RunAll(mult_cod_cfg_dict)
    


def test():
    #test_run_mult_cod()
    #testCombineBarSeq()
    #testHTMLReport()
    print("Not doing anything")


def main():
    test()


if __name__ == "__main__":
    main()



# OLD
# both inputs str paths
def unzip_fastq(zipped_fastq_fp, fastq_dest_path):
    with gzip.open(zipped_fastq_fp, 'rb') as f_in:
        with open(fastq_dest_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)



def testCombineBarSeq():
    CBS_d = {

            'out_prefix_fp': '/usr2/people/omreeg/RBTNSEQ_PROG/PoolCount/tmp/CBS_Test',
            'pool_fp': '/usr2/people/omreeg/RBTNSEQ_PROG/PoolCount/DATA/CL/pool.n10',
            'codes_fp_l': [
                            '/usr2/people/omreeg/RBTNSEQ_PROG/PoolCount/tmp/Test1FEBA_BS_214_IT068_AGGCCG_S68_L006_R1_001.codes'
                ],
            'save_ignore': True
            }

    ctg_d = RunCombineBarSeq(CBS_d)

    with open("tmp/CBS_HTML.json", "w") as f:
        f.write(json.dumps(ctg_d, indent=4))
    print("Wrote CBS HTML to tmp/CBS_HTML.json")



def testHTMLReport():
    
    HTML_op = "tmp/Test_Html_Report.html"

    with open("tmp/MC_HTML.json", "r") as f:
        MC_rep_l = json.loads(f.read())

    with open("tmp/CBS_HTML.json", "r") as f:
        CBS_rep_d = json.loads(f.read())


    inp_HTML_d = {
            "MultiCodes_reports_list": MC_rep_l,
            "CombineBarSeqReport_d": CBS_rep_d
    }

    HTML_str = CreateHTMLString(inp_HTML_d)

    with open(HTML_op, "w") as f:
        f.write(HTML_str)




def create_test_vars():
    mult_cod_cfg_dict = {
            "MC_config_fp": "/usr2/people/omreeg/RBTNSEQ_PROG/PoolCount/SampleMC_inp.json",
            "fq_index_list": [
                {
                    "fq_fp": "/usr2/people/omreeg/RBTNSEQ_PROG/PoolCount/DATA/CL/FEBA_BS_214_IT068_AGGCCG_S68_L006_R1_001.fastq",
                    "index_name": "IT068"
                    }
                ]
            }

    return mult_cod_cfg_dict
