#python3

import subprocess
import gzip
import shutil
import sys
import os
import json
from MultiCodes import RunMultiCodes
from CombineBarSeq import RunCombineBarSeq
from HTMLReport import CreateHTMLString



def RunAll(inp_d):
    """
    inp_d: (d) Must contain

        MC_config_fp: (s) MultiCodes config json file path. MC_cgf_d must contain:
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
        CBS_config_fp: (s) CombineBarSeq config json file path. CBS_cfg_d contains:
            out_prefix_fp: (s) Output PoolCount/Colsum/Ignore File to write to
            pool_fp: (s) Input pool file to parse
            codes_fp_l: list<code_fp> List of all codes filepaths
                code_fp: (str) Path to codes file
            save_ignore: (b) If True, we save ignored lines to out_prefix_fp.ignored 

    """
   

    with open(inp_d["MC_config_fp"], "r") as f:
        MC_cfg_d = json.loads(f.read())

    pre_HTML_d = {}
    MultiCodesHTML_list = []
    codes_fp_list = []
    for fq_ind_d in inp_d["fq_index_list"]:
        MC_cfg_d["fastq_fp"] = fq_ind_d["fq_fp"]
        op_prefix = os.path.dirname(MC_cfg_d["out_prefix"]) \
                + os.path.basename(fq_ind_d["fq_fp"]).split(".")[0]
        MC_cfg_d["out_prefix"] = op_prefix
        codes_fp_list.append(op_prefix + ".codes")
        if "index_name" in fq_ind_d:
            MC_cfg_d.pop("indexfile_fp", None)
            MC_cfg_d["index_name"] = fq_ind_d["index_name"]
             
        elif "indexfile_fp" in fq_ind_d:
            MC_cfg_d.pop("index_name", None)
            MC_cfg_d["indexfile_fp"] = fq_ind_d["indexfile_fp"]
        else:
            raise Exception("Neither index_name nor indexfile_fp found in info dictionary")

        MC_report_d = RunMultiCodes(MC_cfg_d)
        MultiCodesHTML_list.append(MC_report_d)

    # For debugging and writing out to one HTML Report
    with open("tmp/MC_HTML.json", "w") as f:
        f.write(json.dumps(MultiCodesHTML_list, indent=4))

    print("Wrote MC HTML to tmp/MC_HTML.json")



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


def test_run_mult_cod():
    mult_cod_cfg_dict = create_test_vars()
    RunAll(mult_cod_cfg_dict)
    


def test():
    #test_run_mult_cod()
    #testCombineBarSeq()
    testHTMLReport()


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


def run_multi_codes_from_dict(mult_cod_cfg_dict):
    """
    multi codes config dict must contain following keys:

    out_fp_prefix,
    fastq_fp,
    index,

    Optional keys:
    debug
    """
    out_prefix = mult_cod_cfg_dict['out_fp_prefix']
    commands = [
            '-out',
            out_prefix,
            '-fastq_fp',
            mult_cod_cfg_dict['fastq_fp'],
            '-index',
            mult_cod_cfg_dict['index'],
            "-n25"
            ]

    if mult_cod_cfg_dict["debug"]:
        commands.append("-debug")

    #Note: first arg in commands is just a placeholder.
    vars_dict = multi_codes_imported_run(commands)

    out_dict = {
        "codes_file": out_prefix + ".codes",
        "counts_file": out_prefix + ".counts",
        "close_file": out_prefix + ".close",
        "mc_report_dict": vars_dict["report_dict"],
        "mc_report_str": vars_dict["report_str"]
            }

    return out_dict 
