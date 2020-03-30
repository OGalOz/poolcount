#python3

import json
import logging
from run_multi_codes import run_multi_codes_from_dict
from run_combine_barseq import run_combine_barseq_from_dict

def main():

    logging.basicConfig(level=logging.DEBUG)
    
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

def tester():

    return 0

if __name__ == "__main__":
    main()
