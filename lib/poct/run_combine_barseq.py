#python3
import sys,os
sys.path.append(os.path.dirname(__file__))
from combineBarSeq import combine_barseq_imported_run

def run_combine_barseq_from_dict(cmb_bs_dict):

    """
    cmb_bs_dict should include:
        "pool": pool_file_path (str)
        "codes_filepaths_list": a list of (str) each leading to a code file
        "main_out_prefix": (str) prefix for .poolcount & .colsum files
    """
    #Codes files are given in one string separated by commas
    codes_list_str = ",".join(cmb_bs_dict['codes_filepaths_list'])
    combine_barseq_args = [
        cmb_bs_dict["main_out_prefix"],
        cmb_bs_dict["pool"],
        codes_list_str
            ]

    #Where combineBarseq runs
    vars_dict = combine_barseq_imported_run(combine_barseq_args)

    cmb_bs_out = {
    "poolcount": cmb_bs_dict["main_out_prefix"] + ".poolcount",
    "colsum": cmb_bs_dict["main_out_prefix"] + ".colsum",
    "cbs_report_dict": vars_dict["report_dict"],
    "cbs_report_str": vars_dict["report_str"]

    }
    
    return cmb_bs_out


