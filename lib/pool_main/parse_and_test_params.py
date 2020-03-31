#python3
#This file is imported by poolcountImpl.py and checks params and poolfile to 
#   see if they're ready for.
import sys,os
sys.path.append(os.path.dirname(__file__))
from combineBarSeq import init_pool_dict


#params is a dict, as provided by the SDK
def parse_and_check_params(params):


    """
    fastq_files_refs_list should be a list with just refs
        as strings ['1/2/3', '3/2/1', ...]
    """
    parsed_params_dict = {
            'fastq_files_refs_list': params['fastq_refs'],
            'output_name': params['output_name']
            }
    

    return parsed_params_dict


def check_pool_file(pool_fp):

    #Parse pool file and check for errors
    test_vars_dict = {
            "poolfile": pool_fp,
            "report_dict": {
                "warnings": []
                }
    }
   
    try:
        init_pool_dict(test_vars_dict)
    except Exception:
        logging.warning("Pool file seems to have errors - " \
                + "Please check and reupload.")
        raise Exception




    return 0


