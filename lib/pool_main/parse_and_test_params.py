#python3
#This file is imported by poolcountImpl.py and checks params and poolfile to 
#   see if they're ready for.
import sys,os,logging
import re
sys.path.append(os.path.dirname(__file__))
from combineBarSeq import init_pool_dict


#params is a dict, as provided by the SDK
def parse_and_check_params(params):

    logging.warning(params)

    for p in ["poolfile_ref", "fastq_files", "KB_PoolCount_Bool",
            "poolcount_description", "genome_ref", "output_name"]:
        if p not in params:
            raise Exception(p + " not found in params")

    for ref in [params['poolfile_ref'], params['genome_ref']] + params['fastq_files']:
        if len(ref.split('/')) != 3:
            raise Exception('ref format not A/B/C as expected: ' + ref)


    if params["KB_PoolCount_Bool"] in ["yes","Yes"]:
        kb_pc_bool = True
    else:
        logging.info("Will not create KBase pool object due to user param.")
        kb_pc_bool = False 

    parsed_params_dict = {
            "poolfile_ref": params["poolfile_ref"],
            'fastq_files_refs_list': params['fastq_files'],
            "genome_ref": params["genome_ref"],
            'output_name': check_output_name(params['output_name']),
            "KB_PoolCount_Bool": kb_pc_bool,
            "poolcount_description": params["poolcount_description"]
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

# op_name string, (output_name)
def check_output_name(op_name):
    op_name = op_name.replace(' ', '_')
    rgx = re.search(r'[^\w]', op_name)
    if rgx:
        logging.warning("Non-alphanumeric character in output name: " + rgx[0])
        op_name = "Default_Name_Check_Chars"
    return op_name
