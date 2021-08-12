#python3
#This file is imported by poolcountImpl.py and checks params and poolfile to 
#   see if they're ready for run.
import sys, os, logging
import re
sys.path.append(os.path.dirname(__file__))
#from combineBarSeq import init_pool_dict


#params is a dict, as provided by the SDK
def parse_and_check_params(params):
    """
    Args:
        params: (d)
            All 'refs' look like 'A/B/C' e.g. '63336/2/1'
            "poolfile_ref": pool_ref (str),
            "fastq_files": list<fastq_refs (str)>,
            "genome_ref": genome_ref (str), 
            "KB_PoolCount_Bool": "yes"/"no" - create a poolcount file?
            "poolcount_description": (str) A text description of the pool file,
            "output_name": (str),
            "test_local_bool": test_local_bool
            ['workspace_name']: self.wsName,
            "save_ignore_bool": bool,
            "max_Reads": int or -1 or None,
            "minQuality": int,
            "debug": bool,
            "protocol_type": str,
            "doOff1": bool 
    """

    logging.warning(params)

    for p in ["poolfile_ref", "fastq_files", "KB_PoolCount_Bool",
            "poolcount_description", "genome_ref", "output_name",
            "test_local_bool", "save_ignore_bool", "max_Reads",
            "minQuality", "debug", "protocol_type", "doOff1"]:
        if p not in params:
            raise Exception(p + " not found in params")
        elif params[p] in ["true", "True"]:
            params[p] = True
        elif params[p] in ["false", "False"]:
            params[p] = False

    for ref in [params['poolfile_ref'], params['genome_ref']] + params['fastq_files']:
        if len(ref.split('/')) != 3:
            raise Exception('ref format not A/B/C as expected: ' + ref)


    if params["KB_PoolCount_Bool"] in ["yes","Yes"]:
        kb_pc_bool = True
    else:
        logging.info("Will not create KBase pool object due to user param.")
        kb_pc_bool = False 



    # Duplicating params dict
    parsed_params_dict = {x:params[x] for x in params.keys()}
    # Updating certain keys
    parsed_params_dict["KB_PoolCount_Bool"] =  kb_pc_bool
    parsed_params_dict["max_Reads"] = params["max_Reads"] if params["max_Reads"] not in [-1, "-1"] else None
    for x in ["preseq", "postseq"]:
        if x in parsed_params_dict:
            if parsed_params_dict[x] == "None":
                parsed_params_dict[x] = None


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
    if ' ' in op_name:
        raise Exception("Output Set Name Cannot Contain Spaces: {}".format(
                                                                    op_name))
    op_name = op_name.replace(' ', '_')
    rgx = re.search(r'[^\w]', op_name)
    if rgx:
        logging.warning("Non-alphanumeric character in output name: " + rgx[0])
        op_name = "Default_Name_Check_Chars"
    return op_name


#pool is a dict
def init_pool_dict(vars_dict):
    pool = {} # pool dict is rcbarcode to [barcode, scaffold, strand, pos]
    with open(vars_dict['poolfile'], "r") as f:
        poolfile_str = f.read()
        poolfile_lines = poolfile_str.split("\n")
        for pool_line in poolfile_lines:
            pool_line.rstrip()
            pool = check_pool_line_and_add_to_pool_dict(pool_line, pool,
                    vars_dict)
    if len(pool.keys()) == 0:
        raise Exception("No entries in pool file")

    return pool


def check_pool_line_and_add_to_pool_dict(pool_line, pool, vars_dict):

    #We get first 7 columns of pool_line (out of 12)
    split_pool_line = pool_line.split("\t")[:7]

    #We remove spaces:
    for x in split_pool_line:
        x.replace(' ', '')

    if len(split_pool_line) >= 7:
        #We unpack
        barcode, rcbarcode, undef_1, undef_2, scaffold, strand, pos = split_pool_line
    else:
        warning_text = "pool file line with less than 7 tabs:\n{}".format(
            pool_line)
        vars_dict["report_dict"]["warnings"].append(warning_text)
        logging.warning(warning_text)
        barcode = "barcode"

    if barcode == "barcode": #Header line
        pass
    else:
        if not re.search(r"^[ACGT]+$", barcode):
            logging.debug(len(barcode))
            raise Exception("Invalid barcode: |{}|".format(barcode))
        if not re.search( r"^[ACGT]+$",rcbarcode ):
            raise Exception("Invalid rcbarcode: |{}|".format(rcbarcode))
        if not (pos == "" or re.search( r"^\d+$", pos)):
            raise Exception("Invalid position: |{}|".format(pos))
        if not (strand == "+" or strand == "-" or strand == ""):
            raise Exception("Invalid strand: |{}|".format(strand))
        if rcbarcode in pool:
            raise Exception("Duplicate rcbarcode.")
        pool[rcbarcode] = [barcode, scaffold, strand, pos]
    return pool
