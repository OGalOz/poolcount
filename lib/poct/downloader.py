#python3

#File used to download
#We use DataFileUtil: 
# https://github.com/kbaseapps/DataFileUtil/blob/master/DataFileUtil.spec
import sys,os,logging,re

#parsed_params_dict is a variable dict - must have keys in func
# dfu is data file util object
# scratch_dir is str path to scratch file
#outputs_dir path to where outputs stored 
def download_fastq_and_prepare_mc(parsed_params_dict, dfu, scratch_dir, 
        outputs_dir):
    """
    Args:
        parsed_params_dict:
            fastq_files: list<ref>
                    ref: str 'A/B/C'
    """
    new_fastq_fps = []

    fastq_file_refs = parsed_params_dict['fastq_files']
    dfu_dict = {
        "object_refs": fastq_file_refs,
        "ignore_errors": False
            }
    #The following gets us the shock i
    get_objects_results = dfu.get_objects(dfu_dict)
    file_data_list = get_objects_results['data']
    #file_data is a dict
    for i in range(len(file_data_list)):
        file_data = file_data_list[i]
        logging.info("filedata")
        logging.info(file_data)
        orig_fq_fn = file_data['data']['lib']['file']['file_name']
        fq_shock_id = file_data['data']['lib']['file']['id']
        #Downloading the fastq itself
        result_filename = convert_fastq_filename(orig_fq_fn)
        logging.info("Downloading FASTQ File")
        fq_fp = os.path.join(scratch_dir, result_filename)
        new_fastq_fps.append(fq_fp)
        fastq_download_params = {"shock_id": fq_shock_id, 
                "file_path": fq_fp, 'unpack': 'unpack'}
        fastq_file_info = dfu.shock_to_file(fastq_download_params)

    fastq_dicts_list = convert_fastq_fp_list_to_add_index(new_fastq_fps, 
            outputs_dir)

    return fastq_dicts_list 

#original_fastq_fn: (str)
def convert_fastq_filename(original_fastq_fn):
    return ".".join(original_fastq_fn.split(".")[:2])


# Gets poolfile path
def download_poolfile(poolfile_ref, poolfile_path, dfu, genome_ref=None):

    GetObjectsParams = {
            'object_refs': [poolfile_ref]
            }

    # We get the handle id
    PoolFileObjectData = dfu.get_objects(GetObjectsParams)['data'][0]['data']
    logging.info("DFU Pool File Get objects results:")
    logging.info(PoolFileObjectData)

    poolfile_handle = PoolFileObjectData['poolfile']
    related_genome_ref = PoolFileObjectData['related_genome_ref']
    if genome_ref is not None:
        if related_genome_ref != genome_ref:
            raise Exception("The genome ref associated with the poolfile is "
                            f"'{related_genome_ref}', whereas the current genome "
                            f"ref is '{genome_ref}'. The two must match.")
        else:
            logging.info("The pool file's genome ref and the current genome ref given match.")


    # Set params for shock to file
    ShockToFileParams = {
            "handle_id": poolfile_handle,
            "file_path": poolfile_path,
            "unpack": "uncompress"
            }
    ShockToFileOutput = dfu.shock_to_file(ShockToFileParams)
    logging.info(ShockToFileOutput)
    # Poolfile is at location "poolfile_path"

    return None 


#Input is a list of fastq fps strings
#outputs dir str (not in use)
def convert_fastq_fp_list_to_add_index(new_fastq_fps, outputs_dir):


    fq_fp_dicts = []

    for fq_fp in new_fastq_fps:
        fq_fp_dicts.append(get_index_val(fq_fp, outputs_dir))

    return fq_fp_dicts


def get_index_val(fq_fp, outputs_dir ):
    """

    Output needs to look like:
        fq_ind_d: (d) FASTQ INDEX DICT
            fq_fp: (str) Fastq file path
            [index_name] OR (str)
            [indexfile_fp] (str) FOR NOW ONLY index_name
            debug:

            Last 3 keys not in use
    """

    fn = fq_fp.split('/')[-1]
    out_fp_prefix = fn.split('.')[0]
    
    match_it = re.search(r'_IT\d\d\d_', fn)
    match_Index = re.search(r'_Index\d+_', fn)
    match_s = re.search(r'_S\d+_', fn)
    if match_it:
        index_type = "IT"
        index = match_it[0][1:-1]
    elif match_Index:
        index_type = "Index"
        index = match_Index[0][1:-1]
    elif match_s:
        index_type = "S"
        index = match_s[0][1:-1]
    else:
        raise Exception("Could not recognize index: not in IT### or S# form: " \
                + "{}".format(fn))


    fq_fp_dict = {
            "fq_fp": fq_fp,
            "index_name": index,
            "debug": False,
            "index_type": index_type,
            "index_val": index,
            "out_fp_prefix": os.path.join(outputs_dir, out_fp_prefix),
    }

    return fq_fp_dict
    


def get_IT_S_indexes():

    index_suffix_rough = [str(x) for x in range(1000)]
    IT_indexes = []
    S_indexes = []
    for ruf in index_suffix_rough:
        if len(ruf) == 1:
            IT_ind_suffixes.append('IT00' + ruf)
            S_indexes.append('S' + ruf)
        elif len(ruf) == 2:
            IT_ind_suffixes.append('IT0' + ruf)
            S_indexes.append('S' + ruf)
        elif len(ruf) == 3:
            IT_ind_suffixes.append('IT' + ruf)
            S_indexes.append('S' + ruf)
        else:
            raise Exception("Index nums out of range")

    return [IT_indexes, S_indexes]




