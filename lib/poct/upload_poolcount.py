# python3
# This file is to upload a pool file to KBasePoolTSV.PoolFile

import logging
import shutil
import os
import datetime


def upload_poolcount_to_KBase(up):
    """
    upload params (up) must include the following keys:

        username (str),
        genes_table_ref (str),
        fastq_refs (list),
        poolcount_description (str),
        protocol_type (str)
        workspace_id (str),
        ws_obj (obj),
        poolcount_fp (str),
        poolcount_name (str),
        dfu (obj),
        scratch_dir (str)
        set_name (str)
    """

    # We check if the poolcount file is in scratch as it should be
    up["poolcount_fp"] = check_poolcount_in_scratch(
        up["scratch_dir"], up["poolcount_fp"], up["poolcount_name"]
    )

    # We check correctness of poolcount file
    column_header_list, num_lines = check_poolcount_file(up["poolcount_fp"])

    # We create the KBase handle for the object:
    file_to_shock_result = up["dfu"].file_to_shock(
        {"file_path": up["poolcount_fp"], "make_handle": True, "pack": "gzip"}
    )
    # The following var res_handle only created for simplification of code
    res_handle = file_to_shock_result["handle"]

    # We create an updated description (prefix) with username and time
    date_time = datetime.datetime.utcnow()
    updated_description = "Created by {} on {}\n".format(up["username"], str(date_time))

    # We create the data for the object
    poolcount_data = {
        "file_type": "KBaseRBTnSeq.RBTS_MutantPoolCount",
        "poolcount": res_handle["hid"],
        # below should be shock
        "handle_type": res_handle["type"],
        "shock_url": res_handle["url"],
        "shock_node_id": res_handle["id"],
        "compression_type": "gzip",
        "column_header_list": column_header_list,
        "column_headers_str": ", ".join(column_header_list),
        "fastqs_used": up["fastq_refs"],
        "fastqs_used_str": "NA - use fastqs_used",
        "file_name": res_handle["file_name"],
        "utc_created": str(date_time),
        "set_name": up["set_name"],
        "num_lines": str(num_lines),
        "protocol_used": up["protocol_type"],
        "related_genome_ref": up["genome_ref"],
        "related_organism_scientific_name": get_genome_organism_name(
            up["genome_ref"], up["ws_obj"]
        ),
        "description": updated_description + up["poolcount_description"],
    }

    # To get workspace id:
    ws_id = up["workspace_id"]
    save_object_params = {
        "id": ws_id,
        "objects": [
            {
                "type": "KBaseRBTnSeq.RBTS_MutantPoolCount",
                "data": poolcount_data,
                "name": up["poolcount_name"],
            }
        ],
    }
    # save_objects returns a list of object_infos
    dfu_object_info = up["dfu"].save_objects(save_object_params)[0]
    print("dfu_object_info: ")
    print(dfu_object_info)
    return {
        "Name": dfu_object_info[1],
        "Type": dfu_object_info[2],
        "Date": dfu_object_info[3],
    }


def check_poolcount_file(poolcount_fp):
    """
    We check the pool file by initializing into dict format

    The function "init_pool_dict" runs the tests to see if the file is
    correct.
    """
    # Expected fields
    exp_f = "barcode rcbarcode scaffold strand pos".split(" ")

    with open(poolcount_fp, "r") as f:
        f_str = f.read()
    f_list = f_str.split("\n")
    num_lines = len(f_list)
    header_line = f_list[0]

    if header_line == "":
        raise Exception("File format incorrect: " + poolcount_fp)

    fields = header_line.split("\t")

    if not (len(fields) >= 6):
        raise Exception("Too few fields in " + poolcount_fp)
    for i in range(len(exp_f)):
        if not fields[i] == exp_f[i]:
            raise Exception("Expected {} but field is {}".format(exp_f[i], fields[i]))

    return [fields, num_lines]


def check_poolcount_in_scratch(scratch_dir, poolcount_fp, poolcount_name):
    # We make sure the poolcount file is in the scratch dir
    if scratch_dir not in poolcount_fp:
        logging.info("poolcount file not in scratch directory- moving it.")
        new_pool_fp = os.path.join(scratch_dir, poolcount_name)
        shutil.copyfile(poolcount_fp, new_pool_fp)
    else:
        new_pool_fp = poolcount_fp

    return new_pool_fp


def get_genome_organism_name(genome_ref, ws):
    # Getting the organism name using WorkspaceClient
    res = ws.get_objects2(
        {
            "objects": [
                {
                    "ref": genome_ref,
                    "included": ["scientific_name"],
                }
            ]
        }
    )
    scientific_name = res["data"][0]["data"]["scientific_name"]
    return scientific_name



def get_genome_organism_name_from_genes_table(self, gene_table_ref):
    # Getting the organism name using WorkspaceClient
    ws = self.params['ws_obj'] 
    res = ws.get_objects2(
        {
            "objects": [
                {
                    "ref": gene_table_ref,
                    "included": ["related_organism_scientific_name"],
                }
            ]
        }
    )
    scientific_name = res["data"][0]["data"]["related_organism_scientific_name"]
    return scientific_name




