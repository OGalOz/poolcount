# python3
# This file is to upload an exps file to KBasePoolTSV.Experiments

import logging
import re
import shutil
import os


def upload_expsfile_to_KBase(up):
    """
    upload params (up) must include the following keys:
        genome_ref (str),
        expsfile_description (str),
        workspace_id (str),
        ws_obj (obj),
        expsfile_fp (str),
        expsfile_op_name (str),
        dfu (obj),
        scratch_dir (str)
    """

    # We check if the expsfile file is in scratch as it should be
    up["expsfile_fp"] = check_expsfile_in_scratch(
        up["scratch_dir"], up["expsfile_fp"], up["expsfile_op_name"]
    )

    # We check correctness of expsfile
    column_header_list, num_rows = check_exps_file(up["expsfile_fp"])

    # We create the KBase handle for the object:
    file_to_shock_result = up["dfu"].file_to_shock(
        {"file_path": up["expsfile_fp"], "make_handle": True, "pack": "gzip"}
    )

    # The following var res_handle only created for simplification of code
    res_handle = file_to_shock_result["handle"]

    # We create the data for the object
    exps_data = {
        "file_type": "KBasePoolTSV.Experiments",
        "expsfile": res_handle["hid"],
        # below should be shock
        "handle_type": res_handle["type"],
        "shock_url": res_handle["url"],
        "shock_node_id": res_handle["id"],
        "compression_type": "gzip",
        "file_name": res_handle["file_name"],
        "column_header_list": column_header_list,
        "num_lines": str(num_rows),
        "related_genome_ref": up["genome_ref"],
        "related_organism_scientific_name": get_genome_organism_name(
            up["genome_ref"], up["ws_obj"]
        ),
        "description": up["expsfile_description"],
    }

    # To get workspace id:
    ws_id = up["workspace_id"]

    save_object_params = {
        "id": ws_id,
        "objects": [
            {
                "type": "KBasePoolTSV.Experiments",
                "data": exps_data,
                "name": up["expsfile_op_name"],
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


def check_exps_file(expsfile_fp):

    required = ["SetName", "Index", "Description", "Date_pool_expt_started"]

    cols, num_rows = read_table(expsfile_fp, required)

    return [cols, num_rows]


def read_table(fp, required):
    """
    Following function takes a filename and a list of required fields i
    (file is TSV)
    returns list of headers
    Does not return header line
    """
    with open(fp, "r") as f:
        file_str = f.read()
    file_list = file_str.split("\n")
    header_line = file_list[0]
    # Check for Mac Style Files
    if re.search(r"\t", header_line) and re.search(r"\r", header_line):
        raise Exception(
            (
                "Tab-delimited input file {} is a Mac-style text file "
                "which is not supported.\n"
                "Use\ndos2unix -c mac {}\n to convert it to a Unix "
                "text file.\n"
            ).format(fp, fp)
        )
    cols = header_line.split("\t")
    cols_dict = {}
    for i in range(len(cols)):
        cols_dict[cols[i]] = i
    for field in required:
        if field not in cols_dict:
            raise Exception(
                "No field {} in {}. Must include fields".format(field, fp)
                + "\n{}".format(" ".join(required))
            )
    rows = []
    for i in range(1, len(file_list)):
        line = file_list[i]
        # if last line empty
        if len(line) == 0:
            continue
        line = re.sub(r"[\r\n]+$", "", line)
        split_line = line.split("\t")
        if not len(split_line) == len(cols):
            raise Exception(
                "Wrong number of columns in:\n{}\nin {} l:{}".format(line, fp, i)
            )
        new_dict = {}
        for i in range(len(cols)):
            new_dict[cols[i]] = split_line[i]
        rows.append(new_dict)

    return [cols, len(file_list)]


def check_expsfile_in_scratch(scratch_dir, expsfile_fp, expsfile_op_name):
    # We make sure the expsfile file is in the scratch dir
    if scratch_dir not in expsfile_fp:
        logging.info("expsfile file not in scratch directory- moving it.")
        new_exps_fp = os.path.join(scratch_dir, expsfile_op_name)
        shutil.copyfile(expsfile_fp, new_exps_fp)
    else:
        new_exps_fp = expsfile_fp

    return new_exps_fp


def get_genome_organism_name(genome_ref, ws_obj):
    res = ws_obj.get_objects2(
        {"objects": [{"ref": genome_ref, "included": ["scientific_name"]}]}
    )
    scientific_name = res["data"][0]["data"]["scientific_name"]
    return scientific_name
