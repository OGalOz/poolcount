#python3

import os
import logging
import shutil

def clean_output_dir(output_dir, scratch_dir):
    """moves .codes, .close, .counts to scratch dir"""
    
    outputs_dir_files = os.listdir(output_dir)

    logging.info("Outputs Dir: \n{}".format("\n".join(outputs_dir_files)))

    files_to_remove = []

    for f in outputs_dir_files:
        if f.split('.')[-1] in ["codes","counts","close"]:
            files_to_remove.append(f)

    for f in files_to_remove:
        filepath = os.path.join(output_dir, f)
        new_filepath = os.path.join(scratch_dir, f)
        shutil.move(filepath, new_filepath)

    logging.info("files removed:\n{}".format("\n".join(files_to_remove)))

    return None



def convert_from_poolfile_to_sequence_set_and_back(inp_fp_path, 
        op_path, conversion_type, description="", run_id=None):
    """
    In this function we take either pool file or Sequence Set
        and convert from one to the other. Sequence Set is output
        and input as a JSON file, pool file is output and input
        as a TSV file. Conversion_type is an int:
        0 -> going from pool file to Sequence Set
        1 -> going from Sequence Set to pool file

    Args:
        inp_fp_path: Path to either pool file (TSV) or Sequence Set (JSON)
        op_path: Output path to the conversion output
        conversion_type: (int) from [0,1]. 0 -> Pool file to Sequence Set
            1 -> Sequence Set to Pool File.
        description: (string) Optional string describing the set

    Poolfile:
           barcode(str):rcbarcode(str):nTot(int):n(int):scaffold(str):strand(str +/-):pos(int):
           n2(int):scaffold2(int):strand2(str +/-):pos2(int):nPastEnd(int)

    """

    if conversion_type not in [0,1]:
        raise Exception("Cannot recognize conversion type: Must be int." + \
                    "Val: {}".format(conversion_type))

    if conversion_type == 0:
        # Going from poolfile to Sequence Set
        
        if run_id is None:
            run_id = "MapTnSeq_Barcodes_run_on_" + \
                        str(datetime.now()).replace(' ', '_'),

        # output dict
        sequence_set = {
                "sequence_set_id": run_id,
                "description": "MapTnSeq (RBTNSEQ) mapping of barcodes to a " + \
                        "genome. Explanations of values: 'nTot' is total " + \
                        "number of times this barcode was counted." + \
                        " 'n' is number of times this barcode was counted" + \
                        " at this location. 'scf' is scaffold name." + \
                        " 'strand' is strand (+ or -). 'pos' is nucleotide" + \
                        " position. 'n2' is number of times found at second" + \
                        " highest counted location. 'scf2' is second highest" + \
                        " location scaffold, 'strand2' similarly, etc." + \
                        " 'nPastEnd' means number of times this barcode was" + \
                        " found next to the next part of the plasmid (barcode" + \
                        " wasn't inserted into the genome without the rest " + \
                        " of the plasmid).\n" + \
                        " User Description (optional): {}".format(description)
        }

        sequence_list = []

        pool_FH = open(inp_fp_path, "r")

        header = pool_FH.readline().rstrip()

        c_line = pool_FH.readline()
        i = 1
        while c_line != "":

            c_lst = c_line.split('\t')
            nPastEnd = c_lst[-1].rstrip()
            barcode, rcbarcode, nTot, n, scf, strnd, pos = c_lst[:7]
            n2, scf2, strnd2, pos2 = c_lst[7:-1]

            # desc_str holds all the information needed to reconstruct pool file
            desc_str = "nTot:{};n:{};scf:{};strand:{};pos:{};".format(
                            nTot, n, scf, strnd, pos)
            desc_str += "n2:{};scf2:{};strnd2:{};pos2:{};".format(
                    n2, scf2, strnd2, pos2)
            desc_str += "nPastEnd:" + nPastEnd

            sequence_list.append(
                    {
                        "sequence_id": "MapTnSeq_barcode_" + str(i),
                        "description": desc_str,
                        "sequence": barcode
            })
            c_line = pool_FH.readline()
            i += 1



        pool_FH.close()

        sequence_set["sequences"] = sequence_list

        with open(op_path, "w") as g:
            g.write(json.dumps(sequence_set, indent=2))

        logging.info("Wrote Sequence Set JSON file to " + op_path)

    elif conversion_type == 1:
        # Going from Sequence Set to Pool File

        if inp_fp_path.split(".")[-1] != "json":
            raise Exception("Sequence Set indicated but not JSON file")

        sequence_set_d = json.loads(open(inp_fp_path).read())


        out_pool_FH = open(op_path, "w")

        out_pool_FH.write("barcode\trcbarcode\tnTot\tn\tscaffold\tstrand\tpos\t" + \
           "n2\tscaffold2\tstrand2\tpos2\tnPastEnd\n")

        seq_list = sequence_set_d["sequences"]

        for seq in seq_list:
            desc_str = seq["description"]
            barcode = seq["sequence"]
            tsl_d = {"A":"T", "T":"A", "G":"C", "C":"G"}
            rcbc1 = [tsl_d[x] for x in list(barcode)]
            rcbc1.reverse()
            rcbarcode = "".join(rcbc1)
            out_pool_FH.write(barcode + "\t" + rcbarcode + "\t")
            items = [x.split(":")[1] for x in desc_str.split(";")]
            out_pool_FH.write("\t".join(items) + "\n")

        out_pool_FH.close()

        logging.info("Wrote Pool File from Sequence Set at " + op_path)


    return None


