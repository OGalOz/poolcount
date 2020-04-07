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
