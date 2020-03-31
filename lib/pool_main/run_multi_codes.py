#python3


#We pipe in a fastq file and give the following args:
# MultiCodes.pl -minQuality m -n25 -index ind -out outfp
# Questions: is it n25 or bs3?

import subprocess
import gzip
import shutil
import sys,os
sys.path.append(os.path.dirname(__file__))
from MultiCodes import multi_codes_imported_run

# Deciding Index:
# Run Barseqlocal tries to figure out index from the filename
# Require that it has IT(ddd) within the name - RunBarSeqLocal automatically assigns index to file.
# 
# !RunBarSeqLocal is a wrapper for MultiCodes and CombineBarSeq
# RunBarSeqLocal in a directory and a pool file already produced used design random file.
# MultiCodes
# CombineBarSeq is last program to run.
# 2 files: A huge table with one row per barcode, col per experiement
# A table for how many reads were there for each sample and did they match the pool.
# File that contains 
# 3 programs: Last one takes metadata and produces website from tables and metadata.
#
#


# both inputs str paths
def unzip_fastq(zipped_fastq_fp, fastq_dest_path):
    with gzip.open(zipped_fastq_fp, 'rb') as f_in:
        with open(fastq_dest_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)




def run_multi_codes_from_dict(mult_cod_cfg_dict):
    out_prefix = mult_cod_cfg_dict['out_fp_prefix']
    commands = [
            '-out',
            out_prefix,
            '-fastq_fp',
            mult_cod_cfg_dict['fastq_fp'],
            '-index',
            mult_cod_cfg_dict['index'],
            "-n25"
            ]

    if mult_cod_cfg_dict["debug"]:
        commands.append("-debug")

    #Note: first arg in commands is just a placeholder.
    vars_dict = multi_codes_imported_run(commands)

    out_dict = {
        "codes_file": out_prefix + ".codes",
        "counts_file": out_prefix + ".counts",
        "close_file": out_prefix + ".close",
        "mc_report_dict": vars_dict["report_dict"],
        "mc_report_str": vars_dict["report_str"]
            }

    return out_dict 



def test_run_mult_cod():
    mult_cod_cfg_dict = create_test_vars()
    run_multi_codes_from_dict(mult_cod_cfg_dict)
    
    return 0


def create_test_vars():
    mult_cod_cfg_dict = {
            "multi_codes_fp" : "",
            'min_quality' : "0",
            'index' : "IT001",
            'out_fp_prefix' : "Test_Files/tmp/test_out_1",
            'n25_or_bs3_val' :  "-n25",
            'fastq_fp' : 'Test_Files/FEBA_BS_152_IT001_S97_L006_R1_001.fastq',
            "log_fp" : "Test_Files/tmp/test_log_1",
            "debug": True
            }

    return mult_cod_cfg_dict


def test():
    test_run_mult_cod()


def main():
    test()



if __name__ == "__main__":
    main()
