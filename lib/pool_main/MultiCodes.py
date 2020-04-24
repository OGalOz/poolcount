#'indexseqi This file is a translation of Morgan Price's MultiCodes.pl into python.
# Note re.search only works on strings
import json
import os
import sys
import math
import argparse
import statistics
import logging
import re
import copy



#multicodes_args is a list, written as actual arguments besides the first index,
#   which is normally the name of this file.
def multi_codes_imported_run(multicodes_args):

    info_dict = load_config()

    info_dict["current_dir"] = os.path.dirname(os.path.abspath(__file__))

    vars_dict = info_dict["vars"]["multi_codes_start_vars"] 

    usage = get_usage_str(info_dict)
    vars_dict["usage"] = usage

    #Parsing input arguments here:
    vars_dict = init_options(vars_dict, multicodes_args)     
   
    
    vars_dict = test_and_update_options(vars_dict, info_dict) #n25 or bs3 here.
    vars_dict = get_nPreExpectedMinMax(vars_dict)
    vars_dict = get_index_two_len(vars_dict)
    vars_dict = parse_index_file_or_iname(vars_dict) 
            
    vars_dict = parse_fastq_input(vars_dict) 
    vars_dict = print_report(vars_dict)
    vars_dict = write_to_out_codes_and_get_nPerCount(vars_dict)
    vars_dict = write_to_out_counts(vars_dict)
    vars_dict = off_by_one_cases_and_write_to_out_close(vars_dict)
    vars_dict = estimate_diversity(vars_dict)
    vars_dict = estimate_bias(vars_dict)

    return vars_dict


def main_run():
    #NOT RUN IN KBASE - ONLY STAND ALONE
    logging.basicConfig(level=logging.DEBUG)

    info_dict = load_config()


    vars_dict = info_dict["vars"]["multi_codes_start_vars"] 

    usage = get_usage_str(info_dict)
    vars_dict["usage"] = usage

    #Parsing input arguments here:
    vars_dict = init_options(vars_dict, None)     
   
    
    vars_dict = test_and_update_options(vars_dict, info_dict) #n25 or bs3 here.
    vars_dict = get_nPreExpectedMinMax(vars_dict)
    vars_dict = get_index_two_len(vars_dict)
    vars_dict = parse_index_file_or_iname(vars_dict) 
            
    vars_dict = parse_fastq_input(vars_dict) 
    vars_dict = print_report(vars_dict)
    vars_dict = write_to_out_codes_and_get_nPerCount(vars_dict)
    vars_dict = write_to_out_counts(vars_dict)
    vars_dict = off_by_one_cases_and_write_to_out_close(vars_dict)
    vars_dict = estimate_diversity(vars_dict)
    vars_dict = estimate_bias(vars_dict)

    return 0


def load_config():
    #The config has base var inputs:
    info_config_fp = os.path.join(os.path.dirname(os.path.realpath(__file__)),"info.json")
    with open(info_config_fp, "r") as f:
        info_dict = json.loads(f.read())
    return info_dict



def estimate_bias(vars_dict):

    if vars_dict['minQuality'] > 0 and vars_dict['nUniq'] >= 5000:
        #What fraction of reads are accounted for by the top 1% of strains?
        f = 0.01
        nReadsSofar = 0
        nCodesSofar = 0
        countSofar = -1
        nPerCount = vars_dict['nPerCount']
        sorted_npc_keys = sorted(nPerCount.keys())
        for count in sorted_npc_keys:
            nReadsSofar += count * nPerCount[count]
            nCodesSofar += nPerCount[count]
            countSofar = count
            if nCodesSofar >= f * vars_dict['nUniq']:
                break
        vars_dict["report_dict"]["estimate_bias"] = {
            "countSofar": countSofar,
            'percent_codes': (100 * nCodesSofar)/vars_dict['nUniq'],
            'k_codes': nCodesSofar/1000,
            'percent_reads': (100 * nReadsSofar)/vars_dict['nOff'][20],
            'k_reads': nReadsSofar/1000
                }
        report_str = """Barcodes with >= {} reads each: 
                {}% of codes ({} K)
                {}% of reads ({} K)\n""".format(
                   countSofar,
                   (100 * nCodesSofar)/vars_dict['nUniq'],
                   nCodesSofar/1000,
                   (100 * nReadsSofar)/vars_dict['nOff'][20],
                   nReadsSofar/1000
                    )
        vars_dict["report_str"] += report_str
        logging.info(report_str)

    return vars_dict
    



def estimate_diversity(vars_dict):
    if vars_dict['minQuality'] > 0 and vars_dict['nOff'][20] >= 1000 \
            and vars_dict['nPerCount'][2] >= 10:
        for fNoise in [0, 0.005, 0.01, 0.02]:
            nNoise = int(0.5 + vars_dict['nOff'][20] * fNoise)
            if nNoise > vars_dict['nPerCount'][1]:
                continue

            vars_dict["report_dict"]["estimate_diversity"] = {
        "percent_noise": fNoise * 100,
        "diversity": ((vars_dict['nUniq'] - nNoise + \
                        (vars_dict['nPerCount'][1] - nNoise)**2)/ \
                        (2 * vars_dict['nPerCount'][2]))/1000,
        "seen_once" : (vars_dict['nUniq'] - nNoise)/1000,
        "seen_twice": vars_dict['nPerCount'][2]/1000
                    }

            report_str =  """If {}% of reads are noise: diversity {} K 
            from total barcodes {} K 
            seen once {} K 
            seen twice {} K\n""".format((fNoise * 100),
                        ((vars_dict['nUniq'] - nNoise + \
                        (vars_dict['nPerCount'][1] - nNoise)**2)/ \
                        (2 * vars_dict['nPerCount'][2]))/1000,
                        (vars_dict['nUniq'] - nNoise)/1000,
                        (vars_dict['nPerCount'][1] - nNoise)/1000,
                        vars_dict['nPerCount'][2]/1000)

            vars_dict["report_str"] += report_str
            logging.info(report_str)
    
    if (vars_dict['minQuality']>0 and vars_dict['nOff'][20]> 1000 \
            and vars_dict['doOff1']):
        nGoodCodes = 0
        nGoodReads = 0
        for code in vars_dict['codes'].keys():
            count = vars_dict['codes'][code]
            tot = sum(count)
            if (code in vars_dict['offby1'] != None and tot > 1 ):
                nGoodCodes += 1
                nGoodReads += tot

        vars_dict["report_dict"]["estimate_diversity"]["good"] = {
            "k_good_codes": nGoodCodes/1000,
            "percent_good_reads":nGoodReads * 100/vars_dict['nOff'][20]
                }
        report_str = """Aside from singletons and off-by-1s, 
        see {} K barcodes ({}% of reads)\n""".format(
                    nGoodCodes/1000,
                    nGoodReads * 100/vars_dict['nOff'][20]
                    )
        vars_dict["report_str"] += report_str
        logging.info(report_str)
    return vars_dict


def off_by_one_cases_and_write_to_out_close(vars_dict):
    offby1 = {} # barcode => 1 for likely off-by-1 errors
    #doOff1 is a bool
    if vars_dict['doOff1']:
        out_close_str = "\t".join(['code1', 'count1', 'code2', 'count2']) + "\n"
        nCases = 0
        nOff1Reads = 0
        codes = vars_dict['codes']
        for code in codes.keys():
            count = codes[code]
            variants = get_variants(code) # variants is a list
            for variant in variants:
                if ((code > variant) - (code < variant)) > 0 and \
                        variant in codes:  # first operator a substitute for cmp
                    n1 = sum(count)
                    n2 = sum(codes[variant])
                    out_close_str += "\t".join([str(code),
                        str(n1), str(variant), str(n2)]) + "\n"
                    if n1 < n2:
                        offby1[code] = 1
                        nOff_val = n1
                    else:
                        offby1[variant] = 1
                        nOff_val = n2
                    nCases += 1
                    nOff1Reads += nOff_val

        out_close_fp = vars_dict["out"] + ".close"
        with open(out_close_fp, "w") as f:
            f.write(out_close_str)

        if 20 in vars_dict['nOff']:
            denominator = vars_dict['nOff'][20]
        else:
            denominator = 1
        fOff1 = str(nOff1Reads/denominator)

        vars_dict["report_dict"]["off_by_one"] = {
            "nCases": nCases,
            "nOff1Reads": nOff1Reads,
            "fOff1": fOff1,
            "out_close_fp": out_close_fp
                }
        report_str = """Wrote {} off-by-1 pairs ({} reads, fraction {}) 
                to {}\n""".format(nCases, nOff1Reads, fOff1, out_close_fp)
        vars_dict["report_str"] += report_str
        logging.info(report_str)

        vars_dict['nCases'] = nCases
        vars_dict['nOff1Reads'] = nOff1Reads
        vars_dict['fOff1'] = fOff1

    vars_dict['offby1'] = offby1

    return vars_dict
    






def write_to_out_counts(vars_dict):
    out_counts_str = "\t".join(["Count","nCodes", "Frac"]) + "\n"
    nPerCount_keys = vars_dict['nPerCount'].keys()
    sorted_npc_keys = sorted(nPerCount_keys)
    for count in sorted_npc_keys:
        out_counts_str += "\t".join([str(count),
            str(vars_dict['nPerCount'][count]),
            str((vars_dict['nPerCount'][count]/vars_dict['nUniq']))
            ]) + "\n"
    out_counts_fp = vars_dict["out"] + ".counts"
    with open(out_counts_fp, "w") as f:
        f.write(out_counts_str)

    vars_dict["report_dict"]["out_counts_info"] = {
            "out_counts_fp": out_counts_fp
            }
    report_str = "Wrote the number of barcodes seen a certain number of times" \
                  +  " to {}\n".format(out_counts_fp)
    if 1 in vars_dict['nPerCount']:
        report_str += "nOnce = {}\n".format(vars_dict['nPerCount'][1])
        vars_dict["report_dict"]["out_counts_info"]["nOnce"] = vars_dict[
                "nPerCount"][1]
    vars_dict["report_str"] += report_str
    logging.info(report_str)

    return vars_dict
        
    


def write_to_out_codes_and_get_nPerCount(vars_dict):
    #Will add this to vars_dict later
    nPerCount = {} # number of codes with that count
    header_list = ["barcode"] + vars_dict['prefixNames']
    #initializing codes tmp string
    out_codes_str = "\t".join(header_list) + "\n"

    nPrefix = len(vars_dict['prefix']) # will equal 1 if iname given
    for code in vars_dict['codes'].keys():
        count_list = vars_dict['codes'][code]
        counts = []
        for i in range(nPrefix):
            if  i < len(count_list):
                counts.append(count_list[i])
            else:
                counts.append(0)
        current_row_list = [code] + counts
        out_codes_str += "\t".join([str(x) for x in current_row_list]) + "\n"
        nAll = sum(counts)
        if nAll in nPerCount:
            nPerCount[nAll] += 1
        else:
            nPerCount[nAll] = 1
    out_codes_fp = vars_dict['out'] + ".codes"
    with open(out_codes_fp, "w") as f:
        f.write(out_codes_str)

    vars_dict["report_dict"]["out_codes_info"] = {
        "unique_barcodes": len(vars_dict['codes'].keys())
            }
    report_str = "Wrote {} unique barcodes to {}".format(
        len(vars_dict['codes'].keys()),
        out_codes_fp
        )
    vars_dict["report_str"] += report_str
    logging.info(report_str)
    vars_dict['nPerCount'] = nPerCount
    return vars_dict

# Adds nUniq to vars_dict
def print_report(vars_dict):
    nOffTot = sum(vars_dict['nOff'].values())
    vars_dict['nUniq'] = len(vars_dict['codes'].keys())
    if vars_dict['nReads'] > 0:

        vars_dict["report_dict"]["print_report"] = {
                "base" : {
                "reads": vars_dict['nReads'],
                "multiplexed": vars_dict['nMulti'],
                "fraction": vars_dict['nOff'][20],
                "unique_codes": vars_dict['nUniq']
                }
        }
        report_str = """Reads {} Multiplexed {} Usable(20) {} 
        fraction unique codes  {} \n""".format(
                vars_dict['nReads'],
                vars_dict['nMulti'],
                vars_dict['nOff'][20],
                vars_dict['nUniq']
            )
        vars_dict["report_str"] += report_str
        logging.info(report_str)

    if vars_dict['nReads'] > 0 and "index2" in vars_dict:

        vars_dict["report_dict"]["print_report"]["index2"] = {
                "n_wrong_index2": vars_dict['nWrongIndex2'],
                "reads_percent": (100*vars_dict['nWrongIndex2']/vars_dict[
                    'nReads'])
                }
        report_str = "Failed to match index2 : {} reads ({}%)\n".format(
            vars_dict['nWrongIndex2'],
            (100*vars_dict['nWrongIndex2']/vars_dict['nReads'])
            )
        vars_dict["report_str"] += report_str
        logging.info(report_str)
    if vars_dict['nReads'] > vars_dict['nWrongIndex2']:

        vars_dict["report_dict"]["print_report"]["wrong_presequence"] = {
            "wrong_preseq_pos" : vars_dict['nWrongPrePos'],
            "reads_num": (100*vars_dict['nWrongPrePos']/(
               vars_dict['nReads'] - vars_dict['nWrongIndex2'] ))
                }
        report_str = 'Wrong presequence position: {} reads ({}).\n'.format(
            vars_dict['nWrongPrePos'],
            (100*vars_dict['nWrongPrePos']/(
               vars_dict['nReads'] - vars_dict['nWrongIndex2'] ))
            )
        vars_dict["report_str"] += report_str
        logging.info(report_str)
    for off in range(18,23):
        if nOffTot > 0:

            vars_dict["report_dict"]["print_report"]["n_off_tot"] = {
                "off": off,
                "nOff": vars_dict['nOff'][off],
                "ratio": (vars_dict['nOff'][off]/nOffTot)
                    }
            report_str = "Off\t{}\t{}\t{}\n".format(
                off,
                vars_dict['nOff'][off],
                (vars_dict['nOff'][off]/nOffTot)
                )
            vars_dict["report_str"] += report_str
            logging.info(report_str)

    return vars_dict


def parse_fastq_input(vars_dict):
    vars_dict['nWrongPrePos'] = 0
    vars_dict['nWrongIndex2'] = 0

    input_fastq_fp = vars_dict['fastq_fp']
    with open(input_fastq_fp, "r") as f:
        fastq_str = f.read()

    fastq_lines = fastq_str.split('\n')


    for j in range(int(math.floor((len(fastq_lines)/4)))):
        i = j*4
        vars_dict['nReads'] += 1
        seq = (fastq_lines[i+1]).rstrip()
        quality = fastq_lines[i+3] # Why no rstrip here?
        if vars_dict['nLimit'] != None and (
                vars_dict['nReads'] >= vars_dict['nLimit']):
            break

        if "index2" in vars_dict:
            part = seq[vars_dict['index2At'] -1: vars_dict['index2len']]
            if part != vars_dict['index2']:
                vars_dict['nWrongIndex2'] += 1
                continue
        if vars_dict["iname"] != None:
            vars_dict['iPrefix'] = 0
        else:
            vars_dict['iPrefix'] = find_prefix(seq, vars_dict['prefix'], 
                    vars_dict['debug']) # returns match (int) or -1
        if vars_dict['iPrefix'] < 0:
            continue # No prefixes found, skip this iteration of loop.
        vars_dict['nMulti'] += 1

        #In most cases, iname is None, so iPrefix is 0, so we get nLeading = 0 
        nLeading, indexseq, prefixName = vars_dict['prefix'][vars_dict['iPrefix']]
        offset = nLeading + len(indexseq) #offset will likely be 0
        barcode, off = find_barcode(seq, quality, offset, vars_dict) 
        if (barcode == None and off != None and off >=0):
            vars_dict['nWrongPrePos'] += 1
            if vars_dict['nWrongPrePos'] >= 200 and (vars_dict[
                'nWrongPrePos'] >= (0.1 * vars_dict['nReads'])):
                raise Exception("Over 10% of reads have the wrong spacing ( \
                        not {}:{}) to the pre-sequence ({} \
                        of {} so far).\n Perhaps you forgot to specify the \
                        protocol (i.e., -n25 or -bs3)?\n".format(
                            vars_dict["nPreExpectedMin"],
                            vars_dict["nPreExpectedMax"],
                            vars_dict['nWrongPrePos'],
                            vars_dict['nReads']))
        if barcode == None:
            continue
        vars_dict['nOff'][off] += 1
        if off != 20:
            continue
        if (not (barcode in vars_dict['codes'])):
            vars_dict['codes'][barcode] = []

        vars_dict = update_codes_iprefix(vars_dict['codes'], barcode, 
                vars_dict['iPrefix'], vars_dict) # line 226 perl

        #End fastq loop
    
    return vars_dict    


#barcode str, iPrefix int
#All this to substitute the creative properties of Perl
def update_codes_iprefix(codes_dict, barcode, iPrefix, vars_dict):

    if barcode in codes_dict:
        barcode_list = codes_dict[barcode]
        if len(barcode_list) > iPrefix:
            barcode_list[iPrefix] += 1
        else:
            dif = (iPrefix - len(barcode_list)) + 1
            barcode_list += [0]*dif
            barcode_list[iPrefix] += 1
    else:
        barcode_list = [0] * (len(iPrefix) + 1) 
        barcode_list[iPrefix] += 1

    codes_dict[barcode] = barcode_list
    vars_dict['codes'] = codes_dict 
    return vars_dict
        





        



#indexfile must be None or a true path
def parse_index_file_or_iname(vars_dict):
    #New variables:
    vars_dict['nReads'] = 0
    vars_dict['nMulti'] = 0 #number with prefix identified
    nOff = {} # only spacings of 20 are considered, count totals
    for i in range(18,23):
        nOff[i] = 0
    vars_dict['nOff'] = nOff
    vars_dict['codes'] = {} #barcode maps to number of times barcode is seen,
    # with an array of one entry per multiplex tag

    prefix = [] # list of [nLeading, indexseq (no leading Ns), name]
    
    indexfile = vars_dict['indexfile']
    if indexfile != None:
        with open(indexfile, "r") as f:
            indexfile_str = f.read()
        indexfile_lines = indexfile_str.split("\n")
        for line in indexfile_lines:
            line.rstrip()
            name, index = line.split('\t')
            if (re.search(r'name', name ,re.IGNORECASE)):
                    continue
            nLeading = None
            indexseq = None
            if (name == "none" and index == ""):
                nLeading = 0
                indexseq = ""
            else:
                match = re.search(r'^([nN]*)([ACGT]+)$',index)
                if not match:
                    raise Exception("Invalid index sequence {}".format(index))
                else:
                    nLeading = len(match[0])
                    indexseq = match[1]
            if (nLeading == None ) or (indexseq == None) or (name == ''):
                raise Exception(line)
            prefix.append([nLeading, indexseq, name])
            report_str = "Read {} indices from {}\n".format(
                len(prefix),indexfile)
            vars_dict["report_str"] += report_str
            logging.info(report_str)
            prefixNames = [x[2] for x in prefix]
    elif (vars_dict['iname'] != None):
        #Usually goes here
        prefix = [[0, "", vars_dict['iname']]]
        prefixNames = [vars_dict['iname']]
    else:
        debug_varprint(vars_dict)
        raise Exception("Tested previously in program, indexfile or iname " \
                + "must be present.")

    #updating vars_dict
    vars_dict['prefix'] = prefix
    vars_dict['prefixNames'] = prefixNames

    return vars_dict


def get_usage_str(info_dict):

    usage_fp = os.path.join(info_dict["current_dir"],
            info_dict["texts"]["multi_codes_usage_fp"])
    with open(usage_fp, "r") as f:
        usage = f.read()
    return usage


def init_options(vars_dict, args_list):



    #Necessary args: -out, -fastq_fp, -index/ -primers 
    parser = argparse.ArgumentParser(description=vars_dict['usage'])

    parser.add_argument("-primers", type=str) # goes to indexfile - this xor 
    # -index, not both
    parser.add_argument("-out", type=str)
    parser.add_argument("-index", type=str) # goes to iname. This xor -primers.
    # index could be: 
    parser.add_argument("-fastq_fp", type=str) # This is different from usage.
    parser.add_argument("-minQuality", type=int) #Goes to minQuality
    parser.add_argument("-limit", type=int) # Goes to nLimit
    parser.add_argument("-nPreExpected", type=str) # Goes to nPreExpected
    parser.add_argument("-preseq", type=str) # Goes to preseq
    parser.add_argument("-postseq", type=str) # Goes to postseq
    parser.add_argument("-off1", type=int) # Goes to doOff1
    #FLAGS
    #if dntag present its value will be 1
    parser.add_argument("-dntag", action="store_const", const='1', default='0')
    parser.add_argument("-debug", action="store_const", const='1', default='0')
    parser.add_argument("-n25", action="store_const", const='1', default='0')
    parser.add_argument("-bs3", action="store_const", const='1', default='0')
    
    #PARSING COMMAND LINE ARGS:
    if args_list == None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args_list)

    #Updating vars_dict:
    vars_dict["indexfile"] = args.primers
    vars_dict["out"] = args.out
    vars_dict["iname"] = args.index 
    vars_dict["fastq_fp"] = args.fastq_fp 
    if args.minQuality != None:
        #Normally we just use minQuality from the config
        vars_dict["minQuality"] = args.minQuality 
    vars_dict["nLimit"] =  args.limit 
    vars_dict["nPreExpected"] = args.nPreExpected  
    vars_dict["preseq"] =  args.preseq 
    vars_dict["postseq"] = args.postseq 
    vars_dict["doOff1"] = args.off1 
    
    num_to_truth = {'0': False, '1':True}
    vars_dict["dntag"] = num_to_truth[args.dntag] 
    vars_dict["debug"] = num_to_truth[args.debug] 
    vars_dict["n25"] =  num_to_truth[args.n25]
    vars_dict["bs3"] =  num_to_truth[args.bs3]



    if (vars_dict['indexfile'] !=None and vars_dict['iname'] != None):
        raise Exception("There can be only one of primers (indexfile) or " \
                + "index (iname). Please see usage: {}".format(vars_dict["usage"]))

    return vars_dict

# This function runs through variables and checks if program can run
def test_and_update_options(vars_dict, info_dict):
    """
    Defining preseq if not defined.
    

    """
    if vars_dict['doOff1'] == None:
        if vars_dict['minQuality']:
            vars_dict['doOff1'] = True
        else:
            #Usually go here
            vars_dict['doOff1'] = False

    if vars_dict['preseq'] != None:
        if vars_dict['postseq'] == None:
            raise Exception("Missing -postseq: {}".format(vars_dict['usage']))
        if vars_dict['nPreExpected'] == None:
            raise Exception("Missing -nPreExpected: {}".format(vars_dict['usage']))
    else:
        #Usually go here
        seq_dict = info_dict["sequences"]["MultiCodes"]
        if vars_dict['dntag']:
            if vars_dict['iname'] != None:
                raise Exception("-index with -dntag not supported")
            else:
                vars_dict['preseq'] = seq_dict["dntag_preseq"]
                if vars_dict['nPreExpected'] == None:
                    vars_dict['nPreExpected'] = 8
                vars_dict['postseq'] = seq_dict["dntag_postseq"]
        else:
            #Usually go here
            vars_dict['preseq'] = seq_dict["base_preseq"]
            if vars_dict['nPreExpected'] == None:
                if vars_dict['iname'] != None:
                    #Usually here
                    vars_dict['nPreExpected'] = 14
                else:
                    vars_dict['nPreExpected'] = 9
            vars_dict['postseq'] = seq_dict["base_postseq"]

    vars_dict = test_n_25_bs3(vars_dict)
    return vars_dict


def get_index_two_len(vars_dict):
    index2len = None
    if "index2" in vars_dict:
        index2 = vars_dict["index2"]
        match = re.search(r'^[ACGT]+$', index2)
        if not match:
            raise Exception("Invalid index2 sequence {}.".format(index2))
        index2At = vars_dict['index2At']
        match = re.search(r'^\d+$', index2At)
        if not (match and index2At >= 1):
            raise Exception("Invalid index2At: {}.".format(index2At))
        index2len = len(index2)
        report_str = "Checking for index2 {} at position \
                {}\n".format(index2, index2At)
        vars_dict["report_str"] += report_str
        logging.info(report_str)
    vars_dict['index2len'] = index2len
    return vars_dict



def test_n_25_bs3(vars_dict):

    if vars_dict['n25'] and vars_dict['bs3']:
        raise Exception("Cannot specify multiple protocols. Both bs3 and \
                n25 defined.\n")
    if vars_dict['n25']:
        #usually go here
        if vars_dict['iname'] == None:
            raise Exception("Cannot specify -n25 unless -index is set\n")
        #Overwrite nPreExpected
        vars_dict['nPreExpected'] = "11:14"

    if vars_dict['bs3']:
        if vars_dict['iname'] == None:
            raise Exception("Cannot specify -bs3 unless -index is set\n")
        vars_dict['ifile'] = vars_dict['bin_dir'] + "/../primers/barseq3.index2"
        if not os.path.isfile(vars_dict['ifile']):
            raise Exception("No such file: {}".format(vars_dict['ifile']))
        #This is a special function to get the table fom file somehow
        tab_base = ["index_name", "index2", "nN"]
        #tab is a list of dicts with each line from ifile
        tab = read_table(vars_dict['ifile'], tab_base)
        #tabMatch is a subset of tab with dicts whose index name matches iname
        tabMatch =  mygrep(tab, "index_name", vars_dict['iname']) 
        if len(tabMatch) == 0: #?
            
            warning_str = "Warning! Ignoring the second index -- index_name " \
                 + "{} does not appear in {}\n".format(vars_dict['iname'], 
                        vars_dict['ifile'])
            vars_dict["report_dict"]["warnings"].append(warning_str)
            logging.warning(warning_str)
            vars_dict['nPreExpected'] = '16:19'
        else:
            if len(tabMatch) > 1: #?
                raise Exception("Multiple matches for index {} in {}".format(
                    vars_dict['iname'], vars_dict['ifile']))
            #length of tabMatch is exactly 1 
            row = tabMatch[0]
            vars_dict['index2'] = row['index2']
            vars_dict['index2At'] = row['nN'] + 1 
            vars_dict['nPreExpected'] = vars_dict['index2At'] + 6 + 9
    return vars_dict

#ifile is filepath, required is a list of required headers
# returns a list of hashmaps 
def read_table(ifile, required):
    with open(ifile, "r") as f:
        f_str = f.read()
    f_list = f_str.split("\n")
    
    header_line = f_list[0].rstrip()
    header_list = header_line.split("\t")
    cols_map = {}
    for i in range(header_list):
        cols_map[header_list[i]] = i
    for field in tab_base:
        if field not in cols_map:
            raise Exception("No field {} in {}".format(
                field, ifile))
    rows = []
    for line in f_list[1:]:
        line = line.rstrip()
        new_line_list = line.split("\t")
        if len(new_line_list) != len(header_list):
            raise Exception("Wrong number of columns in:\n{} \
                    \n in {}".format(line, ifile))
            row_dict = {}
            for i in range(len(new_line_list)):
                row_dict[header_list[i]] = new_line_list[i]
            rows.append(copy.deepcopy(row_dict))

    return rows


def mygrep(hash_list, index_name, iname):
    new_list = []
    for h in hash_list:
        if h[index_name] == iname:
            new_list.append(h)
    return new_list


def get_nPreExpectedMinMax(vars_dict):
    match = re.search(r'^(\d+):(\d+)$',str(vars_dict['nPreExpected']))
    if match:
        npre_list = str(vars_dict['nPreExpected']).split(":")
        vars_dict['nPreExpectedMin'] = int(npre_list[0])
        vars_dict['nPreExpectedMax'] = int(npre_list[1])
    else:
        #Instead of being m:n, nPreExpected is an int (in string format)
        try:
            npE = int(vars_dict["nPreExpected"])
        except:
            raise Exception(vars_dict['usage'])
        vars_dict['nPreExpectedMin'] = npE - 2
        vars_dict['nPreExpectedMax'] = npE + 2
   
    return vars_dict


#seq is DNA sequence (str) 
#indexes is a list of lists, internal lists are [nLeading, indexseq, name],
#   where nLeading is an int, indexseq is str, name is str 
#debug in vars_dict, True or False 
def find_prefix(seq, indexes, debug):
    matches = []
    for i in range(len(indexes)):
        nLeading, indexseq, name = indexes[i]
        if (seq[nLeading:(nLeading + len(indexseq))] == indexseq):
            matches.append(i)
    if len(matches) == 1:
        return matches[0]
    if debug:
        logging.info("No prefix for {}\n".format(seq))
        return -1

#seq str
#quality str
#offset int, usually 0
def find_barcode(seq, quality, offset, vars_dict):

    """
    Requires more work
    vars_dict keys needed:
        preseq, nPreExpectedMin, nPreExpectedMax, minQuality, debug

    preseq often: 'CAGCGTACG'
    postseq often: 'AGAGACCTC'

    How function works:
        We get true sequence and quality by taking original sequence and 
        starting from offset. Then we find the index of presequence defined 
        earlier in the program within the true sequence and assign it's value
        to preseq_pos. preseq_pos must be within the nPreExpectedMin and
        nPreExpectedMax index values of the true sequence, or it's not viable
        and we return [None, None].
        Then we assign barcode to where the presequence ends to 20 indeces
        after that. If the presequence ends too late or the entire sequence
        is too short and we don't get a full barcode, then we return 
        [None, None]. 
        Then we check for the postsequence in the true sequence. If the length 
        of the true sequence is less than the sum of the barcode and the
        presequence and the index of the presequence and the postsequence, we 
        try a shorter postsequence of length 4, if that's not found we
        claim it's too short, make postsequence "" and continue with function.
        foundEnd: foundEnd is the distance between the end of the presequence 
        and the beginning of the postsequence if it's found. If the 
        postsequence is nothing then foundEnd is automatically 20, 
        if it's not nothing and still not found, then function returns 
        [None, None].


    """

    seq2 = seq[offset:]
    quality2 = quality[offset:]
    preseq = vars_dict['preseq']
    preseq_pos = seq2.find(preseq) # (finding 1st occurence of preseq in seq2) 

    if not (preseq_pos >=0 and preseq_pos >= vars_dict['nPreExpectedMin'] and \
            preseq_pos <= vars_dict['nPreExpectedMax']):
        if vars_dict['debug']:
            warning_str = "seq2 {} has invalid index-of-preseq {}\n".format(
                seq2, preseq_pos)
            vars_dict["report_dict"]["warnings"].append(warning_str)
            logging.warning(warning_str)
        if preseq_pos >= 0:
            return [None, preseq_pos] # report that spacing was wrong
        return [None, None]
    
    #Note: Below line does not throw error if end index > len(seq2)
    barcode = seq2[preseq_pos + len(preseq): preseq_pos + len(preseq) + 20]
    if not (len(barcode) == 20 and re.search(r'^[A-Z]+$',barcode)):
        #note barcodes with ambiguous sequences are ignored
        # we might want to allow one N in the future
        if vars_dict['debug']:
            warning_str = "seq2 {} has invalid barcode sequence {}\n".format(
                seq2, barcode)
            vars_dict["report_dict"]["warnings"].append(warning_str)
            logging.warning(warning_str)
            return [None, None]
    
    postseqUsed = vars_dict['postseq']
    if vars_dict['minQuality'] == 0:
        postseqUsed = ""


    if len(seq2) < preseq_pos + len(preseq) + 20 + len(postseqUsed):
        postseqUsed = postseqUsed[:4] #first 4 characters.
        if vars_dict['debug']:
            warning_str = "Using postseq {}\n".format(postseqUsed)
            vars_dict["report_dict"]["warnings"].append(warning_str)
            logging.warning(warning_str)
        #There might be redundant code here in MultiCodes.pl which checks
        # above length (16 lines above) Redundant code left out of this file.
        #We check again with a shorter postseqUsed 
        if len(seq2) < preseq_pos + len(preseq) + 20 + len(postseqUsed):
            if vars_dict['debug']:
                warning_str = "Ignoring postseq, too short\n"
                vars_dict["report_dict"]["warnings"].append(warning_str)
                logging.warning(warning_str)
            postseqUsed = ""

    foundEnd = -1
    if postseqUsed == "":
        foundEnd = 20
    else:
        #check 20 first in case end of barcode matches post-seq (AGAG issue)
        for off in [20,19,21,18,22]:
            start_slice = preseq_pos + len(preseq) + off
            end_slice = start_slice + len(postseqUsed)
            if (seq2[start_slice:end_slice] == postseqUsed):
                foundEnd = off
                break
    if foundEnd == -1:
        if vars_dict['debug']:
            warning_str = "seq2 {} barcode {} postseq not found\n".format(
                seq2, barcode)
            vars_dict["report_dict"]["warnings"].append(warning_str)
            logging.warning(warning_str)
        return [None,None]

    if vars_dict['minQuality'] > 0:
        # the sanger code for !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ  
        barqual = quality2[preseq_pos + len(preseq): preseq_pos + len(preseq) \
                + 20]

        # quality score is encoded as 33+
        scores = [ord(c) - 33 for c in barqual]
        for score in scores:
            if (score < 0 or score > 100):
                raise Exception("Invalid score {} from barcode {} " \
                       + "quality {}".format(score, barcode, barqual))
            if score < vars_dict['minQuality']:
                if vars_dict['debug']:
                    warning_str = "Low quality {} for barcode {} in " \
                         + "{}\n".format(score, barcode, seq2)
                    vars_dict["report_dict"]["warnings"].append(warning_str)
                    logging.warning(warning_str)
                return [None,None]
    return [barcode, foundEnd]

#baseseq (str)
# returns all variants on sequence
def get_variants(baseseq):
    out = []
    baseseq = baseseq.upper()
    for i in range(len(baseseq)):
        pre = baseseq[0:i]
        char = baseseq[i:i+1]
        post = baseseq[i+1:]
        if not (char in ["A","C","G","T"]):
            continue
        for newchar in ["A","C","G","T"]:
            if not newchar == char:
                out.append(pre + newchar + post)
    return out


def debug_varprint(vars_dict):
    no_usage_dict = vars_dict.deepcopy()
    del no_usage_dict["usage"]
    debug_str = ""
    for key in no_usage_dict.keys():
        debug_str += "{}: {}\n".format(key, str(no_usage_dict[key]))

    logging.info(debug_str)
        
    

def unit_test():
    
    #Testing various arguments:
    args_options = [
    "-out mytest -index myindex -fastq myfastq",
    "",
    ""
    ]
    main_run()

    return 0

def main():
    unit_test()
    #main_run()

if __name__ == "__main__":
    main()
        

def test_find_barcode():
    """
    For Find Barcode func, we need the following inputs:
    vars_dict keys needed:
        preseq, postseq nPreExpectedMin, nPreExpectedMax, minQuality, debug
    Non- vars_dict:
    seq, quality, offset

    Incomplete test
    """
    base_vars_dict_1 = {
            "debug": True,
            "nPreExpectedMin": 11,
            "nPreExpectedMax": 14,
            "minQuality": 10,
            "preseq": "CAGCGTACG"
            } 


    inputs = [{"seq": "",
    "quality": "",
    "offset": "",
    "vars_dict": base_vars_dict_1
    }
            ]
