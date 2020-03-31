#python3
# This program is a translation of Morgan Price's combineBarSeq.pl to 
#   python. Important inputs to this program are the pool_file and the 
#   codesfiles. Outputs are poolcount and colsum
# The oneperfile code is not clean
# "perl will fill in data structures whenever possible"

import os
import sys
import math
import argparse
import statistics
import logging
import re
import json
from functools import cmp_to_key


def combine_barseq_imported_run(args_list):
    info_config_fp = os.path.join(os.path.dirname(os.path.abspath(__file__)),
            "info.json")
    with open(info_config_fp, "r") as f:
        info_dict = json.loads(f.read())
    info_dict["current_dir"] = os.path.dirname(os.path.abspath(__file__))
    usage = get_usage_str(info_dict)
    vars_dict = info_dict["vars"]["combine_barseq_start_vars"]
    vars_dict["usage"] = usage
    vars_dict = init_options(vars_dict, args_list)

    pool = init_pool_dict(vars_dict)

    counts = {} # rcbarcode to vector of counts

    for code_fp in vars_dict['codesFiles']:
        vars_dict = process_codefile(code_fp, vars_dict, pool, counts)


    vars_dict["report_dict"]["pool_info"] = {
            "entries": len(pool.keys()),
            "lines_seen": vars_dict['nUsed'],
            "lines_ignored": vars_dict['nIgnore'],
            "num_files": len(vars_dict['codesFiles'])
            }

    report_str = """Pool {} entries {}.
            Saw {} lines. Ignored {} lines from {} files.\n""".format(
                vars_dict['poolfile'], len(pool.keys()),
                vars_dict['nUsed'], vars_dict['nIgnore'], 
                len(vars_dict['codesFiles']) )
    vars_dict["report_str"] += report_str
    logging.info(report_str)

    vars_dict = categorize_reads(vars_dict)

    vars_dict = write_to_poolcount(pool, counts, vars_dict['nSamples'], 
            vars_dict['out_name'], vars_dict['indexes'],
            vars_dict)

    vars_dict = write_to_colsum(vars_dict['indexes'], vars_dict['colSums'], 
            vars_dict['colSumsUsed'], vars_dict['out_name'],
            vars_dict)


    return vars_dict




def main_run():
    logging.basicConfig(level=logging.INFO)
    info_config_fp = "info.json"
    with open(info_config_fp, "r") as f:
        info_dict = json.loads(f.read())
    usage = get_usage_str(info_dict)
    vars_dict = info_dict["vars"]["combine_barseq_start_vars"]
    vars_dict["usage"] = usage
    vars_dict = init_options(vars_dict, None)

    pool = init_pool_dict(vars_dict)

    counts = {} # rcbarcode to vector of counts

    for code_fp in vars_dict['codesFiles']:
        vars_dict = process_codefile(code_fp, vars_dict, pool, counts)

    report_str = """Pool {} entries {}.\n
            Saw {} lines.\n Ignored {} lines from {} files.\n""".format(
                vars_dict['poolfile'], len(pool.keys()),
                vars_dict['nUsed'], vars_dict['nIgnore'], 
                len(vars_dict['codesFiles']) )
    vars_dict["report_str"] += report_str
    logging.info(report_str)

    #if save ignore write to out.codes.ignored

    vars_dict = categorize_reads(vars_dict)


    vars_dict = write_to_poolcount(pool, counts, vars_dict['nSamples'], 
            vars_dict['out_name'], vars_dict['indexes'],
            vars_dict)

    vars_dict = write_to_colsum(vars_dict['indexes'], vars_dict['colSums'], 
            vars_dict['colSumsUsed'], vars_dict['out_name'],
            vars_dict)



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

def init_options(vars_dict, args_list):
    parser = argparse.ArgumentParser(description=vars_dict['usage'])
    out_help = "the name of the output - need not be a full file name"
    parser.add_argument("out_file_name", help=out_help)
    parser.add_argument("pool_file_name", help ="the name of the poolfile" )
    codes_files_help = "the paths to the codes files- separated by commas " \
            + "and no spaces. e.g:\n file1,file2,file3"
    parser.add_argument("code_files_list", help = codes_files_help)

    if args_list == None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args_list)

    vars_dict['out_name'] = args.out_file_name
    poolfile = args.pool_file_name #poolfile should be path to poolfile
    if not os.path.isfile(poolfile):
        raise Exception("poolfile '{}' is not a file.".format(poolfile))
    vars_dict['poolfile'] = poolfile
    codesFiles = args.code_files_list
    codes_files_list = codesFiles.split(",")
    for code_file in codes_files_list:
        if not os.path.isfile(code_file):
            raise Exception("code_file '{}' is not a file.".format(code_file))

    vars_dict['codesFiles'] = codes_files_list



    return vars_dict



"""
Vars in use:
    indexes (array)
    colSums (array)
    colSumsUsed
"""
def categorize_reads(vars_dict):
    lowcountI = [] #indexes with < 200,000 reads
    lowcounts = [] # those actual counts
    fractions = [] # fraction used for indexes with plenty of reads only
    fractionsS = [] #fraction used for succesful indexes (enough reads & f > 0.25)
    lowhitI = [] #indexes with fraction < 0.25
    okI = [] #indexes with enough reads and fraction above 0.25
    totalReads = 0
    for i in range(len(vars_dict['indexes'])):
        totalReads += vars_dict['colSums'][i]
        if (vars_dict['colSums'][i] < 200*1000):
            lowcountI.append(vars_dict['indexes'][i])
            lowcounts.append(vars_dict['colSums'][i])
        else:
            fraction = (vars_dict['colSumsUsed'][i])/(vars_dict['colSums'][i])
            fractions.append(fraction)
            if fraction < 0.25:
                lowhitI.append(vars_dict['indexes'][i])
            else:
                okI.append(vars_dict['indexes'][i])
                fractionsS.append(fraction)

    vars_dict["report_dict"]["categorize_reads_info"] = {
            "indexes": vars_dict["indexes"],
            "success": len(okI),
            "low_count": len(lowcountI),
            "low_fraction": len(lowhitI),
            "total_reads": totalReads
            }
    report_str = """Indexes {}. Success {}. LowCount {}. LowFraction {}.
            Total Reads (millions): {}.\n""".format(
                len(vars_dict['indexes']), len(okI), len(lowcountI),
                len(lowhitI), (totalReads/(10**6))
                )
    vars_dict["report_str"] += report_str
    logging.info(report_str)
    #stats info:
    if len(fractions) > 0:
        LCE_median = statistics.median(fractions)
    else:
        LCE_median = "NaN"
    if len(fractionsS) > 0:
        sccs_median = statistics.median(fractionsS)
    else:
        sccs_median = "NaN"

    
    vars_dict["report_dict"]["categorize_reads_info"]["median_info"] = {
            "median_fraction_no_lowcount": LCE_median,
            "median_for_succes": sccs_median
            }

    report_str = """Median Fraction (LowCount Excluded) {}.
            Median for Success {}.\n""".format(
            LCE_median, sccs_median)

    vars_dict["report_str"] += report_str
    logging.info(report_str)

    if len(okI) > 0:

        vars_dict["report_dict"]["categorize_reads_info"][
                "okI_short_list"] = short_list(okI, vars_dict)

        report_str = " Success {} \n".format(short_list(okI, vars_dict))
        vars_dict["report_str"] += report_str
        logging.info(report_str)
    if len(lowcountI) > 0:
        vars_dict["report_dict"]["categorize_reads_info"][
                "lowcount_median"] = statistics.median(lowcounts)
        vars_dict["report_dict"]["categorize_reads_info"][
                "lowcountI_short_list"] = short_list(lowcountI, vars_dict)

        report_str = "  LowCount (median {}) {}\n".format(
            statistics.median(lowcounts),
            short_list(lowcountI, vars_dict)
            )
        vars_dict["report_str"] += report_str
        logging.info(report_str)    
    
    #categorize reads vars move to total vars 
    vars_dict['fractions'] = fractions,
    vars_dict['fractionsS'] = fractionsS,
    vars_dict['lowcounts'] = lowcounts,
    vars_dict['lowhitI'] = lowhitI,
    vars_dict['okI'] = okI,
    vars_dict['totalReads'] = totalReads,
    vars_dict['lowcountI'] = lowcountI, 
    
    return vars_dict




"""
Used vars:
    cols,
    nSamples,

"""
def process_codefile(code_fp, vars_dict, pool, counts):
    with open(code_fp, "r") as f:
        code_file_str = f.read()
    code_file_lines = code_file_str.split("\n")
    header = code_file_lines[0]
    header.rstrip()
    header_cols = header.split("\t")

    #We perform tests on the header columns
    # vars_dict is only set first time we run this code
    vars_dict = process_header_cols(header_cols, vars_dict)


    #We process the rest of the lines: (Each a TSV string)
    remaining_lines = code_file_lines[1:]
    vars_dict = process_remaining_lines(remaining_lines, pool, 
            vars_dict, counts)

    return vars_dict



def process_body_line(body_cl_list, pool, vars_dict, counts):

    """
    body_cl_list = code_line split - so [rcbarcode (str), Index number (str of int)]
    pool (dict) keys: barcodes, values: (list) [ barcode, scaffold, strand, pos]
    counts: (dict)
    """
    rcbarcode = body_cl_list[0] #reverse complement barcode
    cl_list = body_cl_list[1:] #Often length one - str of int

    #Checking if rcbarcode string:
    if not re.match(r"^[ACGTN]+$", rcbarcode):
        raise Exception("Invalid rcbarcode: {}".format(rcbarcode))

    #We create a value to check correctness of length of cl_list
    if vars_dict['oneperfile']:
        column_testing_value = 1
    else:
        column_testing_value = vars_dict['nSamples']
    if not (len(cl_list) == column_testing_value):
        raise Exception("Wrong number of columns in {}".format(code_fp))


    #CHECK HERE
    if rcbarcode in pool:
        vars_dict['nUsed'] += 1
        if vars_dict['oneperfile']:
            update_colsums_used(vars_dict, cl_list)
            update_counts_dict(counts, rcbarcode, vars_dict['thisIndex'], cl_list)
        else:
            for i in range(len(vars_dict['nSamples'])):
                vars_dict['colSumsUsed'][i] += cl_list[i]
            if rcbarcode in counts:
                #In the perl code- row (below) is a reference to an array @row
                #In this python version, row is an array 
                #length of row is supposed to be the same as nSamples
                row = counts[rcbarcode]
                for i in range(len(vars_dict['nSamples'])):
                    row[i] += cl_list[i]
                #We update the dict with new row (if not already connected)
                counts[rcbarcode] = row
            else:
                counts[rcbarcode] = cl_list
    else:
        if vars_dict["save_ignore"]:
            logging.info("IGNORE {}\n".format("\t".join([rcbarcode] + cl_list)))
        vars_dict['nIgnore'] += 1
    if vars_dict['oneperfile']:
        vars_dict['colSums'][vars_dict['thisIndex']] += int(cl_list[0])
    else:
        for i in range(len(vars_dict['nSamples'])):
            vars_dict['colSums'][i] += cl_list[i]

    return vars_dict

#counts - dict, rcbarcode (str), thisIndex (int), cl_list may be list of
# length 1 with a string of an int. thisIndex refers to index in cl_list
# to which we are adding this value?
def update_counts_dict(counts, rcbarcode, thisIndex , cl_list):
    #cl_list replaces "F" in Perl
    if rcbarcode in counts:
        counts_barcode_list = counts[rcbarcode]
        if len(counts_barcode_list) > thisIndex:
            counts[rcbarcode][thisIndex] += int(cl_list[0])
        else:
            counts_barcode_list += [0] * (thisIndex - len(counts_barcode_list) + 1)
            counts_barcode_list[thisIndex] += int(cl_list[0])
            counts[rcbarcode] = counts_barcode_list
    else:
        counts_barcode_list = [0] * (thisIndex + 1)
        counts_barcode_list[thisIndex] += int(cl_list[0])
        counts[rcbarcode] = counts_barcode_list


def update_colsums_used(vars_dict, cl_list):
    colSumsUsed = vars_dict['colSumsUsed']
    if len(colSumsUsed) > vars_dict['thisIndex']:
        vars_dict['colSumsUsed'][vars_dict['thisIndex']] += int(cl_list[0])
    else:
        colSumsUsed += [0] * (vars_dict['thisIndex'] - len(colSumsUsed) + 1)
        colSumsUsed[vars_dict['thisIndex']] += int(cl_list[0])
        vars_dict['colSumsUsed'] = colSumsUsed

#remaining lines is a list of strings with TSVs
# pool is a dict
def process_remaining_lines(remaining_lines, pool, vars_dict, counts):
    #We count rows
    nThisFile = 0
    for i in range(len(remaining_lines)):
        code_line = remaining_lines[i]
        nThisFile += 1
        
        code_line = code_line.rstrip()
        cl_list = code_line.split("\t")

        if len(cl_list) > 1:
            #This is where the majority of the work is done:
            vars_dict = process_body_line(cl_list, pool, vars_dict, counts)
        else:

            warning_text = "{}: {}".format(i, code_line)
            vars_dict["report_dict"]["warnings"].append(warning_text)
            logging.warning(warning_text)
    if nThisFile == 0:
        warning_text = "No entries in {}\n".format(code_fp)
        vars_dict["report_dict"]["warnings"].append(warning_text)
        logging.warning(warning_text)
         

    return vars_dict


#header_cols is a list of strings, headers of code file (first row)
def process_header_cols(header_cols, vars_dict):
    if len(header_cols) < 2:
        raise Exception("Not enough columns in {}".format(code_fp))
    first = header_cols[0]
    cols = header_cols[1:]
    if first != "barcode":
        raise Exception("Not a barcode counts file")
    if vars_dict['nSamples'] == 0: #first file
        vars_dict['indexes'] = cols 
        vars_dict['nSamples'] = len(cols)
        vars_dict['colSums'] = [0] * vars_dict['nSamples'] #length of colSums is num samples
        vars_dict['colSumsUsed'] = [0] * vars_dict['nSamples'] #length of colSumsUsed is num samples
        if len(cols) == 1:
            vars_dict['oneperfile'] = True
            vars_dict['thisIndex'] = 0
    elif vars_dict['oneperfile']:
        if len(cols) > 1:
            raise Exception("More than one data column in {} despite having \
                    oneperfile = True".format(code_fp))
        #Creating an oldcol dict with index association
        oldcol = {}
        for i in range(len(vars_dict['indexes'])):
            oldcol[vars_dict['indexes'][i]] = i
        vars_dict['thisCol'] = cols[0]
        if vars_dict['thisCol'] in oldcol:
            vars_dict['thisIndex'] = oldcol[vars_dict['thisCol']]
        else:
            vars_dict['nSamples'] += 1
            vars_dict['indexes'].append(vars_dict['thisCol'])
            vars_dict['thisIndex'] = len(vars_dict['indexes']) - 1
            vars_dict['colSums'].append(0)
            vars_dict['colSumsUsed'].append(0)
    else: #We already have all the variables, checking if they align with new file.
        if not (len(cols) == len(vars_dict['indexes'])):
            raise Exception("Wrong number of columns in {}".format(
                vars_dict['code_fp']))
        for i in range(len(cols)):
            if not (cols[i] == vars_dict['indexes'][i]):
                raise Exception("Index mismatch in {} \
                        vs. {} -- {} vs. {} ".format(
                            vars_dict['code_fp'], vars_dict['codesFiles'][0],
                            vars_dict['cols'][i], vars_dict['indexes'][i]))
    return vars_dict


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


#inp_arr is input array - normally okI "ok Indexes",
def short_list(inp_arr, vars_dict):
    my_list = sorted(inp_arr)

    #Check that the list is greater than length 0
    if len(my_list) == 0:
        warning_text = "The current input to short list contains no ok " \
                + "Values."
        vars_dict["report_dict"]["warnings"].append(warning_text)
        return " "

    my_match = re.search(r'^([a-zA-Z]+)\d+$', my_list[0])
    if not my_match:
        return " ".join(my_list)
    else:
        prefix = my_match[0]
        prelen = len(prefix)
        numbers = [ short_check[x] for x in my_list]
        lastno = None
        inrun = 0
        sofar = ""
        for i in range(len(my_list)):
            crnt_string = my_list[i]
            crnt_num = numbers[i]
            if (crnt_num != None and lastno != None and crnt_num == lastno +1 ):
                if (len(my_list) -1) == i:
                    sofar += ":" + str(crnt_num)
                inrun = 1
                lastno = crnt_num
            else:
                if lastno != None:
                    if inrun > 0:
                        sofar += ":" + str(lastno)
                    lastno = None
                inrun = 0
                if sofar != "":
                    sofar += " "
                sofar += crnt_string
                lastno = crnt_num
        return sofar

        

# If prefix is equal to inp_str pre_len and only digits follow then return digits
def short_check(inp_str, prefix, prelen):
    if inp_str[0:prelen] == prefix and (re.search(r'^\d+$', inp_str[prelen:])):
        return inp_str[prelen:]
    else:
        return None


#pool is a dict
# returns a list of sorted arrays from pool
def special_sorted_pool_keys(pool):
    all_keys = pool.keys()
    list_of_lists_and_keys = []
    for k in all_keys:
        if len(pool[k]) == 4:
            list_of_lists_and_keys.append(pool[k] + [k])

    mergeSort(list_of_lists_and_keys, 0, len(list_of_lists_and_keys)-1)
    sorted_keys = []
    for l in list_of_lists_and_keys:
        sorted_keys.append(l[4])

    return sorted_keys



#key input is a key from pool
#list_a and list_b should be lists
def compare_two_lists(list_a, list_b):
    if not (isinstance(list_a, list) and isinstance(list_b, list)):
        raise Exception("function compares lists")
    a = cmp(list_a[1], list_b[1]) #Sort by scaffold first 
    if a == 0: # if scaffold is the same, we sort by...
        if list_a[1] == "pastEnd":
            b = cmp(list_a[2], list_b[2]) #sort by strand
            if b == 0:
                d = cmp(list_a[0], list_b[0]) #sort by scaffold last
                return d
            else:
                return b
        else:
            c = cmp(int(list_a[3]), int(list_b[3])) # sort by position
            if c == 0:
                b = cmp(list_a[2], list_b[2]) #sort by strand
                if b == 0:
                    d = cmp(list_a[0], list_b[0]) #sort by scaffold last
                    return d
                else:
                    return b
            else:
                return c
    else:
        return a


def cmp(a, b):
    return (a > b) - (a < b)


def merge(arr, l, m, r): 
    n1 = int(m - l + 1)
    n2 = int(r- m)
  
    # create temp arrays 
    L = [0] * (n1) 
    R = [0] * (n2) 
  
    # Copy data to temp arrays L[] and R[] 
    for i in range(0 , n1): 
        L[i] = arr[l + i] 
  
    for j in range(0 , n2): 
        R[j] = arr[m + 1 + j] 
  
    # Merge the temp arrays back into arr[l..r] 
    i = 0     # Initial index of first subarray 
    j = 0     # Initial index of second subarray 
    k = l     # Initial index of merged subarray 

  
    while i < n1 and j < n2 :
        if compare_two_lists(L[i], R[j]) == -1: 
            arr[k] = L[i] 
            i += 1
        else: 
            arr[k] = R[j] 
            j += 1
        k += 1
  
    # Copy the remaining elements of L[], if there 
    # are any 
    while i < n1: 
        arr[k] = L[i] 
        i += 1
        k += 1
  
    # Copy the remaining elements of R[], if there 
    # are any 
    while j < n2: 
        arr[k] = R[j] 
        j += 1
        k += 1
  
# l is for left index and r is right index of the 
# sub-array of arr to be sorted 
#arr is array
def mergeSort(arr,l,r): 
    if l < r: 
  
        # Same as (l+r)/2, but avoids overflow for 
        # large l and h 
        m = (l+(r-1))//2
  
        # Sort first and second halves 
        mergeSort(arr, l, m) 
        mergeSort(arr, m+1, r) 
        merge(arr, l, m, r) 

"""
Vars in use:
    pool (dict)
    counts (dict)
    nSamples (int)
"""
def write_to_poolcount(pool, counts, nSamples, out_name, indexes, vars_dict):

    #First line of poolcount file: barcode, etc. tab separated.
    poolcount_str = "\t".join(["barcode","rcbarcode", "scaffold", "strand", 
        "pos", "\t".join(indexes)]) + "\n"

    sorted_keys = special_sorted_pool_keys(pool) #MUST TEST SORTING
    
    for rcbarcode in sorted_keys:

        barcode_list = pool[rcbarcode]

        if rcbarcode in counts:
            count_num_list = counts[rcbarcode]
        else:
            count_num_list = None

        out_list = [ barcode_list[0], rcbarcode, barcode_list[1], 
                barcode_list[2], barcode_list[3] ]

        if count_num_list != None:
            if len(count_num_list) < nSamples:
                warning_text = "{}:\n Adding empty elements".format(rcbarcode)
                vars_dict["report_dict"]["warnings"].append(warning_text)
                count_num_list += [0] * (nSamples-len(count_num_list))
            out_list.append("\t".join([str(x) for x in count_num_list]))
        else:
            out_list.append("\t".join(["0"]*nSamples))

        poolcount_str += "\t".join(out_list) + "\n"

    out_fp = out_name + ".poolcount"
    with open(out_fp, "w") as g:
        g.write(poolcount_str)
    report_str = " Wrote {}.\n".format(out_fp)
    vars_dict["report_str"] += report_str
    logging.info(report_str)

    return vars_dict
 

def write_to_colsum(indexes, colSums, colSumsUsed, out_name, vars_dict):
    colsum_str = "\t".join(["Index", "nReads", "nUsed", "fraction"]) + "\n"
    for i in range(len(indexes)):
        colsum_str += "\t".join([str(indexes[i]),str(colSums[i]),
            str(colSumsUsed[i]), str((colSumsUsed[i]/(0.5 + colSums[i])))]) + "\n"
    colsum_fp = out_name + ".colsum"
    with open(colsum_fp, "w") as g:
        g.write(colsum_str)
    report_str = " Wrote {}.colsum \n".format(out_name)
    vars_dict["report_str"] += report_str
    logging.info(report_str)

    return vars_dict


def get_usage_str(info_dict):

    usage_fp = os.path.join(info_dict["current_dir"],
            info_dict["texts"]["combine_barseq_usage_fp"])
    with open(usage_fp, "r") as f:
        usage_str = f.read()
    return usage_str

def unit_test():
    sample_inputs = ["myout my_pool codefile1,codefile2,codefile3",
            "",
            ""]
    main_run()


def main():
    unit_test()


if __name__ == "__main__":
    main()
