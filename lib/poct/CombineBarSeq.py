
#python3
# This program is a loose translation of Morgan Price's combineBarSeq.pl to 
#   python. Important inputs to this program are the pool_file and the 
#   codesfiles. Outputs are poolcount and colsum [and ignore for ignored barcodes] 

import os
import sys
import math
import argparse
from statistics import median
import logging
import re
import json


def RunCombineBarSeq(CBS_d):
    """
    Inputs:
        CBS_d: (d)
            out_prefix_fp: (s) Output PoolCount/Colsum/Ignore File to write to
            pool_fp: (s) Input pool file to parse
            codes_fp_l: list<code_fp> List of all codes filepaths
                code_fp: (str) Path to codes file
    Outputs:
        ctg_d: (d)
            Indexes: i
            Success: i
            LowCount: i
            LowFraction: i
            Total_Reads_in_mil: f
            [MedianFraction]: f
            [MedianSuccess]: f or str
            [LowCountMedian]: i
            codes_report_dict: (d)
                codes file path (str)-> report_d
                    report_d: (d)
                        nThisFile: i Number of lines in this file
   
    Description:
        Note all .codes files are barcode -> number of times it was seen.
        It has as many rows as it has different barcodes (almost always a different number
        than the pool file, not computed in the same way).
        We essentially create a simple dict called "counts" which maps barcodes to
        indexes and for each index how many times that barcode was found.
        Within the function ProcessAllCodeFiles we are returned the counts
        dict, which is inside the dict "cds_d".



    """


    # PREPARATION PHASE ---------------------------
    vrs = CheckInputs(CBS_d)

    # pool_d is {rcbarcode (s) => [barcode (str), scaffold (str), strand (str), pos (i)]}
    pool_d, rcbc_list = GetPoolDictFromPoolFile(vrs["pool_fp"])
    vrs["pool_d"] = pool_d

    # PROCESSING PHASE ----------------------------
    cds_d, codes_reports_dict = ProcessAllCodeFiles(vrs)

    # PRINTING OUT PHASE
    WriteToPoolCount(pool_d, 
                    rcbc_list,
                    cds_d["counts"], 
                    cds_d["nSamples"], 
                    vrs["out_prefix_fp"],
                    cds_d["indexes"])
    
    WriteToColSum(vrs["out_prefix_fp"], 
                  cds_d["indexes"], 
                  cds_d["colSums"], 
                  cds_d["colSumsUsed"])


    # REPORT PHASE ------------------------
    ctg_d = CategorizeResults(cds_d)
    ctg_d["codes_reports_dict"] = codes_reports_dict


    return ctg_d


def WriteToColSum(out_prefix, indexes, colSums, colSumsUsed):
    """
    Args:
        out_prefix: (str) Filepath to output minus .
        indexes: list<str> List of index names
        colSums: list<int> samples to total number of counts
        colSumsUsed: list<int> samples to total number of counts for used barcodes

    Description:
        Note that above 'samples to' means that the indexes of the lists are
        matched to the location of the indexes in the list indexes.

    """

    cs_fp = out_prefix + ".colsum"
    CS_FH = open(cs_fp, "w") 
    CS_FH.write("index_name\tcolSum\tcolSumsUsed\tfrac\n")
    for i in range(len(indexes)):
        CS_FH.write("\t".join([
            indexes[i],
            str(colSums[i]),
            str(colSumsUsed[i]),
            str(float(colSumsUsed[i])/(0.5 + float(colSumsUsed[i])))]
            ) + "\n")

    CS_FH.close()
    logging.info("Wrote colsum file to " + cs_fp)

    


def WriteToPoolCount(pool_d, rcbc_list,
                    counts, nSamples, out_prefix, indexes):

    """
    Inputs:
        pool_d (dict)
            rcbarcode (s) => [barcode (str), scaffold (str), strand (str), pos (i)]
        rcbc_list (list<str>): We get a list of all the rcbarcodes as they are
                             in the pool file (dicts might lose ordering) and
                             this way we can print them out in the same order later.
        counts (dict)
            rcbarcode to vector of counts
        nSamples (int)
        out_prefix: (str) Filepath to out file prefix (poolcount, colsum, ignore)
        indexes: (list<str>) List of index names

    Description:
        Note that no Thresholds are employed within this function.
    """
    
    out_fp = out_prefix + ".poolcount"
    PC_FH = open(out_fp, "w")

    #First line of poolcount file: barcode, etc. tab separated.
    PC_FH.write("\t".join(["barcode","rcbarcode", "scaffold", "strand", 
        "pos", "\t".join(indexes)]) + "\n")

    #sorted_keys = special_sorted_pool_keys(pool_d) 
    
    for rcbarcode in rcbc_list:

        barcode_list = pool_d[rcbarcode]

        if rcbarcode in counts:
            # list of counts related to indexes
            count_num_list = counts[rcbarcode]
        else:
            count_num_list = None

        out_list = [ barcode_list[0], rcbarcode, barcode_list[1], 
                barcode_list[2], str(barcode_list[3]) ]

        if count_num_list is not None:
            if len(count_num_list) < nSamples:
                count_num_list += [0] * (nSamples-len(count_num_list))
            out_list += [str(x) for x in count_num_list]
        else:
            out_list += ["0"]*nSamples

        
        PC_FH.write("\t".join(out_list) + "\n")

    PC_FH.close()
    logging.info("Wrote poolcount to " + out_fp)

    return None




def CategorizeResults(cds_d, minCountsPerIndex = 200000,
                    minRatio=0.25):
    """
    Args:
        cds_d: (d) Must contain following keys:
            indexes: list<str> vector of names of samples
            colSums: list<int> samples to total number of counts
            colSumsUsed: list<int> samples to total number of counts for used barcodes
                                    (meaning barcodes were in the poolfile)
        minCountsPerIndex (int): How many counts should be found in an index
                                to make it good enough to be used?
        minRatio (float): minimum Ratio between colSumsUsed/colSums for an index
                          to pass
    Outputs:
        ctg_d: (d)
            Indexes: i Number of total 'indexes' e.g. "IT001"
            Success: i Number of indexes with more than minCountsPerIndex instances and 
                        f > minRatio
            LowCount: i Which 'indexes' (names like IT001) have a total count of 
                        found barcodes that is too low. Barcodes don't have
                        to be in poolfile, this is for total counts over all
                        barcodes.
            LowFraction: i Which 'indexes' ("IT001") have a fraction of counts
                        for barcodes that were found in the poolfile as opposed
                        to total counts.
            Total_Reads: i Sum of counts for all barcodes for all indexes.
            [MedianFraction]: f Median of fracs for all indexes
            [MedianSuccess]: f Median of fracs for succesful indexes
            [LowCountMedian]: i Median of colSums for low count indexes.

    Description:
        At first we get important numbers like the Total number of indexes,
        total number of indexes that pass our thresholds,
        the indexes that fail to pass thresholds and why (low Count, low Ratio),
        and the total number of counts over all indexes.
        If the ratio of [counts for barcodes in pool file/ total counts for barcodes]
        is less than minRatio, then the index is considered to have failed the threshold
        test.
    """
    ix = cds_d["indexes"]
    # Total counts for barcodes in an index
    cS = cds_d["colSums"]
    # Counts for barcodes in an index that were also in poolfile
    cSU = cds_d["colSumsUsed"]

    lowcountI = [] # indexes with < minCountsPerIndex reads
    lowcounts = [] # those actual counts
    fractions = [] # fraction used for indexes with plenty of reads only 
    Succesfulfrac = [] # fraction used for successful indexes (enough reads and f > 0.25)
    lowFracIx = [] # indexes with fraction < minRatio
    okI = [] # indexes with enough reads and fraction above minRatio
    totalReads = 0

    # Note that what is considered a good fraction
    # is when the colSumsUsed/colSums >= minRatio
    for i in range(len(ix)):
        totalReads += cS[i]
        if cS[i] < minCountsPerIndex :
            lowcountI.append(ix[i])
            lowcounts.append(cS[i])
        else:
            frac = float(cSU[i])/float(cS[i])
            fractions.append(frac)
            if frac < minRatio:
                lowFracIx.append(ix[i])
            else:
                okI.append(ix[i])
                Succesfulfrac.append(frac)

    ctg_d = {
            "Indexes": len(ix),
            "Success": len(okI),
            "LowCount": len(lowcountI),
            "LowFraction": len(lowFracIx),
            "Total_Reads": totalReads,
    }
    if len(fractions) > 0:
        # Almost always true
        ctg_d["MedianFraction"] = median(fractions)
        if len(Succesfulfrac) > 0:
            mS = median(Succesfulfrac)
        else:
            mS = "NaN" # Not a number
        ctg_d["MedianSuccess"] = mS
    if len(lowcounts) > 0:
        ctg_d["LowCountMedian"] = median(lowcounts)

    return ctg_d






def ProcessAllCodeFiles(vrs):
    """
    Args:
        vrs: (d) must contain keys:
            codes_fp_l: (list<code_fp>) 
                code_fp: (str) Path to .codes file
                .codes files contain 2 columns: barcode and then number of times barcode seen.
                The column header of the second column is the name of the index.
            pool_d: (d)
                rcbarcode (s) => [barcode (str), scaffold (str), strand (str), pos (i)]
            out_prefix_fp: (str) Filepath prefix PoolCount, Colsum, Ignore

    Returns:
        cds_d:
            counts: d
            nSamples: i
            orig_cd_fp: (str) First Input Codes File Path
            indexes: list<str> List of index names
            colSums: list<int> Same length as indexes, Sum of columns of all lines per index
            colSumsUsed: list<int> Same length as indexes, sum of columns per used line
            nUsed: i number of Barcodes from .codes found in pool file
            nIgnore: number of Barcodes from .codes ignored because not found in pool file
        codes_report_d: (d)
            codes file path (str)-> report_d
                report_d: (d)
                    nThisFile: i Number of lines in this file

    
    Description:
        In this function, we initialize the variables to be used later in
        categorizing the results, as well as creating the "counts" dict,
        which is used in creating the poolcount file, which turns into
        the poolcount object.
        For each codes file generated by MultiCodes, we parse it, and 
        add the numbers from the barcodes in the codes file to the 
        counts dict. The codes files contain the barcode, and the 
        number of times it was seen, and each codes file is associated
        with an index (str like "IT002") which we add to the list 
        'indexes' found in cds_d.
        The variables colSums and colSumsUsed refer to the sums over all 
        the barcodes of a .codes file, but where:
        colSums adds the count of each barcode to its sum
        regardless of whether the barcode is in the poolfile,   
        colSumsUsed only adds the count of the barcodes if the barcodes
        are in the poolfile.

    """
    

    # counts: rcbarcode to vector of counts
    # indexes: vector of names of samples
    # colSums: samples to total number of counts
    # colSumsUsed: samples to total number of counts for used barcodes
    # nSamples: number of Samples so far
    # nUsed: number of Barcodes from .codes found in pool file
    # nIgnore: number of Barcodes from .codes not in pool file
    cds_d = {
        "counts" : {} ,
        "nSamples" : 0,
        "orig_cd_fp": vrs["codes_fp_l"][0],

        "indexes" : [],
        "colSums" : [],
        "colSumsUsed" : [],
        "nUsed" : 0,
        "nIgnore" : 0
    }


    # We write ignored barcodes to this file
    ignore_fp = vrs["out_prefix_fp"] + ".ignore"
    IgnoreFH = open(ignore_fp, "w")
    cds_d["IgnoreFH"] = IgnoreFH

    codes_reports_d = {}
    for cd_fp in vrs["codes_fp_l"]:
        cds_d, cf_report_d= ProcessCodeFile(cd_fp, cds_d, vrs["pool_d"])
        codes_reports_d[cd_fp] = cf_report_d

    IgnoreFH.close()
    logging.info("Wrote ignore to " + ignore_fp)

    return [cds_d, codes_reports_d]


def ProcessCodeFile(cd_fp, cds_d, pool_d):
    """
    Args:
        cd_fp: (str) codes filepath
            .codes files contain 2 columns: barcode and then number of times barcode seen.
            The column header of the second column is the name of the index.
        cds_d: (dict) Codes dict
            counts: (d) rcbarcode to vector of counts
            indexes: list<str> vector of names of samples
            colSums: list<int> samples to total number of counts (of what?)
            colSumsUsed: list<int> samples to total number of counts for used barcodes
            nSamples: (i) number of Samples so far
            nUsed: (i) number of Barcodes from .codes found in pool file
            nIgnore: (i) number of Barcodes from .codes not in pool file
            orig_cd_fp: (str) First codes fp in list
            [IgnoreFH]: File Handle to Ignore file (write)
            [oneperfile]: (b) After ProcessCodeFileHeaderLine this is in cds_d

        pool_d: (d)
            rcbarcode (s) => [barcode (str), scaffold (str), strand (str), pos (i)]

    Returns:
        cds_d: (Above)
        cf_report_d: (python dict)
            nThisFile: i (number of lines this file)

    Description:
        We get preliminary info from 'ProcessCodeFileHeaderLine' function
        like which index we're using and if it's the first file or not.
        It is most important we update the dict 'counts', which is
        a mapping of barcode (actually a rcbarcode) -> a list of counts
        in each index. The 'indexes' list is the one this mapping is 
        associated with.
        Note that colSums adds the count of each barcode to its sum
        regardless of whether the barcode is in the poolfile, whereas
        colSumsUsed only adds the count of the barcodes if the barcodes
        are in the poolfile.
    """
    
    CD_FH = open(cd_fp, "r")

    header_line = CD_FH.readline().rstrip()

    # cds_d is dict as inputted but updated by ProcessCodeFileHeaderLine
    # thisIndex is an int/None which represents index (str) from 'indexes' in cds_d,
    # e.g. if indexes is ["A", "B"], and the code file is index "B" then thisIndex 
    # is 1 (0-indexing)
    # thisIndex is None if it's not one index per file
    # Also, cds_d is updated to have the correct number of items in the lists.
    cds_d, thisIndex = ProcessCodeFileHeaderLine(header_line, cds_d, cd_fp)

    # Number of lines for this file (Starts at one because of header line)
    nThisFile = 1

    c_line = CD_FH.readline().rstrip()
    while c_line != "":
        nThisFile += 1
        c_list = c_line.split('\t')
        barcode = c_list[0] # actually rcbarcode
        index_counts = [int(x) for x in c_list[1:]]
        if not re.match(r'^[ACGTN]+$', barcode):
            raise Exception("Invalid Barcode: " + barcode)
        expected_idx_len = 1 
        if not len(index_counts) == expected_idx_len:
            raise Exception("Wrong number of columns in file " + cd_fp)

        if barcode in pool_d:
            cds_d["nUsed"] += 1
            cds_d["colSumsUsed"][thisIndex] += index_counts[0]
            UpdateCountsAndBarcode(cds_d, barcode, thisIndex, index_counts[0])
            """
            else:
                for i in range(len(cds_d["nSamples"])):
                    cds_d["colSumsUsed"][i] += index_counts[i]
                if barcode in cds_d["counts"]:
                    c_row = cds_d["counts"][barcode]
                    for i in range(len(cds_d["nSamples"])):
                        c_row[i] += index_counts[i]
                    cds_d["counts"][barcode] = c_row
                else:
                    cds_d["counts"][barcode] = index_counts
            """

        else:
            if "IgnoreFH" in cds_d:
                cds_d["IgnoreFH"].write("\t".join([barcode, "\t".join([str(x) for x in index_counts])]))
            cds_d["nIgnore"] += 1

        # colSums is updated regardless of whether the barcode is used or not.
        cds_d["colSums"][thisIndex] += index_counts[0]
        '''
        else:
            for i in range(len(cds_d["nSamples"])):
                cds_d["colSums"][i] += index_counts[i]
        '''
    
        c_line = CD_FH.readline().rstrip()


    if nThisFile == 1: 
        raise Exception("No entries in .codes file" + cd_fp)
    CD_FH.close()

    cf_report_d = {
            "nThisFile": nThisFile
    }

    return [cds_d, cf_report_d]


def UpdateCountsAndBarcode(cds_d, barcode, thisIndex, num_reads):
    """ Function updates the Counts Dict with the barcode

    Args: 
        cds_d:
            counts: (d) rcbarcode to vector of counts
            indexes: list<str> vector of names of samples
        barcode: (str) A/C/T/G/N length ~20
        thisIndex: i Location-index within 'indexes' list from cds_d which this
            codes file has info for
        num_reads: (i) Number of reads for this barcode and index
    
    Description:
        Could be redone with counts dict having barcodes pointing to another
        dict with indexes pointing to number of times seen instead of 
        a list. For now we use a list.
    """
    counts = cds_d["counts"]
    if barcode in counts:
        if thisIndex < len(counts[barcode]):
            counts[barcode][thisIndex] += num_reads
        else:
            counts[barcode] += (thisIndex - len(counts[barcode]) + 1) * [0]
            counts[barcode][thisIndex] += num_reads
    else:
        counts[barcode] = (thisIndex + 1) * [0]
        counts[barcode][thisIndex] += num_reads

    # Not essential:
    # cds_d["counts"] = cts

    return None


def ProcessCodeFileHeaderLine(hl, cds_d, cd_fp):
    """Processes header line (str) from a .codes file output from MultiCodes

    Description:
        We make sure the variables for cds_d are ready to process
        another file. The colSums and colSumsUsed lists have another item
        added to them.

    Args:
        hl: (str) TSV Header line without new-line 
        cd_fp: (str) Debugging
        cds_d: (dict) Codes dict
            counts: (d) rcbarcode to vector of counts
            indexes: list<str> vector of names of samples
            colSums: list<int> samples to total number of counts
            colSumsUsed: list<int> samples to total number of counts for used barcodes
            nSamples: (i) number of Samples so far (equivalent to number of codes files so far)
            nUsed: (i) number of Barcodes from .codes found in pool file
            nIgnore: (i) number of Barcodes from .codes not in pool file
            orig_cd_fp: (str) First codes fp in list
            [oneperfile]: (b) Only exists after first file parsed. If one index per file
    Returns:
        cds_d: (As above)
        thisIndex: (i)/None if oneperfile mode, which index name in (indexes) is this file?
    """
    header_cols = hl.split('\t')
    first = header_cols[0]
    index_name = header_cols[1]
    if cds_d["nSamples"] == 0:
        # We are processing the first file
        cds_d["indexes"] = [index_name]
        cds_d["nSamples"] = 1
        cds_d["colSums"] = [0]
        cds_d["colSumsUsed"] = [0]
        thisIndex = 0
    else:
        # Indexes is a list of the 'indexes' from the first input codes file.
        # We make a dict of index value (str) -> location (0-based) within indexes list.
        oldcol = { cds_d["indexes"][i]:i for i in range(len(cds_d["indexes"]))}
        if index_name in oldcol:
            thisIndex = oldcol[index_name]
        else:
            cds_d["nSamples"] += 1
            cds_d["indexes"].append(index_name)
            thisIndex = len(cds_d["indexes"]) - 1
            cds_d["colSums"].append(0)
            cds_d["colSumsUsed"].append(0)

    ''' 
    Deprecated multiplexing
    else:
        if not (len(cols) == len(cds_d["indexes"])):
            raise Exception("Expecting {} columns, but got {} columns in codes file {}".format( 
                len(cds_d["indexes"]), len(cols), cd_fp))
        for i in range(len(cols)):
            if not cols[i] == cds_d["indexes"][i]:
                raise Exception("Index mismatch in {} vs {} -- \n {} vs {}".format(
                    cd_fp, cds_d["orig_cd_fp"], cols[i], cds_d["indexes"][i] ))
    '''

    return [cds_d, thisIndex]


def GetPoolDictFromPoolFile(pool_fp):
    """
    Input:
        pool_fp: (str) Path to pool file
    Output:
        pool_d: (d)
            rcbarcode (s) => [barcode (str), scaffold (str), strand (str), pos (i)]
        rcbc_list (list<str>): We get a list of all the rcbarcodes as they are
                             in the pool file (dicts might lose ordering) and
                             this way we can print them out in the same order later.
    """

    pool_d = {} # pool dict is rcbarcode to [barcode, scaffold, strand, pos]

    # Pool File Handle
    P_FH = open(pool_fp, "r")



    # We check the header line
    h_line = P_FH.readline().rstrip()
    h_list = h_line.split("\t")
    if not h_list[0] == "barcode":
        raise Exception("Header line incorrect in pool file " + pool_fp)

    # We create this list to get the ordered rcbarcodes as they are in pool file.
    rcbc_list = []

    c_line = P_FH.readline().rstrip()
    line_num = 2

    while c_line != "":
        
        pool_d, rcbarcode = CheckPoolLineAddToDict(c_line, pool_d, line_num, pool_fp)
        rcbc_list.append(rcbarcode)

        c_line = P_FH.readline().rstrip()
        line_num += 1

    P_FH.close()

    if len(pool_d.keys()) == 0:
        raise Exception("No entries in pool file")

    return pool_d, rcbc_list



def CheckPoolLineAddToDict(pool_line, pool_d, line_num, pool_fp):
    """
    pool_line: (str) Entire line from pool file 
    pool_d: (dict)
        rcbarcode (str) => [barcode (str), scaffold (str), strand (str), pos (i)]
    line_num: (int) For debugging
    pool_fp: (str) For debugging

    Note: We start after the header line
    Description:
        We check
    """
    #We get first 7 columns of pool_line (out of 12)
    split_pool_line = pool_line.split("\t")[:7]

    #We remove spaces (?) Is this necessary:
    for x in split_pool_line:
        x.replace(' ', '')

    if len(split_pool_line) < 7:
        raise Exception("Pool Line has less than 7 tab separated values at " \
                + "line # {}.\n File {}".format(line_num, pool_fp))

    # We only care about the best location, not second best
    barcode, rcbarcode, ph_1, ph_2, scaffold, strand, pos = split_pool_line

    if not re.search(r"^[ACGT]+$", barcode):
        raise Exception("Invalid barcode: |{}| Line # {} File {}".format(barcode, 
                line_num, pool_fp))
    if not re.search( r"^[ACGT]+$",rcbarcode ):
        raise Exception("Invalid rcbarcode: |{}| Line # {} File {}".format(rcbarcode,
            line_num, pool_fp))
    if not (pos == "" or re.search( r"^\d+$", pos)):
        raise Exception("Invalid position: |{}| Line # {} File {}".format(pos,
            line_num, pool_fp))
    if pos == "":
        pos = "-1"
    if not (strand == "+" or strand == "-" or strand == ""):
        raise Exception("Invalid strand: |{}| Line # {} File {}".format(strand,
            line_num, pool_fp))
    if rcbarcode in pool_d:
        raise Exception("Duplicate rcbarcode: {} Line # {} File {}".format(rcbarcode,
            line_num, pool_fp))

    pool_d[rcbarcode] = [barcode, scaffold, strand, int(pos)]

    return pool_d, rcbarcode




def CheckInputs(CBS_d):
    """
    Inputs:
        CBS_d: (d) (Combine BarSeq dict)
            out_prefix_fp: (s) Output PoolCount/Colsum/Ignore File to write to
            pool_fp: (s) Input pool file to parse
            codes_fp_l: list<code_fp> List of all codes filepaths
                code_fp: (str) Path to codes file

    Description:
        First we check if directory in which we write the poolcount file and others 
        exists. Then we make sure the directory doesn't already have the files in 
        it, in case something really odd happened.
        Then we make sure the pool file exists.
        Finally, we make sure each .codes file exists.
    """
    
    if not os.path.exists(os.path.dirname(CBS_d["out_prefix_fp"])):
        raise Exception("Directory to write files out to does not exist: \n" \
                + os.path.dirname(CBS_d["out_prefix_fp"]))
    for x in ["poolcount", "colsum", "ignored"]:
        op_fp = CBS_d["out_prefix_fp"] + "." + x
        if os.path.exists(op_fp):
            raise Exception("Output path {} already exists.".format(op_fp))
    if not os.path.exists(CBS_d["pool_fp"]):
        raise Exception("Input pool file {} does not exist.".format(CBS_d["pool_fp"]))
    for x in CBS_d["codes_fp_l"]:
        if not os.path.exists(x):
            raise Exception("Input .codes file {} does not exist".format(x))

    return CBS_d




# returns a list of sorted arrays from pool
def special_sorted_pool_keys(pool_d):
    """
    pool_d: (d)
            rcbarcode (s) => [barcode (str), scaffold (str), strand (str), pos (i)]

    Outputs:
        list of sorted
    """
    all_keys = pool_d.keys()
    pool_dbl_list = [pool_d[k] + [k] for k in all_keys]
    mergeSort(pool_dbl_list, 0, len(pool_dbl_list)-1)
    sorted_keys = [l[-1] for l in pool_dbl_list]

    return sorted_keys


def mergeSort(arr,l,r): 
    if l < r: 
  
        # Same as (l+r)/2, but avoids overflow for 
        # large l and h 
        m = (l+(r-1))//2
  
        # Sort first and second halves 
        mergeSort(arr, l, m) 
        mergeSort(arr, m+1, r) 
        merge(arr, l, m, r) 


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
