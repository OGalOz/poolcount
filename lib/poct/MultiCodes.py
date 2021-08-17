"""

Explanation of Code:



"""


# This file is a loose translation of Morgan Price's MultiCodes.pl into python.
import os
import logging
import re
import copy
import time

"""
Program Process:
    
    
"""


def RunMultiCodes(MC_d):
    """
    (fp stands for filepath)
    Inputs as listed in function 'CheckInputs'
    MC_d: (dict):

        out_prefix: (str) writes to this fp + '.codes', '.counts', '.close'
        maxReads: (int) max Reads from input FASTQ file (Also known as nLimit)
        [index_name]: (str) Name of index
        [indexfile_fp]: (str) Not in use* fp to an index file, which has a format as 
                            listed above XOR index_name. Rarely used*
        minQuality: (int) minQuality for FASTQ base pairs
        protocol_type: (str) 'custom' or 'dntag' or 'base' or 'bs3' or 'n25' or 'Unknown'
        [bs3_fp]: (str) Filepath to barseq3.index2 tsv file (Nec. if protocol_type is bs3)
        doOff1: (bool) - Checks whether or not we write to ".close"- performs analysis
            for likelyhood of OffBy1 errors.
        MC_seqs: (dict)
            dnt_pre: (str) dntag pre sequence
            dnt_post: (str) dntag post sequence
            bs_pre: (str) base pre sequence
            bs_post: (str) base post sequence
        debug: (bool)
        [pre_seq]: (str) Optional pre sequence (Nec. if protocol_type is custom)
        [post_seq]: (str) Optional post sequence (Nec. if protocol_type is custom)
        [nPreExpected]: int to which nPreExpectedMax/Min will be int +/- 2
            (Nec. if protocol_type is custom)
        fastq_fp: (str)

    Returns:
        Report_Dict: (d)
            fastq_fp: s 
            nReads: i
            nOff: d
                int -> int
            nUniqBC: i
            nMultiplexed: i
            [nWrongIndex2]: i
            [nWrongI2Perc]: f (percent)
            [nWrongPrePos]: i
            [nWrongPrePosPerc]: f (percent)
            [nOnce]: (i)
            [nCases]: i
            [nOff1Reads]: i
            [fOff1]: f
            EstimateDiversity_d: (d)
                [noise]: (d)
                    percent_noise: (f) (prcnt)
                    diversity: (f)
                    seen_once: (f)
                    seen_twice: (f)
                [good]: (d)
                    k_good_codes: (f)
                    percent_good_reads: (f) (prcnt)
            Estimate_Bias_d: (d)
                countSofar: (i)
                percent_codes: (f) prcnt
                k_codes: (f)
                percent_reads: (f) prcnt
                k_reads: (f)
    Description:
        *This is after removing multiplexing capabilities. 
        First we get the index preseq/postseq sequences
        to find the barcode, and also the index name.
        This part is run once per FASTQ file.
        We count the number of times each barcode
        found was seen. We find barcodes by looking
        for the preseq and postseq indicated in
        the function 'GetProtocolVariables'.
    """
    # PREPARATION PHASE ---------------------------
    vrs = CheckInputs(MC_d)

    # We get nPreExpectedMin/Max
    vrs = GetProtocolVariables(vrs)

    # We get prefix (essential list for downstream analysis)
    # vrs = IndexFileOrInameToPrefixInfo(vrs)

    # DEBUG
    #print(vrs)
    #sys.exit(0)


    # RUNNING ANALYSIS PHASE ----------------------------
    # We parse input fastq
    vrs = ParseFastqInput(vrs)

    #sys.exit(0)


    # WRITING TO OUTPUT -------------------
   
    # We write out to files
    # out.codes - Barcode and number of times it was seen
    vrs = WriteOutCodesGetnPerCount(vrs)

    # out.counts - For a given number of hits, how many barcodes match that? 
    vrs = WriteOutCounts(vrs)

    # This writes the closest variants to barcodes and the differences
    # in counts between them
    if vrs['doOff1']:
        vrs = WriteOutCloseGetOffBy1(vrs)
    else:
        vrs['offby1'] = {}

    # RETURN TO USER --------------------

    # We Start Report
    output_d = PrepareMCOutput(vrs)

    return output_d


def PrepareMCOutput(vrs):
    """
    In this function we run all the estimates and preparations
    Necessary keys in vrs: (d)
        minQuality
        codes
        nPerCount
        nOff
        Rd
        fastq_fp: s
    """

    vrs['nUniq'] = len(vrs["codes"].keys())

    vrs["Rd"]= InitReport(vrs)

    vrs["Rd"]["EstimateDiversity_d"] = EstimateDiversity(vrs)
  
    vrs["Rd"]["Estimate_Bias_d"] = EstimateBias(vrs)

    vrs["fastq_fp"] = vrs["fastq_fp"]

    return vrs["Rd"]
    

def EstimateDiversity(inp_d):
    """
    Args:
        inp_d: (d) Necessary keys:
            minQuality: (i)
            nOff: (d)
            nPerCount: (d)
            nUniq: (i)
            doOff1: (b)
            offby1: (d) 
                Barcode -> 1 if likely offby1 error
    Returns:
        ED_d: (d) Estimate Diversity dict. Optional keys:
            [noise]: (d)
                percent_noise: (f) (prcnt)
                diversity: (f)
                seen_once: (f)
                seen_twice: (f)
            [good]: (d)
                good_codes: (i)
                percent_good_reads: (f) (prcnt)
    """
    ED_d = {}
    if inp_d['minQuality'] > 0 and inp_d['nOff'][20] >= 1000 \
            and inp_d['nPerCount'][2] >= 10:
        ED_d["noise"] = {}
        for fNoise in [0, 0.005, 0.01, 0.02]:
            nNoise = int(0.5 + inp_d['nOff'][20] * fNoise)
            if nNoise > inp_d['nPerCount'][1]:
                continue
            noise_d = {
                "percent_noise": fNoise * 100,
                "diversity": ((inp_d['nUniq'] - nNoise + \
                                (inp_d['nPerCount'][1] - nNoise)**2)/(2 * inp_d['nPerCount'][2])),
                "seen_once" : (inp_d['nUniq'] - nNoise),
            }
            ED_d["noise"][fNoise] = noise_d

        ED_d["seen_twice"] = inp_d['nPerCount'][2]

    if (inp_d['minQuality']>0 and inp_d['nOff'][20]> 1000 \
            and inp_d['doOff1']):
        nGoodCodes = 0
        nGoodReads = 0
        for code in inp_d['codes'].keys():
            count = inp_d['codes'][code]
            if ((code not in inp_d['offby1']) and (count > 1)):
                nGoodCodes += 1
                nGoodReads += count

        ED_d["good"] = {
            "good_codes": nGoodCodes,
            "percent_good_reads":100.0 * float(nGoodReads)/float(inp_d['nOff'][20])
            }

    return ED_d


def EstimateBias(inp_d):
    """
    Args:
        inp_d: d Necessary keys:
            minQuality
            codes
            nPerCount
            nOff
    Returns:
        countSofar: (i)
        percent_codes: (f) prcnt
        k_codes: (f)
        percent_reads: (f) prcnt
        k_reads: (f)

    """
    EB_d = {"containsValues": False}

    nUniq = len(inp_d["codes"].keys())

    if inp_d['minQuality'] > 0 and nUniq >= 5000:

        # What fraction of reads are accounted for by the top 1% of strains?
        f = 0.01
        nReadsSofar = 0
        nCodesSofar = 0
        countSofar = -1
        nPerCount = inp_d['nPerCount']
        sorted_npc_keys = sorted(nPerCount.keys(), reverse=True)
        for count in sorted_npc_keys:
            nReadsSofar += count * nPerCount[count]
            nCodesSofar += nPerCount[count]
            countSofar = count
            if nCodesSofar >= f * nUniq:
                break
        
        #Estimate Bias dict. Note if nUniq >0 then nOff[20] > 0
        EB_d = {
            "containsValues": True,
            "countSofar": countSofar,
            'percent_codes': (100.0 * float(nCodesSofar))/float(nUniq),
            'k_codes': float(nCodesSofar)/1000.0,
            'percent_reads': (100.0 * float(nReadsSofar))/float(inp_d['nOff'][20]),
            'k_reads': float(nReadsSofar)/1000.0
                }

    return EB_d

def WriteOutCloseGetOffBy1(inp_d):
    """
    Args:
        inp_d: (d)
            out_prefix: (Str)
            doOff1: (b)
            nOff: (int)
            codes: (d)
                barcode -> num times seen 
    """
    offby1 = {} # barcode => 1 for likely off-by-1 errors

    out_close_fp = inp_d["out_prefix"] + ".close"
    OUTCLOSE = open(out_close_fp, "w")
    OUTCLOSE.write("\t".join(['code1', 'count1', 'code2', 'count2']) + "\n")

    nCases = 0
    nOff1Reads = 0
    codes = inp_d['codes']
    for code in codes.keys():
        count = codes[code]
        variants = GetVariants(code) # variants is a list
        for variant in variants:
            if variant in codes: 
                n1 = count
                n2 = codes[variant]
                out_close_str = "\t".join([str(code),
                    str(n1), str(variant), str(n2)]) + "\n"
                OUTCLOSE.write(out_close_str)
                if n1 < n2:
                    offby1[code] = 1
                    nOff_val = n1
                else:
                    offby1[variant] = 1
                    nOff_val = n2
                nCases += 1
                nOff1Reads += nOff_val

    OUTCLOSE.close()


    if 20 in inp_d['nOff']:
        # Almost always should be the case
        denominator = inp_d['nOff'][20]
        if denominator == 0:
            logging.warning("Warning, no barcodes length 20 found.")
            denominator = 1
    else:
        denominator = 1

    fOff1 = str(float(nOff1Reads)/float(denominator))

    logging.info("Wrote .close to " + out_close_fp)
    
    # int, int, float
    inp_d['nCases'] = nCases
    inp_d['nOff1Reads'] = nOff1Reads
    inp_d['fOff1'] = fOff1
    inp_d['offby1'] = offby1

    return inp_d
    
    

def WriteOutCounts(inp_d):
    """
    Args:
        inp_d: (d) Must contain keys
            out_prefix: (str)
            nPerCount: (d)
                sum of counts per index (i) - > number of codes with that sum (i)
            Rd: (d) Report Dict
            codes: (d) Barcode -> list of indexes and numbers of times barcode w/ index

    Description:
        Writes to out.counts file, which holds, for each count,
        (which is a number of times a barcode was seen), how
        many barcodes matched that number. Should be higher
        closer to the smaller numbers.
    """

    out_counts_fp = inp_d["out_prefix"] + ".counts"
    OUTCOUNTS = open(out_counts_fp, "w")


    nUniqBC = len(inp_d['codes'].keys())

    OUTCOUNTS.write("\t".join(["Count","nCodes", "Frac"]) + "\n")

    sorted_npc_keys = sorted(inp_d['nPerCount'].keys())
    for count in sorted_npc_keys:

        out_counts_str = "\t".join([str(count),
            str(inp_d['nPerCount'][count]),
            str((inp_d['nPerCount'][count]/nUniqBC))
            ]) + "\n"
        OUTCOUNTS.write(out_counts_str)

    OUTCOUNTS.close()

    logging.info("Wrote out counts to " + out_counts_fp)


    return inp_d


def WriteOutCodesGetnPerCount(inp_d):
    """
    Args:
        inp_d: (d) Necessary keys:
            index_name (str): Name of index.
            out_prefix: (str) Path to out file
            codes: (d) barcode -> number of times seen 
            Rd: (d) report_dict
    Returns:
        nPerCount: (d)
            sum of counts per index (i) - > number of codes with that sum (i)
    Description:
        We write 'codes' dict out to a file.        
        File name is defined by key 'out_prefix' + '.codes'
        File looks like "barcode\tindex_name\n" Then
        one row per barcode and then number of times seen.
    """

    #Will add this to inp_d later
    nPerCount = {} # sum of counts per index - > number of codes with that count
    header_list = ["barcode", inp_d['index_name']]

    out_codes_fp = inp_d['out_prefix'] + ".codes"
    OUTCODES = open(out_codes_fp, "w")

    #initializing codes tmp string
    OUTCODES.write("\t".join(header_list) + "\n")

    # nPrefix = len(inp_d['prefix']) # will equal 1 if iname given

    for barcode in inp_d['codes'].keys():
        count = inp_d['codes'][barcode]
        current_row_list = [barcode, count]
        OUTCODES.write("\t".join([str(x) for x in current_row_list]) + "\n")
        if count in nPerCount:
            nPerCount[count] += 1
        else:
            nPerCount[count] = 1

    OUTCODES.close()
    
    logging.info("Wrote Out-codes to " + out_codes_fp)
    
    inp_d['nPerCount'] = nPerCount
    return inp_d






def CheckInputs(MC_d):
    """
    MC_d: (dict):

        out_prefix: (str) writes to this fp + '.codes', '.counts', '.close'
        fastq_fp: (str) Path to input fastq file
        maxReads: (int) max Reads from input FASTQ file (Also known as nLimit)
        [index_name]: (str) Name of index
        [indexfile_fp]: (str) fp to an index file, which has a format as listed above XOR index_name
            Rarely used*
        minQuality: (int) minQuality for FASTQ base pairs
        protocol_type: (str) 'custom' or 'dntag' or 'base' or 'bs3' or 'n25' or 'Unknown'
        bs3_fp: (str) Filepath to barseq3.index2 tsv file 
        doOff1: (bool) - Checks whether or not we write to ".close"- performs analysis
            for likelyhood of OffBy1 errors.
        MC_seqs: (dict)
            dnt_pre: (str) dntag pre sequence
            dnt_post: (str) dntag post sequence
            bs_pre: (str) base pre sequence
            bs_post: (str) base post sequence
        debug: (bool)
        [pre_seq]: (str) Optional pre sequence (Nec. if protocol_type is custom)
        [post_seq]: (str) Optional post sequence (Nec. if protocol_type is custom)
        [nPreExpected]: int to which nPreExpectedMax/Min will be int +/- 2
            (Nec. if protocol_type is custom)
    """
    for x in ["out_prefix", "maxReads", "minQuality", "protocol_type", "doOff1",
                "MC_seqs", "debug"]:
        if x not in MC_d:
            raise Exception("Input to MultiCodes must contain " + x + ".")
    
    if not isinstance(MC_d["out_prefix"], str):
        raise Exception("File prefix to write out to must be string")
    if not os.path.exists(os.path.dirname(MC_d["out_prefix"])):
        raise Exception("Directory for MultiCodes output does not exist: {}".format(
            os.path.dirname(MC_d["out_prefix"])))
    if not isinstance(MC_d["fastq_fp"], str):
        raise Exception("FASTQ file path must be string")
    if not os.path.exists(MC_d["fastq_fp"]):
        raise Exception("FASTQ file path does not exist: {}".format(MC_d["fastq_fp"]))
    if not isinstance(MC_d["maxReads"], int):
        if MC_d["maxReads"] is not None:
            raise Exception("max Reads must be an int or null (None)")
    if not ( MC_d["maxReads"] is None or MC_d["maxReads"] > 0  ):
        raise Exception("max Reads must be greater than 0")
    if not "indexfile_fp" in MC_d:
        if not "index_name" in MC_d:
            raise Exception("indexfile_fp or index_name must be input to MultiCodes")
        elif not isinstance(MC_d["index_name"], str):
            raise Exception("index_name must be a string when input to MultiCodes")
    elif not "index_name" in MC_d:
        if not "indexfile_fp" in MC_d:
            raise Exception("indexfile_fp or index_name must be input to MultiCodes")
        elif not os.path.exists(MC_d["indexfile_fp"]):
            raise Exception("Index file does not exist")
    if not (isinstance(MC_d["minQuality"], int) and MC_d["minQuality"] >= 0):
        raise Exception("minQuality must be an integer that is equal to " \
                        + "or greater than 0. {}".format(MC_d["minQuality"]))
    if not isinstance(MC_d["protocol_type"], str):
        raise Exception("Protocol Type must be a string")
    protocol_option_list = ["custom", "dntag", "base", "bs3", "n25", "Unknown"]
    if not MC_d["protocol_type"] in protocol_option_list:
        raise Exception("Protocol Type must be one of:\n {}".format(
                    " ".join(protocol_option_list)))
    if not "bs3_fp" in MC_d:
        raise Exception("bs3_fp must be input to MultiCodes")
    elif not os.path.exists(MC_d["bs3_fp"]):
        raise Exception("BS3 index file does not exist at path " + MC_d["bs3_fp"])
    for x in ["doOff1", "debug"]:
        if not isinstance(MC_d[x], bool):
            raise Exception(x + " must be boolean.")
    for x in ["dnt_pre", "dnt_post", "bs_pre", "bs_post"]:
        if not x in MC_d["MC_seqs"]:
            raise Exception("MC_seqs must contain key " + x)
    if MC_d["protocol_type"] == "custom":
        for x in ["pre_seq", "post_seq"]:
            if x not in MC_d:
                raise Exception(x + " must be in MultiCodes input if protocol type is custom.")
            seq_b = True
            for y in MC_d[x].upper():
                if not y in ["A","C","T","G","N"]:
                    seq_b = False
            if not seq_b:
                raise Exception("custom " + x + " must contain only values 'ACTGN'")
        if not "nPreExpected" in MC_d or not isinstance(MC_d["nPreExpected"], int):
            raise Exception("nPreExpected must be int")
        if not MC_d["nPreExpected"] >= 0:
            raise Exception("nPreExpected must be greater than or equal to 0")
    elif MC_d["protocol_type"] == "dntag":
        if "index_name" in MC_d:
            raise Exception("If protocol is dntag then we "
                            "must have indexfile_fp, not index_name")

    return MC_d




def GetProtocolVariables(inp_d):
    """
    Args: (Essential)
        inp_d: (dict)
            protocol_type: (str)
                'custom' or 'dntag' or 'base' or 'bs3' or 'n25' or 'Unknown'
            [index_name]: (str) Index Name, e.g. ?? or None.
                iname is a necessary input IF protocol_type == "base" or "bs3"
            bs3_fp: (str) Filepath to barseq3 index2 tsv file
                this is a necessary input IF protocol_type == "bs3"
            MC_seqs:
                "dnt_pre": "GTCTCGTAG",
                "dnt_post": "CGATGAATT",
                "bs_pre": "CAGCGTACG",
                "bs_post": "AGAGACCTC"
            nPreExpectedMin (essential if protocol
            nPreExpectedMax


    You get the variables nPreExpectedMin(Max), and preseq and postseq 
        no matter what.
    If bs3 is the protocol type, then you also get index2 and index 2At
        
    """
    # str
    pc = inp_d['protocol_type'] 

    seq_dict = inp_d["MC_seqs"]
    if pc == "bs3":
        # This function fails if index name isn't in the bs3 file
        bs3_info_d = GetBarSeq3Info(bs3_info_d)
        nPreExpectedMin = bs3_info_d['nPreExpectedMin']
        nPreExpectedMax = bs3_info_d['nPreExpectedMax']
        preseq = seq_dict["bs_pre"]
        postseq = seq_dict["bs_post"]
    elif pc == "n25":
        nPreExpectedMin = 11
        nPreExpectedMax = 14
        preseq = seq_dict["bs_pre"]
        postseq = seq_dict["bs_post"]
    elif pc == "dntag":
        preseq = seq_dict["dnt_pre"]
        postseq = seq_dict["dnt_post"]
        nPreExpectedMin = 6
        nPreExpectedMax = 10
    elif pc == "base":
        if "index_name" in inp_d:
            #Usually here
            nPreExpectedMin = 12
            nPreExpectedMax = 16
        else:
            nPreExpectedMin = 7
            nPreExpectedMax = 11
        preseq = seq_dict["bs_pre"]
        postseq = seq_dict["bs_post"]
    elif pc == "custom":
        # if 'preseq' and 'postseq' are empty, then what?
        preseq = inp_d['preseq']
        postseq = inp_d['postseq']
        if inp_d['nPreExpected'] == None:
            raise Exception("Missing -nPreExpected but has preseq and postseq")
        nPreExpectedMin = inp_d['nPreExpected'] - 2
        nPreExpectedMax = inp_d['nPreExpected'] + 2
    else:
        raise Exception("Could not recognize protocol type: '" + pc + "'")

    #Updates
    # DNA sequences:
    inp_d['preseq'] = preseq
    inp_d['postseq'] = postseq 
    # ints:
    inp_d['nPreExpectedMin'] = nPreExpectedMin
    inp_d['nPreExpectedMax'] = nPreExpectedMax

    
    logging.basicConfig(level=logging.DEBUG)
    logging.info("preseq, postseq, nPreExpectedMin, nPreExpectedMax:")
    logging.info(f" {preseq}, {postseq}, {nPreExpectedMin}, {nPreExpectedMax}")

    return inp_d



def GetBarSeq3Info(inp_d):
    """
    Args:
        inp_d: (dict)
            bs3_fp: (str) filepath to barseq3.index2
            index_name: (str) Index name

    Returns:
        new_d (dict):
           nPreExpectedMin (int)
           nPreExpectedMax (int)
           [index2] (str): DNA sequence
           [index2At]:

    Description:
        If the computed index name is not found in the 
        BarSeq3 index file, then the program will halt.

        In the case that the second index is found,
            inp_d will have the keys index2 and index2At
            We will also update the variable index2len

    """
    new_d = {}

    #This is a special function to get the table fom file somehow
    required_bs3_file_columns = ["index_name", "index2", "nN"]
    #tab is a list of dicts with each line from ifile
    bs3_rows = ReadTable(inp_d['bs3_fp'], required_bs3_file_columns)
    #tabMatch is a subset of rows in bs3 with dicts whose index_name matches
    # main index_name
    tabMatch = []
    for row_d in bs3_rows:
        if row_d["index_name"] == inp_d["index_name"]:
            tabMatch.append(row_d)
            break

    if len(tabMatch) == 0: 
        # Could not find index name in the barseq3 file.
        warning_str = "Warning! Ignoring the second index -- index_name " \
             + "{} does not appear in {}\n".format(inp_d['index_name'], 
                    inp_d['bs3_fp'])
        #inp_d["report_dict"]["warnings"].append(warning_str)
        logging.warning(warning_str)
        error_string = f"Index name {inp_d['index_name']} not found in barseq3 " + \
                        f"index file at {inp_d['bs3_fp']}"
        # Getting out of program.
        raise Exception(error_string)
        #new_d['nPreExpectedMin'] = 16
        #new_d['nPreExpectedMax'] = 19
    else:
        #length of tabMatch is exactly 1
        row = tabMatch[0]
        inp_d['index2'] = row['index2'] #str
        inp_d['index2At'] = int(row['nN']) # (Number of ns before)
        nPE_preval = inp_d['index2At'] + 6 + 9
        new_d['nPreExpectedMin'] = nPE_preval - 2
        new_d['nPreExpectedMax'] = nPE_preval + 2
    return new_d


# NOT IN USE
def MyGrep(hash_list, index_name, iname):
    """
    hash_list: (list<subdict>)
        subdict: (dict)
            header -> value
    index_name: (str) key in subdict that maps to index names
    iname: (str) the index name we want.

    Returns:
        list<dict> minimized such that only dicts with wanted index name
            are in it.
    """
    new_list = []
    for h in hash_list:
        if h[index_name] == iname:
            new_list.append(h)
    return new_list



# fp is filepath, required is a list of required headers
# returns a list of hashmaps 
def ReadTable(fp, required):
    """
    fp: (str) Filepath to TSV file
    required: (list<str>) A list of required headers in string format.

    Returns:
        rows: list<dicts>
            each dict is a map of headers to values, which are all strings
    
    In the case of index TSV reading- looks for index_name, index2 and nN
    """

    FH = open(fp, "r")
    
    header_line = FH.readline().rstrip()
    header_list = header_line.split("\t")

    # columns map
    cols_map = {}
    for i in range(len(header_list)):
        cols_map[header_list[i]] = i

    for field in required:
        if field not in cols_map:
            raise Exception("No field {} in {}".format(
                field, fp))

    rows = []
    c_line = FH.readline()
    while c_line != "":

        line = c_line.rstrip()
        new_line_list = line.split("\t")
        if len(new_line_list) != len(header_list):
            raise Exception("Wrong number of columns in:\n{} \
                    \n in {}".format(line, fp))
        row_dict = {}
        for i in range(len(new_line_list)):
            row_dict[header_list[i]] = new_line_list[i]
        rows.append(copy.deepcopy(row_dict))

        c_line = FH.readline()

    return rows


# Deprecated
def IndexFileOrInameToPrefixInfo(inp_d):
    """
    There are two options: Either an index file path to get index information,
        or an index name, in which case prefix is a list with length 1 
    Args:
        inp_d: (dict)
            [indexfile_fp] (str)

            OR

            index_name: (str)

    """


    # Report String doesn't exist yet
    if "index_name" in inp_d:
        logging.info("Recognized index type as index name")
        # nLeading (int), Some sequence (str), index_name (str)
        prefix = [[0, "", inp_d['index_name']]]
        prefixNames = [inp_d['index_name']]
    else:
        raise Exception("Program only works with index_names now,"
                        " index file deprecated.")

    #updating inp_d
    inp_d['prefix'] = prefix
    inp_d['prefixNames'] = prefixNames

    return inp_d



def InitReport(inp_d):
    """
    Args:
        inp_d: (dict)
            fastq_fp: (str)
            nOff: (dict)
                ints (18-23) mapping to number of times (int) this distance
                    occured between preseq and postseq within the reads

            codes: (dict)
                barcode (str) -> prefixlist
                    prefixlist: list<int>
                        such that at the index of a DNAindex from the list 'prefix',
                            the number of times that barcode is found with that DNAindex
            nReads: (int)
            nMulti: (int)
            [index2]: (str)
                [nWrongIndex2]: (int)

    Returns:
        Rd: (d) Report Dict
            fastq_fp: s
            nReads: i
            nOff: d
                int -> int
            nUniqBC: i
            nMultiplexed: i
            [nWrongIndex2]: i
            [nWrongI2Perc]: f (percent)
            [nWrongPrePos]: i
            [nWrongPrePosPerc]: f (percent)
            [nOnce]: (i)
            [nCases]: i
            [nOff1Reads]: i
            [fOff1]: f


    """
    
    # Report Dict
    Rd = {}

    # Total differences between preseq and postseqs found
    nOffTot = sum(inp_d['nOff'].values())
    # Number of unique barcodes
    nUniqBC = len(inp_d['codes'].keys())
  
    stringed_nOff = {str(x):inp_d['nOff'][x] for x in inp_d['nOff'].keys()}

    Rd["fastq_fp"] = inp_d["fastq_fp"]
    Rd['nReads'] = inp_d['nReads']
    Rd['nOff'] = stringed_nOff 
    Rd['nOffTot'] = nOffTot 
    Rd['nUniqBC'] = nUniqBC
    Rd['nMultiplexed'] = inp_d['nMulti']
    offby1vals = False
    if "nCases" in inp_d:
        offby1vals = True
        for x in ["nCases", "nOff1Reads", "fOff1"]:
            if x in inp_d:
                Rd[x] = inp_d[x]
    Rd["offby1vals"] = offby1vals


    if inp_d['nReads'] > 0 and "index2" in inp_d:
        Rd['nWrongIndex2'] = inp_d['nWrongIndex2']
        Rd['nWrongI2Perc'] = (100*inp_d['nWrongIndex2']/inp_d['nReads'])

    if inp_d['nReads'] > inp_d['nWrongIndex2']:
        Rd['nWrongPrePos'] = inp_d['nWrongPrePos']
        Rd['nWrongPrePosPerc'] = 100*(inp_d['nWrongPrePos']/(inp_d['nReads'] - inp_d['nWrongIndex2'] ))

    if 1 in inp_d['nPerCount']:
        Rd["nOnce"] = inp_d['nPerCount'][1]

    return Rd

def ParseFastqInput(inp_d):
    """
    Args:
        inp_d: (dict) Necessary keys: (There are others unused)
            fastq_fp: (str) path to FASTQ file
            maxReads: (int) A limit on number of reads from FASTQ
            index_name: (str) Name of the index
            [index2]: (str) A sequence of length 6 (if protocol_type == 'bs3')
                [index2At]: (int) location of index2 (if protocol_type == 'bs3')
            debug: (bool)
            (For FindBarcode:)
            nPreExpectedMin: (int)
            nPreExpectedMax: (int)
            minQuality: (int)

    Returns:
        nWrongPrePos (int): 
        nWrongIndex2 (int): 
        nPostSeqUsed (int): How many times we used the postSeq while finding
                            the read.
        nReads (int): Total number of reads 
        nMulti (int): Number of reads with the prefix found (should be many). 
        nOff (dict): mapping from number 18-24 to number of times that offset seen.
        codes (dict):  barcode mapped to number of times it was found.

    Description:
        In this function we parse the FASTQ file.
        In the case of bs3 (BarSeq3) we have a 
        key called 'index2' which is mapped 
        to a DNA sequence of length 6, and 'index2At'
        which points to the expected location of index2
        within a read.
        If, within a read, we don't find index2 at
        its expected location, then that's noted as an error, 
        and we add a number to the variable 'nWrongIndex2'.


    """

    nWrongPrePos = 0
    nWrongIndex2 = 0
    nPostSeqUsed = 0
    nReads = 0
    nMulti = 0 # count with prefix identified
    nOff = {i:0 for i in range(18,24)}
    codes = {} #barcode maps to number of times seen

    # FastQ File Handle
    FQ_FH = open(inp_d['fastq_fp'], "r")

    # current read name
    c_readname = FQ_FH.readline()
    line_num = 1
    crnt_time = time.time()
    while c_readname != '':
        # We keep track of the number of reads
        nReads += 1

        # We extract the read info from the file handle
        read_name = c_readname.rstrip()
        seq = FQ_FH.readline().rstrip()
        break_line = FQ_FH.readline().rstrip()
        quality = FQ_FH.readline().rstrip()

        CheckRead(read_name, seq, break_line, quality, 
                  inp_d['fastq_fp'], line_num)
        line_num += 4

        if inp_d['maxReads'] is not None and (
                nReads >= inp_d['maxReads']):
            break

        # Index 2 is a DNA sequence of length 6 - Only exists with bs3
        # protocol under certain conditions
        if "index2" in inp_d:
            part = seq[inp_d['index2At']: inp_d['index2At'] + len(inp_d['index2'])]
            if part != inp_d['index2']:
                nWrongIndex2 += 1
                c_readname = FQ_FH.readline()
                continue

        nMulti += 1
        '''
        if "index_name" in inp_d:
            iPrefix = 0
        else:
            # Only if not multiplexed
            iPrefix = FindPrefix(seq, inp_d['prefix'], 
                        inp_d['debug']) # returns match (int) or -1
        '''


        #In most cases, index_name is used, so iPrefix is 0, 
        # so we get nLeading = 0, indexseq = "", prefixName is index_name

        nLeading, indexseq, prefixName = [0, "", inp_d["index_name"]]

        # Meaning that offset = nLeading + len("") = 0 + 0 = 0
        #offset = nLeading + len(indexseq)

        # barcode is str barcode, off is distance from start of barcode to start of postseq 
        # i.e. distance between end of preseq and beginning of postseq.
        barcode, off, postseqIgnored = FindBarcode(seq, quality, inp_d)
        if not postseqIgnored:
            nPostSeqUsed += 1

        if (barcode is None and off is not None and off >=0):
            nWrongPrePos += 1
            if nWrongPrePos >= 200 and nWrongPrePos >= (0.1 * nReads):
                raise Exception("Over 10% of reads have the wrong spacing ( \
                        not {}:{}) to the pre-sequence ({} \
                        of {} so far).\n Maybe the wrong \
                        protocol indicated (i.e., 'n25' or 'bs3')?\n".format(
                            inp_d["nPreExpectedMin"],
                            inp_d["nPreExpectedMax"],
                            nWrongPrePos,
                            nReads))
        if barcode is None:
            c_readname = FQ_FH.readline()
            continue

        # We count the number of times a distance between preseq and postseq occurs
        nOff[off] += 1

        if off != 20:
            # Barcode length is not 20
            c_readname = FQ_FH.readline()
            continue

        if barcode in codes:
            codes[barcode] += 1
        else:
            codes[barcode] = 1

        if nReads % 1000000 == 0:
            print("Read {} Reads so far".format(nReads))
            new_time = time.time()
            time_dif = int(new_time - crnt_time)
            print("Time: {} seconds".format(time_dif))
            crnt_time = new_time

        # Prepare for next read loop
        c_readname = FQ_FH.readline()



    FQ_FH.close()

    
    inp_d['nWrongPrePos'] = nWrongPrePos
    inp_d['nWrongIndex2'] = nWrongIndex2
    inp_d['nPostSeqUsed'] = nPostSeqUsed
    inp_d['nReads'] = nReads
    inp_d['nMulti'] = nMulti
    inp_d['nOff'] = nOff
    inp_d['codes'] = codes

    #DEBUG
    for x in ["nWrongPrePos", "nWrongIndex2", "nReads", "nMulti"]:
        print(x + ": " + str(inp_d[x]))

    
    return inp_d    



# Deprecated!
# seq is DNA sequence (str) 
# indexes is a list of lists, internal lists are [nLeading, indexseq, name],
#   where nLeading is an int, indexseq is str, name is str 
# debug in inp_d, True or False 
# Deprecated!
def FindPrefix(seq, indexes, debug):
    """
    Note this function is only run if the FASTQ file
    isn't multiplexed and all are in one.
    Args:
        seq: (str) DNA sequence (from FASTQ)
        indexes: list<list> internal lists are [nLeading, indexseq, name]
                                                ([int, str, str])
                the list is the same as the variable 'prefix'
        debug: (bool) 
    Returns:
       An int, the index within the list of indexes which
       matches this current index.
    """
    matches = []
    for i in range(len(indexes)):
        nLeading, indexseq, name = indexes[i]
        if (seq[nLeading:(nLeading + len(indexseq))] == indexseq):
            matches.append(i)
    if len(matches) == 1:
        return matches[0]
    elif len(matches) > 1:
        logging.critical("More than one match for DNA sequence:\n {}\n with indexes".format(
            seq))

    if debug:
        logging.info("No prefix for {}\n".format(seq))
    return -1




def UpdateCodesiPrefix(codes_dict, barcode):
    """
    Args:
        codes_dict: (dict)
            barcode (str) -> prefixlist
                prefixlist: list<int>
                    such that at the index of a DNAindex from the list 'prefix',
                        the number of times that barcode is found with that DNAindex
                        is being updated
                    
        barcode: (str)
    """


    return codes_dict 



def IndexFileToPrefixInfo(index_fp):
    """
    This function does not make sense given extant indexfiles
    Args:
        index_fp: (str) filepath to index file (TSV) with two columns
    Returns:
        ret_d: (dict)
            "report_str": (str)
            "prefixNames": list<str>,
            "prefix": list<[int, str, str]>
                [nLeading, indexseq, name] e.g. [

    """

    IX_FH = open(index_fp, "r")

    header_line = IX_FH.readline()

    c_line = "placeholder"

    # prefix is an important list that holds [[nLeading i, indexseq s, name s],...]
    # nLeading is number of n's before index 
    prefix = []
    line_num = 0

    while c_line != "":
        c_line = IX_FH.readline().rstrip()
        line_num += 1

        line_split = c_line.split('\t')

        if len(line_split) > 2:
            raise Exception("In indexfile, found a line that has more than "\
                    + "2 tsvs.\n Filename: {} Line Number: {}".format(
                        index_fp, line_num))
        #Note name & index are in form  H1, ATCACGAG
        name, index = line_split 

        # What does this account for?
        if (re.search(r'name', name ,re.IGNORECASE)):
            continue

        nLeading = None
        indexseq = None

        match = re.search(r'^([nN]*)([ACGT]+)$',index)
        if not match:
            raise Exception("Invalid index sequence {}".format(index))
        else:
            nLeading = len(match[0])
            indexseq = match[1]

        if (nLeading == None ) or (indexseq == None) or (name == ''):
            raise Exception(line)
        prefix.append([nLeading, indexseq, name])

    IX_FH.close()

    report_str = "Read {} indices from {}\n".format(len(prefix),index_fp)
    prefixNames = [x[2] for x in prefix]

    
    return {
            "report_str": report_str,
            "prefixNames": prefixNames,
            "prefix": prefix
            }


#seq str - 
#quality str
#offset int, usually 0
def FindBarcode(seq, quality, inp_d, offset=0):
    """
    Args:
        seq: (str) Sequence from FASTQ Read
        quality: (str) Quality from FastQ Read
        inp_d: (dict) Required keys:
            preseq: (str) DNA Sequence 
            postseq: (str) DNA Sequence
            nPreExpectedMin: (int)
            nPreExpectedMax: (int)
            minQuality: (int)
            debug: (bool)
        offset: (int) Represents location after nLeading and index sequence

    preseq often: 'CAGCGTACG'
    postseq often: 'AGAGACCTC'

    How function works:
        Note: We assume FASTQs are demultiplexed so that the important
            part of the sequence is the sequence given, meaning that
            offset = 0.
        We get useful part of sequence and quality by taking original sequence and 
        starting from location after leading n's and index sequence.
        Then we find the index of presequence defined 
        earlier in the program within the useful part of the sequence and assign it's value
        to preseq_pos. preseq_pos must be within the nPreExpectedMin and
        nPreExpectedMax index values of the true sequence, or it's not viable
        and we return [None, None].
        Then we assign barcode to where the presequence ends to 20 indeces
        after that. If the entire sequence
        is too short and we don't get a full barcode, then we return 
        [None, None]. 
        Then we check for the postsequence in the true sequence. If the length 
        of the useful sequence is less than the sum of the barcode (20) and the
        presequence and the index of the presequence and the postsequence, we 
        try a shorter postsequence of length 4, if that's not found we
        claim it's too short, make postsequence "" and continue with function.

        foundEnd: foundEnd is the distance between the end of the presequence 
        and the beginning of the postsequence if it's found. If the 
        postsequence is nothing  (will this ever be the case in KBase?) then 
        foundEnd is automatically 20 

    Returns: [barcode, foundEnd, postseqIgnored]
        barcode (str): String of length 20, the barcode
        foundEnd (int): Difference between preseq and postseq
        postseqIgnored (bool): Whether or not we actually used postseq.
                                True means not used.


    """

    seq2 = seq[offset:]
    quality2 = quality[offset:]
    preseq = inp_d['preseq']

    # The following gives -1 if not found, index of beginning if found
    preseq_pos = seq2.find(preseq) 

    #Bad cases 
    if preseq_pos == -1 or preseq_pos < inp_d['nPreExpectedMin'] or \
            preseq_pos > inp_d['nPreExpectedMax']:
        if preseq_pos != -1:
            if inp_d['debug']:
                warning_str = "seq2 {} has invalid index-of-preseq:  {}\n".format(
                    seq2, preseq_pos)
                #inp_d["report_dict"]["warnings"].append(warning_str)
                logging.warning(warning_str)
            return [None, preseq_pos, True] # report that spacing was wrong
        return [None, None, True]
   
    
    #Note: Below line does not throw error if end index > len(seq2)
    barcode = seq2[preseq_pos + len(preseq): preseq_pos + len(preseq) + 20]
    if not (len(barcode) == 20 and re.search(r'^[A-Z]+$',barcode)):
        # Barcode could be too short - meaning read is too short.
        #note barcodes with ambiguous sequences are ignored
        # we might want to allow one N in the future
        if inp_d['debug']:
            warning_str = "seq2 {} has invalid barcode sequence {}\n".format(
                seq2, barcode)
            #inp_d["report_dict"]["warnings"].append(warning_str)
            logging.warning(warning_str)
            return [None, None, True]
    
    crnt_postseq = inp_d['postseq']


    if len(seq2) < preseq_pos + len(preseq) + 20 + len(crnt_postseq):
        crnt_postseq = crnt_postseq[:4] #first 4 characters.
        if inp_d['debug']:
            warning_str = "Using postseq {}\n".format(crnt_postseq)
            #inp_d["report_dict"]["warnings"].append(warning_str)
            logging.warning(warning_str)
        #There might be redundant code here in MultiCodes.pl which checks
        # above length (16 lines above) Redundant code left out of this file.
        #We check again with a shorter crnt_postseq 
        if len(seq2) < preseq_pos + len(preseq) + 20 + len(crnt_postseq):
            if inp_d['debug']:
                warning_str = "Ignoring postseq, too short\n"
                #inp_d["report_dict"]["warnings"].append(warning_str)
                logging.warning(warning_str)
            crnt_postseq = ""

    foundEnd = -1
    if crnt_postseq == "":
        foundEnd = 20
        postseqIgnored = True
    else:
        postseqIgnored = False
        #check 20 first in case end of barcode matches post-seq (AGAG issue)
        for off in [20,19,21,18,22]:
            start_slice = preseq_pos + len(preseq) + off
            end_slice = start_slice + len(crnt_postseq)
            if (seq2[start_slice:end_slice] == crnt_postseq):
                foundEnd = off
                break
    if foundEnd == -1:
        # Could not find postseq
        if inp_d['debug']:
            warning_str = "seq2 \n {} \n barcode \n{}\n postseq not found\n".format(
                seq2, barcode)
            #inp_d["report_dict"]["warnings"].append(warning_str)
            logging.warning(warning_str)
        return [None,None, True]


    if inp_d['minQuality'] > 0:
        # the sanger code for !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ  
        barqual = quality2[preseq_pos + len(preseq): preseq_pos + len(preseq) \
                + 20]

        # quality score is encoded as 33+
        scores = [ord(c) - 33 for c in barqual]
        for score in scores:
            if (score < 0 or score > 100):
                raise Exception("Invalid score {} from barcode {} " \
                       + "quality {}".format(score, barcode, barqual))
            if score < inp_d['minQuality']:
                if inp_d['debug']:
                    warning_str = "Low quality {} for barcode {} in " \
                         + "{}\n".format(score, barcode, seq2)
                    #inp_d["report_dict"]["warnings"].append(warning_str)
                    logging.warning(warning_str)
                return [None,None, True]

    
    return [barcode, foundEnd, postseqIgnored]

def CheckRead(read_name, seq, break_line, quality, fastq_fp, line_num):
    """
    Args are all fastq strings without ending linebreak
    line_num refers to line number of read_name (int)
    """
    if not read_name[0] == "@":
        raise Exception("Read name does not start with @, line # {}\n File: {}".format(
            line_num, fastq_fp))
    for x in seq.upper():
        if x not in ["A","C","T","G","N"]:
            raise Exception("Sequence value {} not recognized. Line # {}\n File: {}".format(
                x, line_num + 1, fastq_fp))
    if not break_line[0] == "+":
        raise Exception("Break line not '+'. Instead '{}'. Line # {}\n File: {}".format(
            break_line[0],line_num + 2, fastq_fp))
    if not len(quality) == len(seq):
        raise Exception("Quality line wrong length. Lines # {}\n File: {}".format(
            line_num + 3, fastq_fp))




#baseseq (str)
# returns all variants on sequence
def GetVariants(baseseq):
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
