#!python3

# This file is a baseline python file

import os
import sys
import logging
import json


def CreateHTMLString(pre_HTML_d):

    # style link-'n'-tag
    stl_lnk_n_tag = getStyle()
    HTML_str = "<!DOCTYPE html>\n<html>\n<head>\n{}\n</head>\n<body>\n".format(stl_lnk_n_tag)

    HTML_str += CreateHTMLStringBody(pre_HTML_d)

    HTML_str += "</body></html>"

    return HTML_str

def CreateHTMLStringBody(pre_HTML_d):
    """
    Inputs:
        pre_HTML_d: (dict)
            MultiCodes_reports_list: (list) Each element is a MC_report_d:
                MC_report_d: (dict)
                    fastq_fp: The path to the fastq file
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
                    EstimateBias_d: (d)
                        countSofar: (i)
                        percent_codes: (f) prcnt
                        k_codes: (f)
                        percent_reads: (f) prcnt
                        k_reads: (f)
            CombineBarSeqReport_d: (d)
                Indexes: i Number of
                Success: i Number of
                LowCount: i
                LowFraction: i
                Total_Reads: i
                [MedianFraction]: f
                [MedianSuccess]: f or str
                [LowCountMedian]: i
                codes_report_dict: (d)
                    codes file path (str)-> report_d
                        report_d: (d)
                            nThisFile: i Number of lines in this file


    We divide the HTML page into three sections:
        Models Report,
        Map Tn Seqs Reports, 
        and the Design Random Pool Report.
    Each will have a brief explanation as to what it did 
    """

    HTML_str = "<h1>Text Report on 'RBTnSeq Counting'</h1>\n"
    # MultiCodes info
    HTML_str += Create_HTML_MC(pre_HTML_d["MultiCodes_reports_list"])
    # Combine BarSeq info
    HTML_str += Create_HTML_CBS(pre_HTML_d["CombineBarSeqReport_d"])

    return HTML_str

    
# We make a nice table for every MapTnSeq Report
def Create_HTML_MC(MultiCodes_reports_list):
    """
    MapTnSeq_reports_list: (list) Each element is a MTS_report_d:
        MC_report_d: (dict)
            fastq_fp: The path to the fastq file

    """

    MC_HTML_str = "<h2>Reading FASTQ files </h2>\n" \
            + "<p>--------------------------------</p>\n"
    for MC_report_d in MultiCodes_reports_list:
        fq_fn = os.path.basename(MC_report_d["fastq_fp"]).split(".")[0]
        MC_HTML_str += "<h4> MultiCodes Report extracting reads with barcodes from {}</h4>\n".format(
                fq_fn)
        MC_HTML_str += Create_MC_Table(MC_report_d)
    MC_HTML_str += "<p>---------------------------------------</p>"

    return MC_HTML_str


def Create_MC_Table(MC_rd):
    """
    Inputs:
        MC_rd:
            fastq_fp: The path to the fastq file
            nReads: i
            nOff: d
                str(int) -> int
            nUniqBC: i
            nMultiplexed: i
            [nWrongIndex2]: i
            [nWrongI2Perc]: f (percent)
            [nWrongPrePos]: i
            [nWrongPrePosPerc]: f (percent)
            [nOnce]: (i)
            offby1vals: (b)
            [nCases]: i Exists if offby1vals True
            [nOff1Reads]: i Exists if offby1vals True
            [fOff1]: f Exists if offby1vals True
            EstimateDiversity_d: (d)
                [noise]: (d) Contains keys [0, 0.005, 0.01, 0.02] which map to data about them
                                each key represents an fNoise value (fraction)
                    float (fnoise) -> fnoise_d
                        fnoise_d:
                            percent_noise: (f) (already prcnt)
                            diversity: (f)
                            seen_once: (f)
                [seen_twice] (i): Exists if key 'noise' is in EstimateDiverity_d
                [good]: (d)
                    good_codes: (f)
                    percent_good_reads: (f) (prcnt)
            Estimate_Bias_d: (d)
                containsValues: (b) - If True then all below values exist
                countSofar: (i)
                percent_codes: (f) prcnt
                k_codes: (f) nCodesSoFar/1000
                percent_reads: (f) prcnt
                k_reads: (f) nReadsSoFar/1000 
    Outputs:
        HTML_str: (str)
        This should be a table with important information displayed neatly
            for user to understand.
    """
    #reads processed
    r = MC_rd["nReads"]
    # What if r == 0?
    if r == 0:
        return ""

    # Preparing Table
    HTML_l = ['<div class="rpt_tbl_div">']
    css_tbl_cls = "dataTable__table table table-striped table-bordered dataTable no-footer"
    HTML_l += ['<table class="{}">'.format(css_tbl_cls)]

    table_list = [
                    ['# Reads Processed ', prep_int(r)],
                    ['# Multiplexed ', prep_int(MC_rd['nMultiplexed'])],
                    ['# Reads with Usable Barcodes (length 20)', prep_int(MC_rd['nOff']["20"])],
                    ['Fraction of Reads with Usable Barcodes', 
                            str(float(MC_rd['nOff']["20"])/float(r))[:5]],
                    ["# Unique Barcodes", prep_int(MC_rd['nUniqBC'])]
                ]
                
    if "nOnce" in MC_rd:
        table_list.append(["# Of Barcodes seen once", prep_int(MC_rd["nOnce"])])

    if MC_rd["offby1vals"]:
        table_list.append(["# Off-By-1 Barcode Cases", prep_int(MC_rd["nCases"])])
        table_list.append(["# Off-By-1 Reads", prep_int(MC_rd["nOff1Reads"])])
        table_list.append(["Off-By-1 Reads / # Reads with Usable Barcodes", str(MC_rd["fOff1"])[:5]])

    if "nWrongIndex2" in MC_rd:
        table_list.append(["% Reads failed to match index-2 ",
                            Prc(float(MC_rd["nWrongIndex2"])/float(r)) + "%"])

        if r > MC_rd["nWrongIndex2"]:
            if "nWrongPrePos" in MC_rd:
                table_list.append(["% Wrong presequence position",
                            Prc(float(MC_rd["nWrongPrePos"])/float(r - MC_rd[
                                                            "nWrongIndex2"])) + "%"])
    for offDistance in MC_rd["nOff"].keys():
        table_list.append(["# reads with distance {} between preseq and postseq".format(
                            offDistance),
                            prep_int(MC_rd["nOff"][offDistance])])


    # Estimate Diversity Values
    if "noise" in MC_rd["EstimateDiversity_d"]:
        nd = MC_rd["EstimateDiversity_d"]["noise"]
        for fNoise in nd.keys():
            table_list.append(["If {}% of reads are noise: ".format(nd[fNoise]["percent_noise"]),
                                "Diversity: {}".format(nd[fNoise]["diversity"])])
            table_list.append(["If {}% of reads are noise: ".format(nd[fNoise]["percent_noise"]),
                                "# Barcodes seen once: {}".format(prep_int(nd[fNoise]["seen_once"]))])
    if "good" in MC_rd["EstimateDiversity_d"]:
        gd = MC_rd["EstimateDiversity_d"]["good"]
        table_list.append(["Besides singletons and off-by-1s, # Barcodes:",
                            prep_int(gd["good_codes"])])

        table_list.append(["Besides singletons and off-by-1s, % Reads:",
            str(gd["percent_good_reads"])[:5]])

    # Estimate Bias Values
    if MC_rd["Estimate_Bias_d"]["containsValues"]:
        ebd = MC_rd["Estimate_Bias_d"]

        table_list.append(["# Barcodes with >= {} reads each:".format(prep_int(ebd["countSofar"])),
                            prep_int(int(ebd["k_codes"]*1000))])
        table_list.append(["% Barcodes with >= {} reads each:".format(prep_int(ebd["countSofar"])),
                            Prc(float(ebd["k_codes"]*1000)/float(MC_rd["nUniqBC"]))])
        table_list.append(["# Reads with barcodes accounting for above barcodes:",
                            prep_int(int(ebd["k_reads"]*1000))])
        table_list.append(["% Reads with barcodes accounting for above barcodes:",
                            Prc(float(ebd["k_reads"]*1000)/MC_rd["nOff"]["20"])])


    
    # We create a table with first column holding primarily text, second column holding value
    for info in table_list:

        html_str = '<tr role="row" class="MC_row">\n' \
                + '<td class="MC_col_1">' + info[0] + '</td>\n' \
                + '<td class="MC_col_2">' + str(info[1]) + '</td>\n' \
                + '</tr>'
        HTML_l.append(html_str)

    HTML_l.append("</table>")
    HTML_l.append("</div>")

    return "\n".join(HTML_l)



def Create_HTML_CBS(CBS_rd):
    """ Returns HTML string of table with info for CombineBarSeq

        Args:
            CBS_rd: (d) Combine BarSeq Report Dict
                Indexes: i Number of
                Success: i Number of
                LowCount: i
                LowFraction: i
                Total_Reads: i
                codes_report_dict: (d)
                    codes file path (str)-> report_d
                        report_d: (d)
                            nThisFile: i Number of lines in this file
                [MedianFraction]: f
                [MedianSuccess]: f or str
                [LowCountMedian]: i
    """
    HTML_l = ["<h2>Pool Count Report</h2>"]

    HTML_l += ['<div id="CBS_table" class="rpt_tbl_div">']
    
    css_tbl_cls = "dataTable__table table table-striped table-bordered dataTable no-footer"
    HTML_l += ['<table class="{}">'.format(css_tbl_cls)]


    table_list = [
            ['# Indexes ', prep_int(CBS_rd['Indexes'])],
            ['# Successes ', prep_int(CBS_rd['Success'])],
            ['# Indexes with less than 2 * 10^5 Reads ', prep_int(CBS_rd['LowCount'])],
            ['# Indexes with fraction of used BC less than 0.25 of total', prep_int(CBS_rd['LowFraction'])],
            ['# Total Reads', prep_int(CBS_rd["Total_Reads"])]
        ]

    if "MedianFraction" in CBS_rd:
        table_list.append(["Median fraction of used Reads over total reads per index",
            str(CBS_rd["MedianFraction"])[:5]])
        table_list.append(["Median fraction like above but for succesful indexes",
            str(CBS_rd["MedianSuccess"])[:5]])
    if "LowCountMedian" in CBS_rd:
        table_list.append(["Median read count for indexes with less than 2 * 10^5 Reads",
            prep_int(CBS_rd["LowCountMedian"])])

    for info in table_list:
        html_str = '<tr role="row" class="CBS_row">\n' \
                + '<td class="CBS_col_1">' + info[0] + '</td>\n' \
                + '<td class="CBS_col_2">' + str(info[1]) + '</td>\n' \
                + '</tr>'
        HTML_l.append(html_str)

    HTML_l.append('</table>')
    HTML_l.append("</div>")

    return "\n".join(HTML_l)

    

def prep_int(inp_i):
    # inp_i is an int or float with no decimal nonzero digits, value will be converted into
    # string and have commas: 1000000 -> '1,000,000'
    # OR inp_i is '?'
    
    # converting floats into ints and ints into strings and strings into lists
    if inp_i != '?':
        x = list(str(int(inp_i)))
    else:
        return '?'

    op_str = ''
    while len(x) > 3:
        c_char = x.pop(0)
        op_str += c_char
        if len(x) % 3 == 0:
            op_str += ","

    op_str += ''.join(x) 

    return op_str


def Prc(flt):
    #flt is a fraction to be turned into percent and cut short to 3 decimals
    # returns string
    if flt <= 1.01 and flt > -0.01:
        flt = str(flt*100)
    else:
        return str("Percent Error? " + str(flt) )

    # We split the float into before decimal and after
    l = flt.split(".")
    # We round the after-decimal part
    ad = l[1][:2]
    if int(ad[-1]) > 4:
        ad = str(int(ad[0]) + 1)
    else:
        ad = str(int(ad[0]))

    op = l[0] + "." + ad

    return op 


def getStyle():
    """
    returns string of style. First links, and then the style tag
    """
    style_link = '<link rel="stylesheet" href="style.css">\n' 
    style_tag = '<style>\n.rpt_tbl_div {margin: auto; width:60%;}\n'
    style_tag += 'h1 {text-align:center;}\nh2 {text-align:center;}\nh3 {text-align:center;}\n'
    style_tag += 'h4 {text-align:center;}\np {text-align:center;}\n'
    style_tag += '</style>'

    opt_HTML_style = """
    table, th, td {
              border: 1px solid black;
                border-collapse: collapse;
                }
    th, td {
              padding: 15px;
                text-align: left;
                }
    #t01 tr:nth-child(even) {
      background-color: #eee;
          }
          #t01 tr:nth-child(odd) {
           background-color: #fff;
               }
               #t01 th {
                 background-color: black;
                   color: white;
                   }\n
                """
    return style_link + style_tag




def main():
    
    args = sys.argv
    if args[1] == "how":
        print("python3 HTMLReport.py HTML_d.json op_fp.html")
        sys.exit(0)
    
    with open(args[1], "r") as f:
        html_d = json.loads(f.read())

    HTML_str = CreateHTMLString(html_d)

    with open(args[2], "w") as f:
        f.write(HTML_str)

    print("Wrote HTML file to " + args[2])



    return None

if __name__ == "__main__":
    main()
