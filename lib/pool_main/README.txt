


    The program runs in the following order:
        MultiCodes, combineBarSeq
       
        MultiCodes
            Inputs: 
                FASTQ files which are separated by "index". The "index" is 
                    contained within the name of the FASTQ file. e.g.
                    "FEBA_BS_148_IT003_S99_L003_R1_001.fastq" would have
                    index "003". The function 'get_index_val' in the file
                    downloader.py in this directory.
                These FASTQ files contain reads which contain the barcodes
                    within them.

                You also need an indicator of what type of inner-indexing is 
                used (within the fastq reads to show barcodes):
                    -n25 or -bs3, which says something about how the barcodes
                    are wrapped within the FASTQ file.

                The names and locations for the outputs


            Outputs:
                3 Different files: 'codes', 'count', 'close'
                    Most important file is the 'codes' file, which has a list of
                        barcodes and their related counts
                    The 'count' file contains for each number of counts, how
                        many barcodes had that value. i.e. if 57 barcodes appeared
                        12 times throughout the FASTQ file, then next to the number 
                        12 you will have 57 in the counts file.
                    The 'close' file contains barcodes that are close to other barcodes
                        within the string. i.e. they have slight variation to
                        the original barcode within their 20 nt.

        combineBarSeq
            
            Inputs:
                The Pool File (A list of valid good barcodes to use)
                All the 'codes' files from the MultiCodes runs.
                A name for the output

                    







                



