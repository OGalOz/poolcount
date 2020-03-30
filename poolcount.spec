/*
A KBase module: poolcount
*/

module poolcount {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_poolcount(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
