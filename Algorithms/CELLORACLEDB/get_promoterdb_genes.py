from optparse import OptionParser
import os
import sys
import pandas as pd
import warnings

def get_promoter_base_grn_genes(db_fname,genes_fname):

    base_grn = pd.read_parquet(db_fname)

    base_grn['gene_short_name'].to_csv(genes_fname,index=False,header=False)

def parseArgs(args):
    parser = OptionParser()

    parser.add_option('', '--promoter_base_grn_fname', type='str',default=None,
                      help='Path of input promoter base GRN file')
    parser.add_option('', '--promoter_base_grn_genes_fname', type='str',default=None,
                      help='Path to output files of genes in the base GRN file')
    
    (opts, args) = parser.parse_args(args)

    return opts, args
            

def main(args):
    opts, args = parseArgs(args)

    get_promoter_base_grn_genes(opts.promoter_base_grn_fname,
                                opts.promoter_base_grn_genes_fname)
                        
if __name__ == "__main__":
    main(sys.argv)

