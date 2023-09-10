from optparse import OptionParser
import os
import sys
import pandas as pd
import warnings

from pyscenic.cli import pyscenic

def get_rkdb_genes(rkdb_fname,rkdb_genes_fname):

    with open(rkdb_fname,'r') as rkdb_f:
        rkdb = pyscenic._load_dbs([rkdb_f])[0]
        rkdb_genes = rkdb.genes
        pd.DataFrame(rkdb_genes).to_csv(rkdb_genes_fname,sep=',',index=False,header=False)

def parseArgs(args):
    parser = OptionParser()

    parser.add_option('', '--rkdb_fname', type='str',default=None,
                      help='Path of input ranking database')
    parser.add_option('', '--rkdb_genes_fname', type='str',default=None,
                      help='Path to output files of genes in the database')
    
    (opts, args) = parser.parse_args(args)

    return opts, args
            

def main(args):
    opts, args = parseArgs(args)

    get_rkdb_genes(opts.rkdb_fname,opts.rkdb_genes_fname)
                        
if __name__ == "__main__":
    main(sys.argv)

