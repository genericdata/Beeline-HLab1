from optparse import OptionParser
import os
import sys
import pandas as pd
import warnings

from pyscenic.utils import load_motifs

def scenicdf_to_rankededges(scenic_motifs_fname,ranked_edges_fname):
    scenic_df = load_motifs(scenic_motifs_fname)

    if scenic_df.shape[0] == 0:
        mat_out = pd.DataFrame([],columns=['TF','Target','Score'])
        
    else:
        scenic_df = scenic_df.copy()
        scenic_df.index = scenic_df.index.droplevel(1)
        scenic_df.columns = scenic_df.columns.droplevel(0)

        mat_dict = {tf:pd.DataFrame(data,columns=['Target','Score']).set_index('Target')
                    for tf,data in scenic_df['TargetGenes'].iteritems()}
        mat_long = pd.concat(mat_dict,names=['TF','Target'])

        mat_out = mat_long.groupby(['TF','Target']).max()

    outDF = mat_out.reset_index()[['TF','Target','Score']].copy()
    outDF = outDF.sort_values(by=['Score','Target','TF'], ascending=[False,True,True])
    outDF.to_csv(ranked_edges_fname,
                 sep='\t',index=False,header=['Gene1','Gene2','EdgeWeight'])


def parseArgs(args):
    parser = OptionParser()

    parser.add_option('', '--scenic_motifs_fname', type='str',default=None,
                      help='Path SCENIC motifs file')
    parser.add_option('', '--ranked_edges_fname', type='str',default=None,
                      help='Path to ranked edges file')
    
    (opts, args) = parser.parse_args(args)

    return opts, args
            

def main(args):
    opts, args = parseArgs(args)

    scenicdf_to_rankededges(opts.scenic_motifs_fname,opts.ranked_edges_fname)
                        
if __name__ == "__main__":
    main(sys.argv)

