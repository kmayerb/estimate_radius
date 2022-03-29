"""
Previously analysis presented with the tcrdist3 manuscript considered computing radii 
from a set of synthetic TCRs that matched the TRBV and TRBJ gene usage frequency 
of a set of epitope-associated receptors. When only single chain data is available 
it may be motr useful to deeply estimate an appropriate distance radii for every 
clone using a V-J matched TCR set. We still use the approach of reweighting the 
radius estimate with priors on the frequency of each V-J pairing. 
"""
import platform as pt
import pandas as pd
import parmap
import os
from tcrsampler.sampler import TCRsampler
from tcrdist.background import get_stratified_gene_usage_frequency
from routine import vj_tr
import argparse


def main(p,
         pout, 
         f,
         n_cpus, 
         n_runs, 
         organism, 
         chains):

    fp = os.path.join(p,f)
    d  = pd.read_csv(fp)

    if 'count' not in d.columns:
        d['count'] = 1

    du = d.groupby(['cdr3_b_aa','v_b_gene','j_b_gene']).size().sort_values(ascending = False).reset_index(drop = False)
    du.columns = ['cdr3_b_aa','v_b_gene','j_b_gene', 'count']
    print("VJ GENE USAGE FREQUENCY PER UNIQUE CLONE")
    print(du.groupby(['v_b_gene','j_b_gene']).size().sort_values(ascending = False))

    # Load a sampler for the purposes of estimating V and J gene usage 
    # frequency in real repertoires, which will be used later for weightings

    ts = TCRsampler(default_background = 'britanova_human_beta_t_cb.tsv.sampler.tsv')
    ts = get_stratified_gene_usage_frequency(ts = ts, replace = True) 

    # Create a list of Data.Frames each grouped by v and j gene usage 
    gs = [g for i,g in d.groupby(['v_b_gene','j_b_gene'])]
    fs = list()
    for sim in range(n_runs):
        print(sim)
        # Use multiple cpus to speed this up, the actual work is being done in the routine.py script
        y2 = parmap.map(vj_tr, gs, ts = ts, organism = organism, chains = chains, pm_processes = n_cpus, pm_pbar = True)
        y_df = pd.concat(y2)
        fs.append(f"{f}_radii_estimate_{sim}.tsv")
        y_df.to_csv(os.path.join(pout, f"{f}_radii_estimate_{sim}.tsv"), sep = "\t", index = False)

    # After thsi stage all of the computation is done. 
    # we just reread all of the saved file and compute summary statistics 
    # on radii estimates per clone.

    #fs = [f for f in os.listdir(pout) if f.find('radii_estimate_') != -1]
    print(fs)
    dfs = [pd.read_csv(os.path.join(pout,f), sep= "\t") for f in fs]
    dfall = pd.concat(dfs)

    x = dfall[['cdr3_b_aa','v_b_gene','j_b_gene','radi']].\
        groupby(['cdr3_b_aa','v_b_gene','j_b_gene']).\
        agg(['min','mean','median','max','std']).\
        reset_index()

    mi  = x.columns
    ind = pd.Index([e[0] +"_"+ e[1] if e[0] not in ['cdr3_b_aa','v_b_gene','j_b_gene'] else e[0] for e in mi.tolist() ])
    x.columns = ind
    x['radi_range']=x['radi_max']-x['radi_min']
    x.to_csv(os.path.join(pout,f'{f}_radii_estimates_across_{n_runs}_runs.tsv'),sep = "\t", index = False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Estimate radius per clone using dense V-J gene matched synthetic neighborhoos')
    parser.add_argument('-p',        help='path to input file', required=True)
    parser.add_argument('-pout',     help='path to store outputs', required=True)
    parser.add_argument('-f',        help='filename clones file', required=True)
    parser.add_argument('-n_cpus',   help='availale cpus to distribute', required=True)
    parser.add_argument('-n_runs',   help='number of times to estimate ideal radius', required=True)
    parser.add_argument('-organism', help='human or mouse', required=True)
    parser.add_argument('-chains',   help='comma separated list of beta  beta,alpha or alpha' , required=False, default = "beta")  
    args = parser.parse_args()
    print(vars(parser.parse_args()))
    
    main(
        p        = args.p, #"data"
        pout     = args.pout, #"output"
        f        = args.f, #"dash_human.csv"
        n_cpus   = int(args.n_cpus),#6
        n_runs   = int(args.n_runs), #10
        organism = args.organism ,#"human"
        chains   = args.chains.split(",") )
    """
    python radii_estimation.py -p data -pout output -f dash_human.csv -n_cpus 2 -n_runs 10 -organism human -chains beta
    """
