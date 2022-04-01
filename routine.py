"""
routine.py - this is a function that is called by radii estimation.
It is in a separate file so that it can be parmapped across many cpus
"""
from tcrdist.repertoire import TCRrep
from tcrdist.public import _neighbors_fixed_radius, _K_neighbors_fixed_radius
from tcrdist.background import _synthesize_mouse_alpha_vj_background
import numpy as np
import pandas as pd
from tcrdist.ecdf import distance_ecdf

def vj_tr(
    df, 
    ts,
    organism,
    ctrl_bkgd = 1E-6, 
    thrs = [x for x in range(0,50,2)],
    chain = 'beta'):
    """
    df : DataFrame 
        Clones all of the same V and J gene usage
    
    ts : TCRsampler 
    
    ctrl_bkgd : float
        1E-6 
    
    thrs : list
        list of TCRdist thresholds to consider e.g., [x for x in range(0,50,2)], 
    
    organism : str
         "human" or "mouse"
    
    chains : list
        ["beta"], ["alpha"], ["alpha","beta"]
    """
    chain_cols = { 'alpha' : ['cdr3_a_aa','v_a_gene','j_a_gene'],
                    'beta' : ['cdr3_b_aa','v_b_gene','j_b_gene']}
    pull_cols = chain_cols.get(chain)

    try:
        """compute_optimal_radii for a subset of TRBV/TRBJ"""
        i = (df[pull_cols[1]].iloc[0], df[pull_cols[2]].iloc[0])
        print(i)
        tr = TCRrep(   
                cell_df =df,
                organism = organism,
                chains = [chain],
                deduplicate = False, 
                compute_distances = False)
        
        if organism == 'mouse' and chain == 'alpha':
            df_vj_background = _synthesize_mouse_alpha_vj_background(ts = ts, df = tr.clone_df)
        else:
            df_vj_background = tr.synthesize_vj_matched_background(ts = ts, chain = chain)

        tr_bkgd = TCRrep(
            cell_df = df_vj_background,
            organism = organism, 
            chains = [chain], 
            compute_distances = False)
        
        tr.compute_rect_distances(df = tr.clone_df, df2 = tr_bkgd.clone_df)
        try:
            inv_weight = ts.v_occur_freq.get(i[0]) * ts.j_occur_freq.get(i[1])

        except TypeError:
            inv_weight = .01
        pV = ts.v_occur_freq.get(i[0])
        pJ = ts.j_occur_freq.get(i[1])
        if chain == 'alpha':
            mat = tr.rw_alpha
        else:
            mat = tr.rw_beta

        ws = [inv_weight] * mat.shape[1]
        thresholds, ecdfs2 = distance_ecdf(pwrect =mat, 
        thresholds = thrs, 
        weights= ws, 
        pseudo_count=0, 
        skip_diag = False, 
        absolute_weight = True)
        
        ecdfs2 = [pd.Series(x, index = thresholds) for x in ecdfs2]
        max_radi2 = [x[x<=ctrl_bkgd].last_valid_index() for x in ecdfs2]
        # TO SOLVE NOTABLE BUG IN THE ABOVE LINE!! If a radius is None (the next line will fail, thus set Nones to 0.
        max_radi2 = [x if (x is not None) else 0 for x in max_radi2]
        print(max_radi2)
        tr.clone_df['radi'] =  max_radi2 
        tr.clone_df['pV']= pV
        tr.clone_df['pJ']= pJ
        tr.clone_df['weight'] = inv_weight 

        for t in thrs:
            tr.clone_df[f'vj_hits_{t}'] = _K_neighbors_fixed_radius(mat, t)

        tr.clone_df['N'] = tr_bkgd.clone_df.shape[0]
        return tr.clone_df
    except:
        print(f"FAILED ON A PARTICUALR V-J SUBGROUP {i}")
        print(df)
        return df