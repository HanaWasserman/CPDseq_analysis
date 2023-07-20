#!/usr/local/bin/python3

import numpy as np
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(prog='calcPredicted',
                                 description='Calculate the predicted CPD signal in TFBS region based on background model')
parser.add_argument('data', nargs='+')
parser.add_argument('-tf', '--transcription_factors')
#args = parser.parse_args(['XPC_12J_NakedDNA_S22_CPD_1bp_sorted', '-tf', 'motif_annotations_subsubset.csv'])
args = parser.parse_args()


def main(tf_txt, bkgrd_opt, wind_sz, roi_wind_sz, dat_wind_sz, Cdip_opt):

    bkgrd = dict.fromkeys(d_files)
    yys, frmrs = dict.fromkeys(['same', 'opp']), dict.fromkeys(['same', 'opp'])

    if Cdip_opt:
        dip_txt = 'Cdipy'
    elif not Cdip_opt:
        dip_txt = 'dipy'

    def get_yys(strnd_type, ori_txt):
        yys[strnd_type] = {}
        yy_data_comp = pd.read_csv(f"Analysis/{tf_txt}_{roi_wind_sz}_seq_dipy_4mer_counts_{ori_txt}_ROI.csv", usecols=[1, 2, 3])
        if Cdip_opt:
            yy_data_comp = yy_data_comp.loc[yy_data_comp['seq'].str[1:3]!='TT']
        yy_data = yy_data_comp[["pos", "counts"]].groupby(["pos"])["counts"].sum().reset_index()
        yy_data["pos"] = np.arange(-wind_sz, yy_data.shape[0] - wind_sz)
        yys[strnd_type], frmrs[strnd_type] = yy_data, yy_data_comp

    get_yys('same', 'plus'), get_yys('opp', 'minus')


    def get_bkgrd(exp, stats, strnd_type):
        if Cdip_opt:
            stats = stats.loc[stats['seq'].str[1:3]!='TT']
        exp_cpds = pd.merge(stats, frmrs[strnd_type], on='seq')
        exp_cpds['mean_tot'], exp_cpds['var_tot'] = exp_cpds['counts'] * exp_cpds['mean'], exp_cpds['counts'] * exp_cpds['mean']
        exp_cpds[['seq', 'pos', 'counts', 'mean_tot', 'var_tot']].to_csv(f"Analysis/{tf_txt}_{exp}_{bkgrd_opt}_{dip_txt}_predicted_{strnd_type}_{dat_wind_sz}_all_seqs.csv", index=False)
        exp_cpds_df = exp_cpds.groupby(['pos'], as_index=False).sum(numeric_only=True)
        exp_cpds_df['std_dev'], exp_cpds_df["pos"] = np.sqrt(exp_cpds_df['var_tot']), np.arange(-wind_sz, exp_cpds_df.shape[0] - wind_sz)
        bkgrd[exp][strnd_type] = exp_cpds_df[3:-4]
        bkgrd[exp][strnd_type][['pos','counts','mean_tot','var_tot','std_dev']].to_csv(f"Analysis/{tf_txt}_{exp}_{bkgrd_opt}_{dip_txt}_predicted_{strnd_type}_{dat_wind_sz}.csv", index=False)  # OUTPUT3

    for b in bkgrd:
        bkgrd[b] = {}
        cpd_stats = pd.read_csv(f'CPDseq_data/background/{bkgrd_opt}_intersect_{b}_cpd_4mers_BKGRND.csv', usecols=[0, 3, 4])
        get_bkgrd(b, cpd_stats, "same"), get_bkgrd(b, cpd_stats, "opp")


fold_txt = 'Analysis'
if not os.path.exists(fold_txt):
    os.system(f"mkdir {fold_txt}")
motifs = pd.read_csv(f'TFs/{args.transcription_factors}', usecols=['Width','Name_str'])
motifs_r = motifs.to_records(index=False)
d_files = args.data
for w, ns in motifs_r:
    print(f'Processing: {ns}')
    #main(ns, f'{ns}_flanks_unmerged_200to20', 200.5, '202', '200', False)
    #main(ns, f'{ns}_flanks_unmerged_200to20', 200.5, '202', '200', True)
    main(ns, f'SKCM', 200.5, '202', '200', False)
    main(ns, f'SKCM', 200.5, '202', '200', True)
