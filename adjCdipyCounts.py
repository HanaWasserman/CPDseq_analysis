#!/usr/local/bin/python3


import pandas as pd
import os
import argparse
parser = argparse.ArgumentParser(prog='adjCdipyCounts',
                                 description='Adjust Cdipy signal for CtoT mutation counts')
parser.add_argument('data', nargs='+')
parser.add_argument('-tf', '--transcription_factors')
args = parser.parse_args()

def group_sum(df, grpby):
    return df.groupby(grpby, as_index=False).sum()

def main(tf_txt, dat_wind_sz, bkgrd_opt, exp_txt):

    def get_data(ori):
        dat_df = pd.read_csv(f"Analysis/{tf_txt}_{exp_txt}_Cdipy_{ori}_seqs_{dat_wind_sz}.csv", skiprows=[0], usecols=[1, 2, 3], names=['pos', 'seq', 'counts'])
        dat_df['counts'] = pd.to_numeric(dat_df['counts'])
        dat_df = group_sum(dat_df, ['pos', 'seq'])
        return dat_df[500:-500]

    dat_same = get_data('same')
    dat_opp = get_data('opp')

    def get_pred(ori):
        pred_df = pd.read_csv(f"Analysis/{tf_txt}_{exp_txt}_{bkgrd_opt}_Cdipy_predicted_{ori}_{dat_wind_sz}_all_seqs.csv", skiprows=[0], usecols=[0, 1, 3, 4], names=['old_seq','pos','counts','vars'])
        pred_df['seq'] = pred_df['old_seq'].str[1:3]
        pred_df = group_sum(pred_df, ['pos','seq'])
        return pred_df[['pos','seq','counts','vars']][500:-500]

    pred_same = get_pred("same")
    pred_opp = get_pred("opp")

    def cnts(df, m, pred, offset, title_txt):
        st = df.reset_index(drop=True)
        st["pos"] = st["pos"] + 0.5
        cc_seqs, i, lng = list(), 0, st.shape[0] - 1
        cols = []
        while i < lng:
            if pred:
                pos, seq, counts, vars = st.loc[i]
                cols = ['pos', 'seq', 'counts', 'vars']
            elif not pred:
                pos, seq, counts = st.loc[i]
                cols = ['pos', 'seq', 'counts']
            if seq in ["CT"]:
                st.at[i, 'pos'] = st.at[i, 'pos'] - m
            elif seq in ["CC"]:
                st.at[i, 'pos'] = st.at[i, 'pos'] - m
                if pred:
                    cc_seqs.append((pos + m, seq, counts, vars))
                elif not pred:
                    cc_seqs.append((pos+m, seq, counts))
            elif seq in ["TC"]:
                st.at[i, 'pos'] = st.at[i, 'pos'] + m
            i += 1
        cc_seqs_df = pd.DataFrame(cc_seqs, columns=cols)
        tot_st = group_sum(pd.concat([st, cc_seqs_df]), ['pos', 'seq'])
        tot_st['pos'] = tot_st['pos'] - dat_wind_sz - offset
        tot_st[1:-2].to_csv(f"{fold_txt}/{tf_txt}_{exp_txt}_Cdipy_{title_txt}_adj.csv", index=False)


    cnts(dat_same, 0.5, False, 0, f'same_seqs_{dat_wind_sz}')
    cnts(dat_opp, -0.5, False, 0, f'opp_seqs_{dat_wind_sz}')
    cnts(pred_same, 0.5, True, 2, f'{bkgrd_opt}_predicted_same_{dat_wind_sz}_all_seqs')
    cnts(pred_opp, -0.5, True, 3, f'{bkgrd_opt}_predicted_opp_{dat_wind_sz}_all_seqs')


fold_txt = 'Analysis'
if not os.path.exists(fold_txt):
    os.system(f"mkdir {fold_txt}")
motifs = pd.read_csv(f'TFs/{args.transcription_factors}', usecols=['Width','Name_str'])
datasets = args.data
motifs_r = motifs.to_records(index=False)
for tf_l, tf in motifs_r: #[(11, 'ETS_2'), (14, 'RFX_1')]
    print(tf)
    for d in datasets:
        main(tf, 200, f'{tf}_flanks_unmerged_200to20', d)

