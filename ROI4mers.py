#!/usr/local/bin/python3

import os
import pandas as pd
import progressbar
import argparse

parser = argparse.ArgumentParser(prog='ROI4mers',
                                 description='Count 4mers at each position in TFBS region')
parser.add_argument('data', nargs='+')
parser.add_argument('-tf', '--transcription_factors')
args = parser.parse_args()

def get_tf(bed_file):
    """
    Load region of interest bed file and determine full window size
    String:param bed_file: Region of interest bed file path
    Dataframe:return: Dataframe containing windows for region of interest
    Int:return: Window size
    """
    bs = pd.read_table(f'{bed_file}.bed', header=None, usecols=[0, 1, 2, 3, 4], names=['chr', 'start', 'end', 'seq', 'rv_seq'])
    bs[["start", "end"]] = bs[["start", "end"]].to_numpy(dtype='int')
    ws = bs["end"][0] - bs["start"][0]
    return bs, ws

def init4mers():
    yys = ('CC', 'CT', 'TC', 'TT')
    ns = ('A', 'C', 'T', 'G')
    qn_l = list()
    for qn in yys:
        for n1 in ns:
            for n2 in ns:
                qn_l.append((f'{n1}{qn}{n2}'))
    qn = dict.fromkeys(qn_l, 0)
    for q in qn:
        qn[q] = []
    return qn

def main(tf_txt, bf):
    """
    Get distribution of all dipyrimidine-containing 4mers occurring across TFBS windows sequences
    String:param tf_txt: Name of TF of interest
    Saved Output:
    7-8. 4mer_counts_by_relative_position.csv (for both the + and - strand sequences)
    """

    bs, window_sz = get_tf(f'TFs/{bf}')

    def count4(windows, ori):
        """
        Dataframe:param windows: TFBS windows
        Boolean:param tf_o: True means TFBS is on the + strand and False means TFBS is on - strand
        List:return: List of dipyrimdine-containing 4mer sequences and their relative position in the TFBS window
        """
        yys = ['TT', 'CT', 'TC', 'CC']
        frs = init4mers()
        for start, end, seq in progressbar.progressbar(windows.iloc):
            assert end > start
            if end-start > window_sz:
                continue
            seq_l = len(seq)
            for i in range(seq_l - 3):
                seq_n = seq[i:i + 4]
                if 'N' in seq_n:
                    print('hi')
                    continue
                if seq_n[1:3] in yys:
                    if ori:
                        frs[seq_n].append(i + 1)
                    elif not ori:
                        frs[seq_n].append(seq_l - i - 2)
        return frs

    print(f'Counting {tf_txt} 4mers...')
    p_4mers = count4(bs[['start', 'end', 'seq']], True)
    m_4mers = count4(bs[['start', 'end', 'rv_seq']], False)

    def prep_4mer_counts(cnts):
        """
        Aggregates counts for each 4mer by relative position in the TFBS windows
        List:param cnts: List of dipyrimdine-containing 4mer sequences and their relative position in the TFBS windows
        Dataframe:return: Dataframe of 4mer counts by relative position
        """

        cnts_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in cnts.items()])).melt(var_name='seq', value_name='pos')
        cnts_df = cnts_df.dropna()
        cnts_df['counts'] = 1
        cnts_df = cnts_df.groupby(['pos', 'seq'], as_index=False).sum()
        cnts_df['pos'] = cnts_df['pos'].astype(int)
        cnts_df_cdip = cnts_df[cnts_df['seq'].str[1:3] != 'TT']
        return cnts_df, cnts_df_cdip

    p_4mer_counts_dipy, p_4mer_counts_cdipy = prep_4mer_counts(p_4mers)
    m_4mer_counts_dipy, m_4mer_counts_cdipy  = prep_4mer_counts(m_4mers)
    p_4mer_counts_dipy.to_csv(f"{fold_txt}/{bf}_dipy_4mer_counts_plus_ROI.csv") #OUTPUT7
    m_4mer_counts_dipy.to_csv(f"{fold_txt}/{bf}_dipy_4mer_counts_minus_ROI.csv") #OUTPUT8

fold_txt = 'Analysis'
if not os.path.exists(fold_txt):
    os.system(f"mkdir {fold_txt}")
d_path = f"CPDseq_data"
#d_files = ['WT_CSB_8J_0hr_CPD_1bp_sorted','XPC_12J_NakedDNA_S22_CPD_1bp_sorted'] #XPC_12J_NakedDNA_S22_CPD_1bp_sorted WT_CSB_8J_0hr_CPD_1bp_sorted
d_files = args.data

motifs = pd.read_csv(f'TFs/{args.transcription_factors}', usecols=['Width','Name_str'])
motifs_r = motifs.to_records(index=False)
for w, n in motifs_r:
    main(n, f'{n}_202_seq')
