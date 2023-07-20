#!/usr/local/bin/python3

import os
import numpy as np
import pandas as pd
import progressbar
import argparse

parser = argparse.ArgumentParser(prog='countCPDs',
                                 description='Aggregates all CPD counts in TFBS windows')
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
    bs = pd.read_table(f'{bed_file}.bed', header=None, usecols=[0, 1, 2, 3, 4, 5], names=['chr', 'start', 'end', 'seq', 'rv_seq', 'ori'])
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

def main(tf_txt, nuc_txt, d_file, bed_file, dat_wind_sz):
    """
    Aggregates all CPD counts in TFBS windows
    String:param tf_txt: Name of TF of interest
    String:param nuc_txt: Analyze all CPDs or only cytosine-containing CPDs
    String:param d_file: Pre-processed CPD data file name
    Saved Output:
    1-2. CPD_counts_per_position_in_TFBS_windows (for both CPDS on same and opposite strand)
    3-6. Info_for_each_CPD_occuring_in_TFBS_windows (for all strand combos for CPD and TFBS location)
    """
    f_txt = f"{d_path}/{bed_file}_intersect_{d_file}_{nuc_txt}_proc"
    chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    bs1, window_sz_1 = get_tf(f'TFs/{bed_file}')

    def get_cpds_from_file(fname):
        """
        Get pre-processed CPD data
        String:param fname: CPD data file
        Dataframe:return: Dataframe of CPD data (with pUC19 omitted)
        """
        print("Loading " + fname)
        cpd = pd.read_table(fname, header=None, usecols=[0, 1, 3, 4], names=['chr','start','counts','seq'])
        return cpd

    def prep_cpds_chrom(all_cpds, chromo):
        """
        Filter CPD data for chromosome of interest
        Dataframe:param all_cpds: CPD data
        String:param chromo: Chromsome of interest
        Dataframe:return: Dataframe of CPD data for chromosome of interest
        """
        chrom_cpds = all_cpds[all_cpds["chr"] == chromo]
        return chrom_cpds

    def count_cpds(windows, wind_sz, cpd_df, tf_o):
        """
        Iterate through each TFBS window and count / record CPDs for each position in the window
        Dataframe:param bs: TFBS windows
        Int:param wind_sz: Window size
        Dataframe:param cpd_df: Dataframe of pre-processed CPD data filtered for chromosome of interest
        String:param chromo: Chromosome of interest
        Boolean:param tf_o: True means TFBS is on the + strand and False means TFBS is on - strand
        Array:return: Array of total CPD counts per position in TFBS windows
        Dataframe:return: Dataframe of CPD info for all CPDs that occur in the TFBS windows
        """
        cpd_df_len = cpd_df.shape[0]
        d_counts = np.zeros(wind_sz, dtype='int')
        c_counts = np.zeros(wind_sz, dtype='int')
        d_seqs = list()
        i = 0
        for bs_start, bs_end, bs_seq in progressbar.progressbar(windows.iloc):
            assert bs_end > bs_start
            if bs_end-bs_start > wind_sz:
                continue
            while True:
                if i > cpd_df_len-1:
                    break
                if cpd_df[i]['start'] > bs_start:
                    if i > 0:
                        i = i - 1
                        continue
                    else:
                        break
                else:
                    break
            while i < cpd_df_len:
                cpd, count, seq = cpd_df[i]
                if cpd <= bs_start:
                    i = i + 1
                    continue
                if cpd >= bs_end - 1:
                    break
                if tf_o:
                    d_counts[cpd - bs_start] = d_counts[cpd - bs_start] + count
                    d_seqs.append((cpd - bs_start, seq, count))
                    if seq != 'TT':
                        c_counts[cpd - bs_start] = c_counts[cpd - bs_start] + count
                elif not tf_o:
                    d_counts[bs_end - cpd - 2] = d_counts[bs_end - cpd - 2] + count
                    d_seqs.append((bs_end - cpd - 2, seq, count))
                    if seq != 'TT':
                        c_counts[bs_end - cpd - 2] = c_counts[bs_end - cpd - 2] + count
                else:
                    print(f'CPD did not match on either strand: {cpd} {seq}')
                i = i + 1
        return d_counts, c_counts, d_seqs

    p_cpd = get_cpds_from_file(f'{f_txt}_plus.bed')
    m_cpd = get_cpds_from_file(f'{f_txt}_minus.bed')

    """Initialize empty output array for total CPD counts and empty output Dataframe for complete CPD data.
    pp = CPD and TFBS is on + strand
    pm = CPD is on + strand and TFBS is on - strand
    mp = CPD is on - strand and TFBS is on + strand
    mm = CPD and TFBS is on - strand"""
    all_psame_counts = np.zeros(window_sz_1, dtype='int')
    all_msame_counts = np.zeros(window_sz_1, dtype='int')
    all_popp_counts = np.zeros(window_sz_1, dtype='int')
    all_mopp_counts = np.zeros(window_sz_1, dtype='int')
    all_psame_counts_c = np.zeros(window_sz_1, dtype='int')
    all_msame_counts_c = np.zeros(window_sz_1, dtype='int')
    all_popp_counts_c = np.zeros(window_sz_1, dtype='int')
    all_mopp_counts_c = np.zeros(window_sz_1, dtype='int')
    all_same_seqs = pd.DataFrame()
    all_opp_seqs = pd.DataFrame()

    print('Counting CPDs...')
    for c in chroms:
        """Iterate through all chromsomes"""
        p_cpd_chr = prep_cpds_chrom(p_cpd, c)
        m_cpd_chr = prep_cpds_chrom(m_cpd, c)
        p_windows = bs1.loc[(bs1["chr"] == c) & (bs1['ori'] == '+'), ['start', 'end', 'seq']]
        m_windows = bs1.loc[(bs1["chr"] == c) & (bs1['ori'] == '-'), ['start', 'end', 'seq']]

        p_same_counts, p_same_counts_c, p_same_seqs = count_cpds(p_windows, window_sz_1, m_cpd_chr[["start", "counts", "seq"]].to_records(index=False), True)
        m_same_counts, m_same_counts_c, m_same_seqs = count_cpds(m_windows, window_sz_1, p_cpd_chr[["start", "counts", "seq"]].to_records(index=False), False) #False
        p_opp_counts, p_opp_counts_c, p_opp_seqs = count_cpds(p_windows, window_sz_1, p_cpd_chr[["start", "counts", "seq"]].to_records(index=False), True)
        m_opp_counts, m_opp_counts_c, m_opp_seqs = count_cpds(m_windows, window_sz_1, m_cpd_chr[["start", "counts", "seq"]].to_records(index=False), False) #False

        """Consolidate CPD count arrays across chromosomes"""
        all_psame_counts = all_psame_counts + p_same_counts
        all_msame_counts = all_msame_counts + m_same_counts
        all_popp_counts = all_popp_counts + p_opp_counts
        all_mopp_counts = all_mopp_counts + m_opp_counts

        all_psame_counts_c = all_psame_counts_c + p_same_counts_c
        all_msame_counts_c = all_msame_counts_c + m_same_counts_c
        all_popp_counts_c = all_popp_counts_c + p_opp_counts_c
        all_mopp_counts_c = all_mopp_counts_c + m_opp_counts_c

        """Consolidate CPD info dataframes across chromosomes"""
        def concat_seqs(s, a_df):
            s_df = pd.DataFrame(s)
            tot_df = pd.concat([a_df, s_df])
            return tot_df
        all_same_seqs = concat_seqs(p_same_seqs, all_same_seqs)
        all_same_seqs = concat_seqs(m_same_seqs, all_same_seqs)
        all_opp_seqs = concat_seqs(p_opp_seqs, all_opp_seqs)
        all_opp_seqs = concat_seqs(m_opp_seqs, all_opp_seqs)

    all_same_seqs_c = all_same_seqs.loc[all_same_seqs[1] != 'TT']
    all_opp_seqs_c = all_opp_seqs.loc[all_opp_seqs[1] != 'TT']
    all_same_counts = all_psame_counts + all_msame_counts
    all_opp_counts = all_popp_counts + all_mopp_counts
    all_same_counts_c = all_psame_counts_c + all_msame_counts_c
    all_opp_counts_c = all_popp_counts_c + all_mopp_counts_c
    #np.savetxt(f"{fold_txt}/{tf_txt}_{d_file}_dipy_psame_{dat_wind_sz}.txt", all_psame_counts, delimiter=",")  # OUTPUT1
    #np.savetxt(f"{fold_txt}/{tf_txt}_{d_file}_dipy_msame_{dat_wind_sz}.txt", all_msame_counts, delimiter=",")  # OUTPUT2
    #np.savetxt(f"{fold_txt}/{tf_txt}_{d_file}_dipy_popp_{dat_wind_sz}.txt", all_popp_counts, delimiter=",")  # OUTPUT1
    #np.savetxt(f"{fold_txt}/{tf_txt}_{d_file}_dipy_mopp_{dat_wind_sz}.txt", all_mopp_counts, delimiter=",")  # OUTPUT2
    np.savetxt(f"{fold_txt}/{tf_txt}_{d_file}_dipy_same_{dat_wind_sz}.txt", all_same_counts, delimiter=",") #OUTPUT1
    np.savetxt(f"{fold_txt}/{tf_txt}_{d_file}_dipy_opp_{dat_wind_sz}.txt", all_opp_counts, delimiter=",") #OUTPUT2
    np.savetxt(f"{fold_txt}/{tf_txt}_{d_file}_Cdipy_same_{dat_wind_sz}.txt", all_same_counts_c, delimiter=",")  # OUTPUT1
    np.savetxt(f"{fold_txt}/{tf_txt}_{d_file}_Cdipy_opp_{dat_wind_sz}.txt", all_opp_counts_c, delimiter=",")  # OUTPUT2
    all_same_seqs.to_csv(f"{fold_txt}/{tf_txt}_{d_file}_dipy_same_seqs_{dat_wind_sz}.csv") #OUTPUT3
    all_opp_seqs.to_csv(f"{fold_txt}/{tf_txt}_{d_file}_dipy_opp_seqs_{dat_wind_sz}.csv") #OUTPUT4
    all_same_seqs_c.to_csv(f"{fold_txt}/{tf_txt}_{d_file}_Cdipy_same_seqs_{dat_wind_sz}.csv")  # OUTPUT3
    all_opp_seqs_c.to_csv(f"{fold_txt}/{tf_txt}_{d_file}_Cdipy_opp_seqs_{dat_wind_sz}.csv")  # OUTPUT4


fold_txt = 'Analysis'
if not os.path.exists(fold_txt):
    os.system(f"mkdir {fold_txt}")
d_path = f"CPDseq_data"
#d_files = ['WT_CSB_8J_0hr_CPD_1bp_sorted','XPC_12J_NakedDNA_S22_CPD_1bp_sorted'] #XPC_12J_NakedDNA_S22_CPD_1bp_sorted
d_files = args.data

#motifs = pd.read_csv('TFs/motif_annotations_proc_filt_count.csv', usecols=['Width','Name_str'])
motifs = pd.read_csv(f'TFs/{args.transcription_factors}', usecols=['Width','Name_str'])
motifs_r = motifs.to_records(index=False)
for w, n in motifs_r:
    for df in d_files:
        print(n)
        main(n, "dipy", df, f'{n}_200_seq', '200')
