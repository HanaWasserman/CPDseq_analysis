import pandas as pd
import pybedtools
import progressbar
import os
import pickle
import numpy as np
import argparse

parser = argparse.ArgumentParser(prog='ProcessBackgroundData',
                                 description='Pre-processes CPD-seq dataset to consolidate damage signal')
parser.add_argument('data', nargs='+')
parser.add_argument('-tf', '--transcription_factors')
parser.add_argument('-p', '--path_to_data')
parser.add_argument('-g', '--path_to_fa')
args = parser.parse_args()

def main(bg, bf, d_file):

    def prep(o_txt):
        """
        BedTools intersect CPD-seq experiment data with background region bed file to get relevant subset of CPD-seq
        signal
        and consolidate signal
        String:param o_txt: Strand orientation
        String:param bg: Region of interest name
        String:param bed_file: Region of interest windows bed file path
        Dataframe:return: Consolidated CPD-seq signal in regions of interest windows
        """

        f_path = f"{fold_txt}/{bg}_intersect_{d_file}_{o_txt}strand"
        print(f'Preparing: {bg}_intersect_{d_file}_{o_txt}')
        os.system(f"bedtools intersect -u -a {d_path}/{d_file}_{o_txt}_reordered.bed -b TFs/{bf}.bed > {f_path}.bed")  # OUTPUT1, OUTPUT2
        df = pd.read_table(f'{f_path}.bed', header=None, usecols=[0, 1, 2], names=['chr', 'st', 'end'])
        """Iterates through CPD-seq experiment data and consolidates CPD-seq signal from two, 1-nucleotide rows per damage 
        to one 2-nucleotide row per damage in a new dataframe"""
        df["counts"] = 1
        df = df.groupby(["chr", "st", "end"])["counts"].count().reset_index()
        agg_d = list()
        df_end = df.shape[0] - 1
        for i in progressbar.progressbar(range(0, df_end)):
            chr, st, end, counts = df.loc[i]
            chr2, st2, end2, counts2 = df.loc[i + 1]
            if end == st2 and i < df_end and counts != 0:
                if counts <= counts2:
                    """Consecutive rows where end_(n) and start_(n+1) are equal should belong to same damage and can 
                    be consolidated"""
                    agg_d.append((chr, st, end, counts))
                    """Subtract n damage count from n+1 damage count to prevent overcounting"""
                    df.at[i + 1, 'counts'] = df.at[i + 1, 'counts'] - df.at[i, 'counts']
                elif counts > counts2:
                    """Consolidates damages in cases where there are two consecutive damages and the first damage has 
                    higher signal than the following damage"""
                    agg_d.append((chr, st, end, counts2))
                    df.at[i + 1, 'counts'] = df.at[i, 'counts'] - df.at[i + 1, 'counts']
            i = i + 1
        agg_df = pd.DataFrame(agg_d)
        """Add strand orientation"""
        if o_txt == 'plus':
            agg_df['ori'] = '+'
        elif o_txt == 'minus':
            agg_df['ori'] = '-'
        agg_df = agg_df.rename(columns={0: "chr", 1: "st", 2: "end", 3: "counts", 4: "ori"})
        return agg_df[['chr', 'st', 'end', 'counts']]

    prep_df_plus = prep('plus')
    prep_df_minus = prep('minus')

    def adjust(df, o_txt, st_l, en_l, s_bool, f_path=f'{fold_txt}/{bg}_intersect_{d_file}'):
        print(f'Adjusting: {bg}_{d_file}_{o_txt}')
        df["ori"] = "-"
        df["misc"] = ""
        df["st"] += st_l
        df["st"] = df["st"].astype(np.int64)  #
        df["end"] += en_l
        df["end"] = df["end"].astype(np.int64)  #
        df_ = df[['chr', 'st', 'end', 'counts', 'misc', 'ori']].copy()
        tf_bed = pybedtools.BedTool.from_dataframe(df_)
        tf_bed_seq = tf_bed.sequence(fi=fasta, bedOut=True, s=s_bool)
        tf_bed_seq.save_seqs(f"{f_path}_{o_txt}_seq_tmp")
        seqs = pd.read_table(f"{f_path}_{o_txt}_seq_tmp", header=None)
        df_["seq"] = seqs
        if o_txt == 'plus':
            df_["ori"] = "+"
        if os.path.exists(f"{f_path}_{o_txt}_seq_tmp"):
            os.remove(f"{f_path}_{o_txt}_seq_tmp")
        return df_[['chr', 'st', 'end', 'counts', 'ori', 'seq']]

    adj_df_plus = adjust(prep_df_plus, 'plus', -1, 2, True)
    adj_df_minus = adjust(prep_df_minus, 'minus', -1, 2, False)

    def filter(df, o_txt, f_path):
        """
        Verifies CPD signal maps to 4mers containing dipyrimdine sequences and filters for dipy and Cdipy conditions
        Dataframe:param df: Consolidated CPD-seq data with mapped fasta sequence
        String:param o_txt: Strand orientation
        String:param f_path: Destination path to save files to
        String:param dipy_txt: All dipyrimidines (dipy) or cytosine-containing dipyrimdine condition (Cdipy)
        """
        all_n = df.shape[0]
        df["seq"] = df["seq"].str.upper()
        df["YY"] = df["seq"].str[1:3]
        df_dipy = df[(df["YY"] == "CC") | (df["YY"] == "CT") | (df["YY"] == "TC") | (df["YY"] == "TT")]
        df_dipy_ = df_dipy[["chr", "st", "end", "counts", "seq", "ori"]].copy()
        dipy_n = str(df_dipy_.shape[0] / all_n)
        tf_dipy_bed = pybedtools.BedTool.from_dataframe(df_dipy_)
        tf_dipy_bed.saveas(f"{f_path}_dipy_proc_{o_txt}.bed")  # OUTPUT3, OUTPUT4, OUTPUT5, OUTPUT6
        print(f"{dipy_n}% dipy")

    print(f'Filtering: {bg}_{d_file}')
    filter(adj_df_plus, 'plus', f'{fold_txt}/{bg}_intersect_{d_file}')
    filter(adj_df_minus, 'minus', f'{fold_txt}/{bg}_intersect_{d_file}')

    os.remove(f"{fold_txt}/{bg}_intersect_{d_file}_minusstrand.bed")
    os.remove(f"{fold_txt}/{bg}_intersect_{d_file}_plusstrand.bed")


#d_path = f"CPDdata"  # Path to CPD-seq data
d_path = args.path_to_data
fold_txt = 'CPDseq_data/background'  # Path output will be saved to
#d_files = ['XPC_12J_NakedDNA_S22_CPD_1bp_sorted']  # Name of CPD-seq experiment XPC_12J_NakedDNA_S22_CPD_1bp_sorted WT_CSB_8J_0hr_CPD_1bp_sorted
d_files = args.data
#fasta = pybedtools.BedTool("../genomes/hg19/hg19.fa")  # Path to reference genome
fasta = pybedtools.BedTool(args.path_to_fa)

if not os.path.exists(fold_txt):
    os.system(f"mkdir {fold_txt}")
motifs = pd.read_csv(f'TFs/{args.transcription_factors}', usecols=['Width', 'Name_str'])
motifs_r = motifs.to_records(index=False)
for df in d_files:
    for w, n in motifs_r:
        main(f'{n}_flanks_unmerged_200to20', f'{n}_flanks_unmerged_200to20_seq', df)
        print("")
    print("")
