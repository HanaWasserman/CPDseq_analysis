import pandas as pd
import pybedtools
import progressbar
import os
import pickle
import numpy as np
import argparse

parser = argparse.ArgumentParser(prog='ReorderBeds',
                                 description='Pre-processes bed files for bedtools steps')
parser.add_argument('data', nargs='+')
parser.add_argument('-p', '--data_path')
parser.add_argument('-f', '--path_to_fa')
args = parser.parse_args()

def main(dat, strnd):

    f_path = f"{d_path}/{dat}_{strnd}"
    print(f'Preparing: {f_path}')
    df = pd.read_table(f'{f_path}.bed', header=None, usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'strand'])
    df['name'], df['score'] = '', ''
    df_reordered = df[['chrom', 'start', 'end', 'name', 'score', 'strand']].copy()
    df_reordered_bed = pybedtools.BedTool.from_dataframe(df_reordered)
    df_reordered_bed.saveas(f'{f_path}_reordered.bed')


#d_path = f"CPDdata"  # Path to CPD-seq data
#d_files = ['WT_CSB_8J_24hr_CPD_1bp_sorted','WT_CSB_noUV_S32_CPD_1bp_sorted','XPC_6J_0hr_S19_CPD_1bp_sorted','XPC_6J_6hr_S20_CPD_1bp_sorted','XPC_6J_24hr_S21_CPD_1bp_sorted','WT_CSB_6J_0hr_CPD_1bp_sorted','WT_CSB_6J_6hr_CPD_1bp_sorted','WT_CSB_6J_24hr_CPD_1bp_sorted']  # Name of CPD-seq experiment XPC_12J_NakedDNA_S22_CPD_1bp_sorted WT_CSB_8J_0hr_CPD_1bp_sorted
#fasta = pybedtools.BedTool("../genomes/hg19/hg19.fa")  # Path to reference genome

d_path = args.data_path
d_files args.data
fasta = pybedtools.BedTool(args.path_to_fa) 

for df in d_files:
    main(df, 'plus')
    main(df, 'minus')
    print("")
