import pandas as pd
import pybedtools
from pybedtools import BedTool
import progressbar
from weblogo import *
import argparse
parser = argparse.ArgumentParser(prog='PrepareFlanks',
                                 description='Prepare TFBS bed files to define flanking regions for background modeling')
parser.add_argument('-g', '--path_to_genome_files')
parser.add_argument('-tf', '--transcription_factors')
args = parser.parse_args()


def main(fold_txt, tf_txt, gen, shft_sz, padding_sz, merge):

    gen_fa_fai, gen_fa = f'{gen}/hg19.fa.fai', f'{gen}/hg19.fa'
    if not os.path.exists(fold_txt):
        os.system(f"mkdir {fold_txt}")

    def slop_it():

        """Option to force strandedness when using bedtools getfasta"""
        if merge:
            s_opt = ' -s'
        elif not merge:
            s_opt = ''

        os.system(f"bedtools getfasta -fi {gen_fa} -bed {fold_txt}/{tf_txt}_flanks_tmp.bed -fo {fold_txt}/{tf_txt}_flanks_tmp_seq{s_opt} -bedOut")

        tf_bed = BedTool(f'{fold_txt}/{tf_txt}_flanks_tmp.bed')
        bed_df = pybedtools.BedTool.to_dataframe(tf_bed)
        seqs = pd.read_table(f"{fold_txt}/{tf_txt}_flanks_tmp_seq", header=None)
        tf_bed_df = bed_df[['chrom', 'start', 'end']].copy()
        if merge:
            mrg_txt = 'merged'
            tf_bed_df['strand'] = bed_df['name']
        elif not merge:
            mrg_txt = 'unmerged'
            tf_bed_df['strand'] = bed_df['strand']
        tf_bed_df['seq'] = seqs[0].str.upper()
        """Add reverse complement of fasta sequence to bedfile dataframe"""
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        tf_bed_df['seq_rv'] = ""
        for i in progressbar.progressbar(range(tf_bed_df.shape[0])):
            tf_bed_df['seq_rv'][i] = "".join(complement.get(base, base) for base in reversed(tf_bed_df['seq'][i]))

        tf_bed_df_reordered = tf_bed_df[['chrom', 'start', 'end', 'seq', 'seq_rv', 'strand']].copy()
        bed_final = pybedtools.BedTool.from_dataframe(tf_bed_df_reordered)

        bed_final.saveas(f"{fold_txt}/{tf_txt}_flanks_{mrg_txt}_{shft_sz}to{padding_sz}_seq.bed") #OUTPUT2

    def prep_flanks():
        """
        Create bedfile of flanking regions around TFBS to be used for background modeling
        Saved Output:
        1. TF_flanks.bed (Required file for subsequent slop_it, then deleted)
        """
        os.system(f"bedtools slop -i {fold_txt}/{tf_txt}.bed -g {gen_fa_fai} -l {shft_sz} -r -{padding_sz} > {fold_txt}/{tf_txt}_flanks_left_tmp_slop.bed")
        os.system(f"bedtools slop -i {fold_txt}/{tf_txt}.bed -g {gen_fa_fai} -l -{padding_sz} -r {shft_sz} > {fold_txt}/{tf_txt}_flanks_right_tmp_slop.bed")
        os.system(f"cat {fold_txt}/{tf_txt}_flanks_left_tmp_slop.bed {fold_txt}/{tf_txt}_flanks_right_tmp_slop.bed > {fold_txt}/{tf_txt}_flanks_cat_tmp.bed")
        os.system(f"bedtools sort -i {fold_txt}/{tf_txt}_flanks_cat_tmp.bed > {fold_txt}/{tf_txt}_flanks_cat_tmp_sorted.bed")
        if merge:
            os.system(f"bedtools merge -i {fold_txt}/{tf_txt}_flanks_cat_tmp_sorted.bed -s -c 6 -o distinct > {fold_txt}/{tf_txt}_flanks_tmp.bed")
            os.remove(f"{fold_txt}/{tf_txt}_flanks_cat_tmp_sorted.bed")
            #tmp_bed = pd.read_table(f"{fold_txt}/{tf_txt}_flanks_tmp.bed", header=None, names=['chrom','start','end','strand'])
            #tmp_bed['name'], tmp_bed['score'] = "",""
        if not merge:
            os.system(f"mv {fold_txt}/{tf_txt}_flanks_cat_tmp_sorted.bed {fold_txt}/{tf_txt}_flanks_tmp.bed")
        os.remove(f"{fold_txt}/{tf_txt}_flanks_left_tmp_slop.bed")
        os.remove(f"{fold_txt}/{tf_txt}_flanks_right_tmp_slop.bed")
        os.remove(f"{fold_txt}/{tf_txt}_flanks_cat_tmp.bed")


    prep_flanks()
    slop_it()
    os.remove(f"{fold_txt}/{tf_txt}_flanks_tmp.bed")
    #os.remove(f"{fold_txt}/{tf_txt}_flanks_tmp_seq")


motifs = pd.read_csv(f'TFs/{args.transcription_factors}', usecols=['Name','Width','Name_str'])
motifs_r = motifs.to_records(index=False)
#g_path, f_path = '../genomes/hg19', 'TFs'
#archetype_site_Skin_intersect
g_path, f_path = args.path_to_genome_files, 'TFs'
for n, w, ns in motifs_r:
    print(f'Processing: {n}')
    main(f_path, ns, g_path, 200, 20, False)

#if not os.path.exists('TFs/SKCM_intersect_seq.bed'):
#    main('TFs', None, 'SKCM', None, '../genomes/hg19', None)
