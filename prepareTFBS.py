import pandas as pd
import pybedtools
from pybedtools import BedTool
import progressbar
from weblogo import *
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(prog='PrepareTFBS',
                                 description='Prepare TFBS bed files')
parser.add_argument('-g', '--path_to_genome_files')
parser.add_argument('-c', '--path_to_tf_clusters')
parser.add_argument('-tf', '--transcription_factors')
args = parser.parse_args()


def main(fold_txt, tf, tf_txt, arch_SKCM_path, gen, shft_sz):

    gen_fa_fai, gen_fa = f'{gen}/hg19.fa.fai', f'{gen}/hg19.fa'
    if not os.path.exists(fold_txt):
        os.system(f"mkdir {fold_txt}")

    """
    Subset binding site windows for TF of interest
    """
    os.system(f"grep -w {tf} {arch_SKCM_path} > {fold_txt}/{tf_txt}_tmp.bed")
    os.system(f"bedtools sort -i {fold_txt}/{tf_txt}_tmp.bed > {fold_txt}/{tf_txt}.bed")
    os.remove(f"{fold_txt}/{tf_txt}_tmp.bed")

    def slop_it(l_sz, r_sz, fl):
        """Create temporary .bed file with modified for coordinates using 'slop' for correct fasta retrieval"""
        os.system(f"bedtools slop -i {fold_txt}/{tf_txt}{fl}.bed -g {gen_fa_fai} -l {l_sz} -r {r_sz} > {fold_txt}/{tf_txt}{fl}_{l_sz}_{r_sz}_tmp_slop.bed")
        """Retrieve fasta sequences for coordinates in modified .bed file"""
        os.system(f"bedtools getfasta -fi {gen_fa} -bed {fold_txt}/{tf_txt}{fl}_{l_sz}_{r_sz}_tmp_slop.bed -fo {fold_txt}/{tf_txt}{fl}_{l_sz}_{r_sz}_tmp_slop_seq -s -bedOut")

        tf_bed = BedTool(f'{fold_txt}/{tf_txt}{fl}_{l_sz}_{r_sz}_tmp_slop.bed')
        bed_df = pybedtools.BedTool.to_dataframe(tf_bed)
        seqs = pd.read_table(f"{fold_txt}/{tf_txt}{fl}_{l_sz}_{r_sz}_tmp_slop_seq", header=None)
        tf_bed_df = bed_df[['chrom', 'start', 'end', 'strand']].copy()
        tf_bed_df['seq'] = seqs[0].str.upper()
        """Add reverse complement of fasta sequence to bedfile dataframe"""
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        tf_bed_df['seq_rv'] = ""
        for i in progressbar.progressbar(range(tf_bed_df.shape[0])):
            tf_bed_df['seq_rv'][i] = "".join(complement.get(base, base) for base in reversed(tf_bed_df['seq'][i]))

        os.remove(f"{fold_txt}/{tf_txt}{fl}_{l_sz}_{r_sz}_tmp_slop.bed")
        os.remove(f"{fold_txt}/{tf_txt}{fl}_{l_sz}_{r_sz}_tmp_slop_seq")
        tf_bed_df_reordered = tf_bed_df[['chrom','start','end','seq','seq_rv','strand']].copy()
        bed_final = pybedtools.BedTool.from_dataframe(tf_bed_df_reordered)

        bed_final.saveas(f"{fold_txt}/{tf_txt}_{l_sz}_seq.bed")


    #slop_it(30, 30, '')
    #slop_it(32, 32, '')
    slop_it(shft_sz, shft_sz, '')
    slop_it(shft_sz+2, shft_sz+2, '')

def makeLogo(fold_txt, tf_txt, w_len, gen):
    """
    Generate TFBSA PWM logo for later plots
    String:param fold_txt: Folder to save file to
    String:param tf_txt: Modified TF name for string manipulation
    Int:param w_len: Length of TFBS
    String:param gen: Path to reference genome files
    Saved Output:
    1. TFBS_PWM_logo.png
    """
    f_txt = f'{fold_txt}/{tf_txt}'
    gen_fa = f'{gen}/hg19.fa'
    os.system(f"bedtools getfasta -fi {gen_fa} -bed {f_txt}.bed -fo {f_txt}_tmp_seqs -s")
    filt_seqs = []
    for record in SeqIO.parse(f'{f_txt}_tmp_seqs', 'fasta'):
        if len(record.seq) == w_len:
            filt_seqs.append(record)
    SeqIO.write(filt_seqs, f'{f_txt}_tmp_filt_seqs', 'fasta')
    os.remove(f'{f_txt}_tmp_seqs')
    os.system(f"weblogo -f {f_txt}_tmp_filt_seqs -o {f_txt}_logo.png -F png -i 0 -s large -t {tf_txt} --errorbars NO --color-scheme classic --resolution 600") #OUTPUT1
    os.remove(f'{f_txt}_tmp_filt_seqs')

motifs = pd.read_csv(f'TFs/{args.transcription_factors}', usecols=['Name','Width','Name_str'])
motifs_r = motifs.to_records(index=False)
#g_path, f_path = '../genomes/hg19', 'TFs'
#archetype_site_Skin_intersect
g_path, f_path, c_path = args.path_to_genome_files, 'TFs', args.path_to_tf_clusters
for n, w, ns in motifs_r:
    print(f'Processing: {n}')
    main(f_path, n, ns, c_path, g_path, 200)
    makeLogo(f_path, ns, w, g_path)
