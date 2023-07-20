import pandas as pd
import os
import argparse
parser = argparse.ArgumentParser(prog='PrepMuts',
                                 description='Pre-processes mutation data')
parser.add_argument('-tf', '--transcription_factors')
parser.add_argument('-f', '--path_to_mutation_data')
args = parser.parse_args()

def main(old_tf, tf_txt, lng, sz):
    """
    Shift SKCM mutation data so that 0th index is at the start of TFBS
    String:param tf: Unmodified TF name
    String:param tf_txt: Modified TF name for string manipulation
    Int:param lng: Width of TFBS in bp
    String:param fold_txt: Path to original TF_skcm_extended_mutations_counts.tsv file
    Int:param sz: Desired TFBS window size to save mutation data for
    Saved Output:
    1. TF_skcm_mutation_counts.csv
    """
    if not os.path.exists('mutations'):
        os.system(f"mkdir mutations")
    """Input mutation file (from Harshit Sahay) contains cummulative SKCM mutation counts at each position in 1kb windows around
    Vierstra TFBS"""
    try:
        tf_d = pd.read_table(f"{fold_txt}/hg19.archetype_site_{old_tf}_intersection_SKCM.bed_extended_500_sequence_mutation_counts.tsv",
                             usecols=[0, 3])
    except:
        print(f'{old_tf}_skcm_extended_mutations_counts.tsv not found')
        return
    tf_d = tf_d.rename(columns={"mutation_position_from_center": "Position_old", "Total_counts": "Total_Mutations"})
    shft = lng//2
    tf_d['Position'] = tf_d['Position_old'] + shft
    tf_d_c = tf_d[["Position", "Total_Mutations"]].copy()
    half = tf_d_c.shape[0]//2
    """Only want mutation counts in TFBS window defined by sz to match future CPD analysis"""
    tf_d_c.sort_values(by=['Position'])[(half-shft-sz):(half+shft+sz)].to_csv(f'mutations/{tf_txt}_skcm_mutation_counts_{sz}.csv', index=False, sep='\t') #OUTPUT1


motifs = pd.read_csv(f'TFs/{args.transcription_factors}', usecols=['Width', 'Name_str'])
motifs_r = motifs.to_records(index=False)
#fold_txt = '/usr/project/xtmp/hs239/mutation_counts_jun_23/skcm'
fold_txt = args.path_to_mutation_data
for w, ns in motifs_r:
    print(f'Processing: {ns}')
    old_n = ns.replace('_', '.')
    main(old_n, ns, w, 200)
