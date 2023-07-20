import pandas as pd
import progressbar
import os
import argparse

parser = argparse.ArgumentParser(prog='background4mers',
                                 description='Count 4mers at each position in region for background modeling')
parser.add_argument('data', nargs='+')
parser.add_argument('-tf', '--transcription_factors')
args = parser.parse_args()

def init4mers():
    '''
    Initialize an empty 4mer dictionary of dipyrimdines with all possible 1mer flanks
    Dictionary:return: Dipyrimidine 4mer dictionary
    '''
    yys = ('CC', 'CT', 'TC', 'TT')
    ns = ('A', 'C', 'T', 'G')
    qn_l = list()
    for qn in yys:
        for n1 in ns:
            for n2 in ns:
                qn_l.append((f'{n1}{qn}{n2}'))
    qn = dict.fromkeys(qn_l, 0)
    return qn

def main(bg, bf, d_file):
    """
    String:param bkgrnd: TF (or region) of interest
    String:param bkgrnd_txt: Flanks for region of interest
    String:param bf: Region of interest windows bed file path
    String:param d_file: CPD data file name
    Saved Output:
    1. Dictionary_of_total_count_of_dipyrimdine_4mers_in_background_region.pkl
    2. Summary_stats_for_CPD_signal_accross_4mers_in_preprocessed_CPDseq
    """

    def count4mers(fr):
        '''
        Total counts for all dipyrimdine 4mers occuring in the regions of interest windows
        Dictionary:param fr: Empty 4mer dictionary
        String:param bed_file: Regions of interest windows bed file path
        Dictionary:return: 4mer dictionary with counts
        '''
        for f in fr:
            fr[f] = 0
        yys = ['CC', 'CT', 'TC', 'TT']
        bkgrd_df = pd.read_table(f"TFs/{bf}.bed", header=None, usecols=[3,4], names=['seq', 'seq_rc'])
        bkgrd_df['seq'] = bkgrd_df['seq'].str.upper()
        bkgrd_df['seq_rc'] = bkgrd_df['seq_rc'].str.upper()
        sites = [bkgrd_df['seq'], bkgrd_df['seq_rc']]
        for s in sites:
            for seqs in progressbar.progressbar(s):
                seq_l = len(seqs)
                for i in range(seq_l - 3):
                    seq_n = seqs[i:i + 4]
                    if 'N' in seq_n:
                        continue
                    if seq_n[1:3] in yys:
                        fr[seq_n] = fr[seq_n] + 1
        return fr

    chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
              "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]


    def prep_cpds_chrom(all_cpds, chromo):
        chrom_cpds = all_cpds[all_cpds["chr"] == chromo]
        return chrom_cpds

    def count_bkgrnd_cpds(binding_sites_, cpds_, cpd_frs):
        cpds_df_length = cpds_.shape[0]
        i = 0
        for binding_site_start, binding_site_end, binding_site_seq in progressbar.progressbar(binding_sites_.iloc):
            assert binding_site_end > binding_site_start
            while True:
                if i > cpds_df_length-1:
                    break
                if cpds_[i]['start'] > binding_site_start:
                    if i > 0:
                        i = i - 1
                        continue
                    else:
                        break
                else:
                    break
            while i < cpds_df_length:
                cpd_start, cpd_count, cpd_seq = cpds_[i]
                if 'N' in cpd_seq:
                    i = i + 1
                    continue
                if cpd_start <= binding_site_start:
                    i = i + 1
                    continue
                if cpd_start >= binding_site_end-1:
                    break
                elif cpd_start>binding_site_start and cpd_start<binding_site_end-1:
                    cpd_frs[cpd_seq].append(cpd_count)
                i = i + 1
        

    def org_CPD_4mers(cp, fr):
        '''
        Create summary statistics for the 4mer distribution of pre-processed CPD-seq datasubsets for region of interest
        Dataframe:param d: 4mer counts from CPD signal
        Dataframe:param q: Total 4mer counts in TFBS windows for region of interest
        Dataframe:return: Summary statistics of CPD signal for each 4mer
        '''
        tot = list()
        for f in fr:
            f_tot = fr[f]
            c_tot = sum(cp[f])
            c_n = len(cp[f])
            if c_n == 0:
                c_n = 1
                print(f'too few cpds for {bg}, frmr: {f}')
            m = c_tot / f_tot
            v = sum([(d - m) ** 2 for d in cp[f]]) / c_n
            tot.append((f, c_tot, f_tot, m, v))
        tot_df = pd.DataFrame(tot, columns=['seq', 'cpd_sum', 'frmr_sum', 'mean', 'var'])
        return tot_df

    print(f'Prepping 4mers: {bg}_{d_file}')
    frmrs_list = init4mers()
    frmrs = count4mers(frmrs_list)
    """Convert regions of interest windows 4mer dictionary to Dataframe"""
    #frmrs_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in frmrs.items()])).melt(var_name='seq', value_name='tot')
    p_cpds = pd.read_table(f"{fold_txt}/{bg}_intersect_{d_file}_dipy_proc_plus.bed", header=None, usecols=[0, 1, 3, 4], names=['chr', 'start', 'counts', 'seq'])
    m_cpds = pd.read_table(f"{fold_txt}/{bg}_intersect_{d_file}_dipy_proc_minus.bed", header=None, usecols=[0, 1, 3, 4], names=['chr', 'start', 'counts', 'seq'])
    tf_binding_sites_ = pd.read_table(f'TFs/{bf}.bed', header=None, usecols=[0, 1, 2, 3], names=['chr', 'start', 'end', 'seq'])
    tf_binding_sites_[["start", "end"]] = tf_binding_sites_[["start", "end"]].to_numpy(dtype='int')
    cpd_frmrs = init4mers()
    #var_frmrs = init4mers()
    for f in cpd_frmrs:
        cpd_frmrs[f] = []
        #var_frmrs[f] = []

    print('Counting CPDs...')
    for chromosome in chroms:
        """Iterate through all chromsomes"""
        plus_strand_cpds = prep_cpds_chrom(p_cpds, chromosome)
        minus_strand_cpds = prep_cpds_chrom(m_cpds, chromosome)
        binding_sites = tf_binding_sites_.loc[(tf_binding_sites_["chr"] == chromosome), ['start', 'end', 'seq']]
        if binding_sites.shape[0] == 0:
            continue
        count_bkgrnd_cpds(binding_sites, plus_strand_cpds[["start", "counts", "seq"]].to_records(index=False), cpd_frmrs)
        count_bkgrnd_cpds(binding_sites, minus_strand_cpds[["start", "counts", "seq"]].to_records(index=False), cpd_frmrs)

    """Convert CPD-seq dataset 4mer dictionary to Dataframe"""
    #cpd_frmrs_all_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in cpd_frmrs_all.items()]))
    cpd_frmrs_all_df_final = org_CPD_4mers(cpd_frmrs, frmrs)
    cpd_frmrs_all_df_final.to_csv(f'{fold_txt}/{bg}_intersect_{d_file}_cpd_4mers_BKGRND.csv', index=False) #OUTPUT2


fold_txt = 'CPDseq_data/background'                     # Path output will be saved to
#d_files = ['WT_CSB_8J_0hr_CPD_1bp_sorted']  # Name of CPD-seq experiment XPC_12J_NakedDNA_S22_CPD_1bp_sorted
d_files = args.data

if not os.path.exists(fold_txt):
    os.system(f"mkdir {fold_txt}")
motifs = pd.read_csv(f'TFs/{args.transcription_factors}', usecols=['Width','Name_str'])
motifs_r = motifs.to_records(index=False)
for df in d_files:
    for w, n in motifs_r:
        main(f'{n}_flanks_unmerged_200to20', f'{n}_flanks_unmerged_200to20_seq', df)
        print("")
    print("")
