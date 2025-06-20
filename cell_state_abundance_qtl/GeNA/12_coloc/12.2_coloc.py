import pandas as pd
import click 
import os
import re
import sys
import logging
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

# 1) get lead csaQTL GWAS loci
# 2) Extract the SNP GWAS data for the cis-window (lead csaQTL GWAS loci +/- window; default: 100kB)
# 3) Run coloc for each lead csaQTL (if the SNP GWAS data has at least one significant variant; default pval <5e-8)
# 4) Write the results to a TSV file

def coloc_runner(
    gwas, gwas_pheno_name, csaqtl_file_path, celltype, outdir,  save = True, common_maf_threshold=0.05
):
    """ Modified function from: https://github.com/populationgenomics/saige-tenk10k/blob/main/coloc/coloc_runner.py
    """

    ro.r('suppressMessages(library(coloc))')
    ro.r('suppressMessages(library(tidyverse))')

    with (ro.default_converter + pandas2ri.converter).context():
        gwas_r = ro.conversion.get_conversion().py2rpy(gwas)
    ro.globalenv['gwas_r'] = gwas_r
    ro.r(
        '''
    gwas_r = gwas_r %>% dplyr::select(beta, varbeta, position,snp)
    gwas_r = gwas_r %>% filter((beta!=0) | (varbeta!=0))
    gwas_r = gwas_r %>% distinct(snp, .keep_all = TRUE)
    gwas_r = gwas_r%>% as.list()
    gwas_r$type = 'cc'

    ''',
    )
    csaqtl = pd.read_csv(
        csaqtl_file_path,
        sep='\t',
    )
    csaqtl['chr'] = csaqtl['#CHROM'].apply(lambda chr: f'chr{chr}')
    csaqtl['beta'] = csaqtl['BETA']
    csaqtl['se'] = csaqtl['SE']
    csaqtl['position'] = csaqtl['POS']
    csaqtl['snp'] = csaqtl['ID'].apply(lambda snp: 'chr' + snp.replace(':', '_'))
    # while I figure out if it's easy to extract sdY, give N and MAF instead
    # https://chr1swallace.github.io/coloc/articles/a02_data.html#what-if-i-dont-have-sdy
    csaqtl['MAF'] = csaqtl['A1_FREQ']
    csaqtl['N'] = csaqtl['OBS_CT'] # n_samples 
    # subset to only results for variants above chosen MAF (0.05 by default)
    csaqtl = csaqtl[csaqtl['MAF'] > common_maf_threshold] 
    with (ro.default_converter + pandas2ri.converter).context():
        csaqtl_r = ro.conversion.get_conversion().py2rpy(csaqtl)

    lead_snp = "_".join(csaqtl_file_path.split('/')[-1].split("_")[2:4]) # Parse lead SNP from the 

    ro.globalenv['csaqtl_r'] = csaqtl_r
    ro.globalenv['lead_snp'] = lead_snp
    ro.globalenv['gwas_pheno_name'] = gwas_pheno_name
    ro.r(
        '''
    csaqtl_r = csaqtl_r %>% filter(!is.na(beta))
    csaqtl_r = csaqtl_r %>% distinct(snp, .keep_all = TRUE)
    csaqtl_r$varbeta = csaqtl_r$se**2
    csaqtl_r$position = csaqtl_r$pos
    csaqtl_r = csaqtl_r %>% select(beta, varbeta, position, snp, N, MAF)

    csaqtl_r = csaqtl_r %>% as.list()
    csaqtl_r$type = 'quant'

    my.res <- coloc.abf(dataset1=gwas_r,
                    dataset2=csaqtl_r)
    p_df <- data.frame(lead_snp,gwas_pheno_name, my.res$summary[1], my.res$summary[2], my.res$summary[3], my.res$summary[4], my.res$summary[5], my.res$summary[6])
    names(p_df) <- c('lead_snp', 'gwas_pheno_name', 'nsnps_coloc_tested','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')
    ''',
    )

    # convert to pandas df
    with (ro.default_converter + pandas2ri.converter).context():
        pd_p4_df = ro.conversion.get_conversion().rpy2py(ro.r('p_df'))

    # add cell type and chrom annotation to df
    pd_p4_df['celltype'] = celltype
    pd_p4_df['chrom'] = csaqtl['chr'].iloc[0]

    if save: 
        coloc_results_file=f'{outdir}/{gwas_pheno_name}_{lead_snp}.tsv'
    # write to GCS
        pd_p4_df.to_csv(
            coloc_results_file,
            sep='\t',
            index=False,
        )
    else:
        return pd_p4_df

@click.option(
    '--gwas-file-path',
    help='Path to the SNP GWAS file',
)
@click.option(
    '--gwas-pheno-name',
    help='name of the gwas phenotype being colocalised',
)
@click.option(
    '--csaqtl-file-path',
    help='Path to the csaQTL fixed-phenotype summary statistics files',
)
@click.option(
    '--celltype',
    help='Name of the cell type'
)
@click.option(
    '--common-maf-threshold',
        help='Common variant MAF threshold - currently not used ',
        default=0.01)
@click.option(
    '--cis-window',
        help='Window around lead csaQTL variant to test for colocalisations',
        default=100000)
@click.option(
    '-gwas-significance-threshold',
    help='Significance threshold for GWAS. Only GWAS with > 1 variant in the cis region will be tested for colocalisation',
    default=5e-8
)
@click.option(
    '--outdir',
    help='Output directory')
@click.command()
def main(
    gwas_file_path: str,
    gwas_pheno_name: str,
    csaqtl_file_path: str,
    celltype: str,
    common_maf_threshold: float,
    cis_window: int,
    gwas_significance_threshold: float,
    outdir: str,

):   
    os.makedirs(outdir, exist_ok = True)
    chr, pos = [int(i) for i in csaqtl_file_path.split('/')[-1].split("_")[2:4]]

    pattern = '_snp_str_gwas_results_hg38_chr[0-9]+'

    # skip if for different chromosome 
    if re.search(pattern, gwas_file_path):
        gwas_chr = int(re.search(pattern, gwas_file_path).group().replace('_snp_str_gwas_results_hg38_chr', ''))
        if gwas_chr != chr:
            logging.warning(f'GWAS file for different chr than target SNP, skipping....')
            return

    # read in GWAS data
    gwas = pd.read_csv(gwas_file_path, sep = "\t")
    
    # Filter GWAS to include only SNPs within cis region around lead-snp 
    gwas = gwas.loc[(gwas['chromosome'] == f'chr{chr}') & (gwas['position'] > pos - cis_window) & (gwas['position'] < pos + cis_window),:]

    # Do some checks on the gwas data 
    if sum((gwas['chromosome'] == f'chr{chr}')) == 0:
        logging.warning(f'No SNPs for chr{chr} in gwas file: {os.path.basename(gwas_file_path)} skipping....')
    elif gwas['p_value'].min() > gwas_significance_threshold:
        logging.warning(f'No significant SNP GWAS data in the cis-window around chr{chr} {pos}: skipping....')
    elif gwas.empty: 
        logging.warning(f'No SNPs present in cis window {cis_window}B around chr{chr} {pos}: skipping....')
    else:
        logging.info(f'Running coloc for lead snp at chr{chr} {pos} and gwas: {gwas_file_path}')

        # run coloc 
        coloc_runner(gwas=gwas, gwas_pheno_name=gwas_pheno_name, 
                    csaqtl_file_path=csaqtl_file_path, celltype=celltype, outdir=outdir)


if __name__ == '__main__':
    main()


    ## TESTING CODE 
    # # # test run coloc for one sample 
    # gwas_file_path = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/gwas/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_eosinophil_count_snp_str_gwas_results_hg38_chr15.tab.gz'
    # gwas_pheno_name = 'lymphocyte_count'
    # csaqtl_file_path = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/plink_fixed_pheno/Monocyte/gwas_local_15_79972101_cis.spheno_15_79972101.glm.linear'
    # celltype = 'B'
    # outdir = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/B/'
    # cis_window = 10000
    # gwas_significance_threshold = 5e-8

    # os.makedirs(outdir, exist_ok = True)
    # chr, pos = [int(i) for i in csaqtl_file_path.split('/')[-1].split("_")[2:4]]

    # pattern = '_snp_str_gwas_results_hg38_chr[0-9]+'

    # # skip if for different chromosome 
    # if re.search(pattern, gwas_file_path):
    #     gwas_chr = int(re.search(pattern, gwas_file_path).group().replace('_snp_str_gwas_results_hg38_chr', ''))
    #     if gwas_chr != chr:
    #         logging.warning(f'GWAS file for different chr than target SNP, skipping....')
    #         # return
    
    # # read in GWAS data
    # gwas = pd.read_csv(gwas_file_path, sep = "\t")
    
    # # Filter GWAS to include only SNPs within cis region around lead-snp 
    # gwas = gwas.loc[(gwas['chromosome'] == f'chr{chr}') & (gwas['position'] > pos - cis_window) & (gwas['position'] < pos + cis_window),:]

    # # Do some checks on the gwas data 
    # if sum((gwas['chromosome'] == f'chr{chr}')) == 0:
    #     logging.warning(f'No SNPs for chr{chr} in gwas file: {os.path.basename(gwas_file_path)} skipping....')
    # elif gwas['p_value'].min() > gwas_significance_threshold:
    #     logging.warning(f'No significant SNP GWAS data in the cis-window around chr{chr} {pos}: skipping....')
    # elif gwas.empty: 
    #     logging.warning(f'No SNPs present in cis window {cis_window}B around chr{chr} {pos}: skipping....')
    # else:
    #     logging.info(f'Running coloc for lead snp at chr{chr} {pos} and gwas: {gwas_file_path}')
    #     # run coloc 
    #     x = coloc_runner(gwas=gwas, gwas_pheno_name=gwas_pheno_name, 
    #                 csaqtl_file_path=csaqtl_file_path, celltype=celltype, outdir=outdir, save = False)

    # x

 
    # gwas.columns = [f'gwas_{col}' for col in gwas.columns]
 
    # # csaqtl.columns = [f'csawtl_{col}' for col in csaqtl.columns]
    # merged_df = pd.merge(gwas, csaqtl, left_on="gwas_position", right_on = "csawtl_position")
    # merged_df.loc[merged_df['gwas_snp'] == merged_df['csawtl_snp'],:]