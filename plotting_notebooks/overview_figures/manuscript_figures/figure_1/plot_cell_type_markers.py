import scanpy as sc 
import pandas as pd
import matplotlib

sc.settings.figdir = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/figure_1/figures/'

adata = sc.read("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_concatenated_harmony_filtered_min1000genes.h5ad")

# Plot some established pbmc marker genes 
# from https://www.sc-best-practices.org/cellular_structure/annotation.html#:~:text=marker_genes%20%3D%20%7B,negative%20marker%0A%7D

# set plot order
adata.obs.columns

adata.obs["celltype"] = adata.obs["wg2_scpred_prediction"].str.replace("_", " ")

marker_genes = {
    # CD4 T cells
    "CD4 TCM": ["IL7R", "TMSB10", "CD4", "ITGB1", "LTB", "TRAC", "AQP3", "LDHB", "IL32", "MAL"],
    "CD4 Naive": ["TCF7", "CD4", "CCR7", "IL7R", "FHIT", "LEF1", "MAL", "NOSIP", "LDHB", "PIK3IP1"],
    "CD4 TEM": ["IL7R", "CCL5", "FYB1", "GZMK", "IL32", "GZMA", "KLRB1", "TRAC", "LTB", "AQP3"],
    "CD4 CTL": ["GZMH", "CD4", "FGFBP2", "ITGB1", "GZMA", "CST7", "GNLY", "B2M", "IL32", "NKG7"],
    "Treg": ["RTKN2", "FOXP3",  "CD4", "IL2RA", "TIGIT", "CTLA4", "FCRL3", "LAIR2", "IKZF2"], # "AC133644.2",
    "CD4 Proliferating": ["MKI67", "TOP2A", "PCLAF", "CENPF", "TYMS", "NUSAP1", "ASPM", "PTTG1", "TPX2", "RRM2"],

    # unconventional T cells 
    "gdT": ["TRDC", "TRGC1", "TRGC2", "KLRC1", "NKG7", "TRDV2", "CD7", "TRGV9", "KLRD1", "KLRG1"],
    "MAIT": ["KLRB1", "NKG7", "GZMK", "IL7R", "SLC4A10", "GZMA", "CXCR6", "PRSS35", "RBM24", "NCR3"],
    "dnT": ["PTPN3", "MIR4422HG", "NUCB2", "CAV1", "DTHD1", "GZMA", "MYB", "FXYD2", "GZMK"], # "AC004585.1"
    "ILC": ["KIT", "TRDC", "TTLL10", "LINC01229", "SOX4", "KLRB1", "TNFRSF18", "TNFRSF4", "IL1R1", "HPGDS"],

    # CD8 T cells
    "CD8 TEM": ["CCL5", "GZMH", "CD8A", "TRAC", "KLRD1", "NKG7", "GZMK", "CST7", "CD8B", "TRGC2"],
    "CD8 Naive": ["CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"],
    "CD8 TCM": ["CD8B", "ANXA1", "CD8A", "KRT1", "LINC02446", "YBX3", "IL7R", "TRAC", "NELL2", "LDHB"],
    "CD8 Proliferating": ["MKI67", "CD8B", "TYMS", "TRAC", "PCLAF", "CD3D", "CLSPN", "CD3G", "TK1", "RRM2"],

    # NK
    "NK": ["GNLY", "TYROBP", "NKG7", "FCER1G", "GZMB", "TRDC", "PRF1", "FGFBP2", "SPON2", "KLRF1"],
    "NK CD56bright": ["XCL2", "FCER1G", "SPINK2", "TRDC", "KLRC1", "XCL1", "SPTSSB", "PPP1R9A", "NCAM1", "TNFRSF11A"],
    "NK Proliferating": ["MKI67", "KLRF1", "TYMS", "TRDC", "TOP2A", "FCER1G", "PCLAF", "CD247", "CLSPN", "ASPM"],

    # B  
    "B naive": ["IGHM", "IGHD", "CD79A", "IL4R", "MS4A1", "CXCR4", "BTG1", "TCL1A", "CD79B", "YBX3"],
    "B intermediate": ["MS4A1", "TNFRSF13B", "IGHM", "IGHD", "AIM2", "CD79A", "LINC01857", "RALGPS2", "BANK1", "CD79B"],
    "B memory": ["MS4A1", "COCH", "AIM2", "BANK1", "SSPN", "CD79A", "TEX9", "RALGPS2", "TNFRSF13C", "LINC01781"],
    "Plasmablast": ["IGHA2", "MZB1", "TNFRSF17", "DERL3", "TXNDC5", "TNFRSF13B", "POU2AF1", "CPNE5",  "NT5DC2"], # "HRASLS2",

    # Monocytes 
    "CD14 Mono": ["S100A9", "CTSS", "S100A8", "LYZ", "VCAN", "S100A12", "IL1B", "CD14", "G0S2", "FCN1"],
    "CD16 Mono": ["CDKN1C", "FCGR3A", "PTPRC", "LST1", "IER5", "MS4A7", "RHOC", "IFITM3", "AIF1", "HES4"],

    # Dendritic 
    "cDC2": ["FCER1A", "HLA-DQA1", "CLEC10A", "CD1C", "ENHO", "PLD4", "GSN", "SLC38A1", "NDRG2", "AFF3"],
    "pDC": ["ITM2C", "PLD4", "SERPINF1", "LILRA4",  "TPM2", "MZB1", "SPIB", "IRF4", "SMPD3"], # "IL3RA",
    "cDC1": ["CLEC9A", "DNASE1L3", "C1orf54", "IDO1", "CLNK", "CADM1", "FLT3", "ENPP1", "XCR1", "NDRG2"],
    "ASDC": ["PPP1R14A", "LILRA4", "AXL", "SCT", "SCN9A", "LGMN", "DNASE1L3", "CLEC4C", "GAS6"], # "IL3RA", 

    # Other 
    "HSPC": ["SPINK2", "PRSS57", "CYTL1", "EGFL7", "GATA2", "CD34", "SMIM24", "AVP", "MYB", "LAPTM4B"],

    # "Eryth": ["HBD", "HBM", "AHSP", "ALAS2", "CA1", "SLC4A1", "IFIT1B", "TRIM58", "SELENBP1", "TMCC2"],
    # "Platelet": ["PPBP", "PF4", "NRGN", "GNG11", "CAVIN2", "TUBB1", "CLU", "RGS18", "GP9"],# "HIST1H2AC", 
}

adata.obs["celltype"]  = pd.Categorical(
    values=adata.obs["celltype"], categories=marker_genes.keys(), ordered=True
)


sc.pl.dotplot(adata, marker_genes, groupby="celltype", standard_scale="var", dendrogram = False, gene_symbols = "gene_name", color_map=matplotlib.cm.get_cmap('RdYlBu_r'), save = "marker_dotplot.pdf")

# sc.pl.stacked_violin(adata, marker_genes, groupby="celltype", standard_scale="var", dendrogram = True, gene_symbols = "gene_name", save = "marker_stacked_violin.pdf")

# sc.pl.heatmap(adata, marker_genes, groupby="celltype", standard_scale="var", dendrogram = False, gene_symbols = "gene_name", color_map="magma", save = "marker_heatmap.pdf")

# find top marker genes for each cell type ??
# may need to down sample prior to finding DEGs
# sc.tl.rank_genes_groups(adata, groupby="wg2_scpred_prediction", method="wilcoxon")

marker_genes_subset = {
    # CD4 T cells
    "CD4 TCM": ["IL7R", "TMSB10", "CD4", "ITGB1", "LTB", "TRAC", "AQP3", "LDHB", "IL32", "MAL"],
    "CD4 Naive": [ "CCR7", "FHIT", "LEF1", "MAL"],
    "CD4 TEM": ["IL7R", "CCL5", "FYB1", "GZMK", "IL32", "GZMA", "KLRB1", "TRAC", "LTB", "AQP3"],
    "CD4 CTL": ["GZMH", "CD4", "ITGB1", "B2M"],
    "Treg": ["RTKN2", "FOXP3", "IL2RA",  "IKZF2"], # "AC133644.2",
    "CD4 Proliferating": ["MKI67", "TOP2A", "PCLAF", "CENPF"],

    # unconventional T cells 
    "gdT": ["TRDC", "TRGC1", "TRGC2", "TRDV2", "TRGV9",  "KLRG1"],
    "MAIT": ["KLRB1", "GZMK", "SLC4A10", "NCR3"],
    "dnT": ["PTPN3", "MIR4422HG", "NUCB2", "DTHD1",  "MYB", "FXYD2"],
    "ILC": ["KIT", "TNFRSF18", "TNFRSF4", "IL1R1"],

    # CD8 T cells
    "CD8 TEM": ["CCL5", "GZMH", "CD8A", "TRAC", "KLRD1", "NKG7", "GZMK", "CST7", "CD8B", "TRGC2"],
    "CD8 Naive": ["CD8B", "S100B", "CCR7", "LINC02446"],
    "CD8 TCM": ["CD8B", "ANXA1", "CD8A", "KRT1", "LINC02446", "YBX3", "IL7R", "TRAC", "NELL2", "LDHB"],
    "CD8 Proliferating": ["CD8B", "CD3D"],

    # NK
    "NK": ["NKG7", "PRF1", "FGFBP2", "SPON2", "KLRF1"],
    "NK CD56bright": ["XCL2", "KLRC1", "XCL1", "SPTSSB", "PPP1R9A", "NCAM1"],
    "NK Proliferating": [],

    # B  
    "B naive": ["IGHM", "IGHD",  "IL4R", "TCL1A" ],
    "B intermediate": ["MS4A1","LINC01857", "RALGPS2", "BANK1"],
    "B memory": ["COCH", "AIM2", "SSPN", "TEX9",  "TNFRSF13C"],
    "Plasmablast": ["IGHA2", "MZB1", "TNFRSF17", "NT5DC2"],

    # Monocytes 
    "CD14 Mono": ["S100A9",  "S100A8", "VCAN", "S100A12", "CD14"],
    "CD16 Mono": ["CDKN1C", "FCGR3A", "MS4A7", "RHOC", "IFITM3", "AIF1", "HES4"],

    # Dendritic 
    "cDC2": ["FCER1A","CLEC10A", "CD1C", "ENHO"],
    "pDC": ["PLD4", "SERPINF1", "TPM2", "SMPD3"], # "IL3RA",
    "cDC1": ["CLEC9A", "IDO1", "CLNK", "CADM1", "FLT3"],
    "ASDC": ["PPP1R14A", "AXL", "SCT", "LGMN"], # "IL3RA", 

    # Other 
    "HSPC": ["SPINK2", "PRSS57", "CYTL1", "EGFL7", "GATA2", "CD34", "SMIM24"],

    # "Eryth": ["HBD", "HBM", "AHSP", "ALAS2", "CA1", "SLC4A1", "IFIT1B", "TRIM58", "SELENBP1", "TMCC2"],
    # "Platelet": ["PPBP", "PF4", "NRGN", "GNG11", "CAVIN2", "TUBB1", "CLU", "RGS18", "GP9"],# "HIST1H2AC", 
}

sc.pl.dotplot(adata, marker_genes_subset, groupby="celltype", standard_scale="var", dendrogram = False, gene_symbols = "gene_name", color_map=matplotlib.cm.get_cmap('RdYlBu_r'), save = "marker_dotplot_subset.pdf")
