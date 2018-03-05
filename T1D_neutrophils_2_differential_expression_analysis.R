#### scripts for analysis of RNA-seq data for Vecchio, Lo Buono, ..., Battaglia. 2018. Abnormal neutrophil signature in the blood and pancreas of pre-symptomatic and symptomatic type 1 diabetes.
### this file includes scripts for running differential expression analyses

##### set up environment: load packages #####

# load general packages
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
theme_set(
  theme_bw(20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, size=1)))
library(gplots)
library(ggthemes)

# load analysis_specific packages
library(edgeR)
library(limma)

# load packages with custom functions (available at https:://github.com/mjdufort)
library(countSubsetNorm)
library(RNAseQC)
library(limmaTools)
library(miscHelpers)
library(geneSetTools)

# load maggritr (load last to ensure default pipe function selection)
library(magrittr)


##### load data if not using results from T1D_neutrophils_1_process_count_data.R (e.g. if using data from GEO) #####
# skip this section if objects were generated using the scripts in T1D_neutrophils_1_load_data.R

# need to read in annotation from create counts.final, master, counts.final.inc_outliers, and master.inc_outliers
# read in annotation from GEO, and standardize and match column names
master.inc_outliers <-
  xlsx::read.xlsx("T1D_neutrophils_metadata_GEO.xlsx") %>%
  standardize_dimnames() %>%
  plyr::rename(
    replace=c(
      "title"="patient_id",
      "libraryid"="library_id",
      "group"="donor_type",
      "group_aab_status"="donor_type_aab_status"))

master.inc_outliers$library_id <-
  master.inc_outliers$library_id %>%
  str_extract("[0-9]+$")

# create version without outliers
outlier_libs.tmp <- c("3474", "4177", "4183")
master <-
  master.inc_outliers[
    !(master.inc_outliers$library_id %in% outlier_libs.tmp),]

## read in counts from GEO, and modify to match upstream processing
counts.final.inc_outliers <-
  read.table("raw_counts_T1D_neutrophils.txt", sep="\t") %>%
  mutate(HGNC.symbol=get_HGNC(rownames(.), type="protein_coding"))
colnames(counts.final.inc_outliers) <-
  colnames(counts.final.inc_outliers) %>%
  str_extract("[0-9]+$")

## sum counts.merged for duplicated HGNC symbols
# this also drops rows with HGNC.symbols==NA, which should include any genes that are not protein-coding
counts.final.inc_outliers <-
  aggregate(counts.final.inc_outliers[, -ncol(counts.final.inc_outliers)],
            by=list(counts.final.inc_outliers$HGNC.symbol), sum)
rownames(counts.final.inc_outliers) <- counts.final.inc_outliers$Group.1
counts.final.inc_outliers <-
  counts.final.inc_outliers[,-which(colnames(counts.final.inc_outliers) == "Group.1")]
colnames(counts.final.inc_outliers) <-
  colnames(counts.final.inc_outliers) %>%
  str_extract("[0-9]+$")

# create version of counts without outliers
counts.final <-
  counts.final.inc_outliers[
    , match(master$library_id, colnames(counts.final.inc_outliers))]


##### load data if using results of T1D_neutrophils_1_process_count_data.R #####
# skip this section if objects were NOT generated using the scripts in T1D_neutrophils_1_load_data.R

load(file="T1D_neutrophils_data_for_analysis.Rdata")

## generate "master" objects for simplicity of adapting downstream code
master <-
  sample_annotation.final %>%
  plyr::rename(replace=c("age_at_sample_coll"="age"))
master$donor_type <-
  factor(master$donor_type,
         levels=c("HC", "at_risk", "T1Dnew"))

master.inc_outliers <- sample_annotation.final.inc_outliers %>%
  plyr::rename(replace=c("age_at_sample_coll"="age"))
master.inc_outliers$donor_type <-
  factor(master.inc_outliers$donor_type,
         levels=c("HC", "at_risk", "T1Dnew"))


##### generate vwts objects for general downstream use #####

# for all samples (excluding outliers)
vwts.all <-
  counts.final %>%
  calc_norm_counts(
    design=master, libID_col="library_id",
    min_cpm = 1, min_libs_perc = 0.15,
    return_DGEcounts=TRUE) %>%
  voomWithQualityWeights(plot=TRUE)

# for all samples (including outliers)
vwts.all.inc_outliers <-
  counts.final.inc_outliers %>%
  calc_norm_counts(
    design=master.inc_outliers, libID_col="library_id",
    min_cpm = 1, min_libs_perc = 0.15,
    return_DGEcounts=TRUE) %>%
  voomWithQualityWeights(plot=TRUE)


##### create color palettes for use in heatmaps and such #####

pal.donor_type_aab_status <-
  ggthemes::colorblind_pal()(4)[c(1,3,4,2)] %>%
  setNames(levels(master$donor_type_aab_status))
pal.donor_type <-
  c(ggthemes::colorblind_pal()(4)[1],
    average_colors(ggthemes::colorblind_pal()(4)[3:4]),
    ggthemes::colorblind_pal()(4)[2]) %>%
  setNames(levels(master$donor_type))


##### Set up and run limma for donor_type_aab_status.inc_outliers #####

condition.tmp <-
  with(master.inc_outliers, !is.na(donor_type_aab_status))
master.donor_type_aab_status.inc_outliers <-
  master.inc_outliers %>%
  filter(condition.tmp) %>%
  fix_factors()
DGECounts.donor_type_aab_status.inc_outliers.tmp <-
  calc_norm_counts(
    counts=counts.final.inc_outliers,
    design=master.donor_type_aab_status.inc_outliers,
    libID_col="library_id",
    min_cpm = 1, min_libs_perc = 0.15,
    return_DGEcounts=TRUE,
    group=master.donor_type_aab_status.inc_outliers$participant_id)
master.donor_type_aab_status.inc_outliers <-
  master.donor_type_aab_status.inc_outliers[
    match(colnames(DGECounts.donor_type_aab_status.inc_outliers.tmp),
          master.donor_type_aab_status.inc_outliers[,"library_id"]),]

DesignMat.donor_type_aab_status.inc_outliers <-
  model.matrix(~ donor_type_aab_status + sex,
               data=master.donor_type_aab_status.inc_outliers)
vwts.donor_type_aab_status.inc_outliers <-
  voomWithQualityWeights(
    DGECounts.donor_type_aab_status.inc_outliers.tmp,
    design=DesignMat.donor_type_aab_status.inc_outliers, plot=TRUE, span=0.1)
vfit.donor_type_aab_status.inc_outliers <-
  lmFit(vwts.donor_type_aab_status.inc_outliers) %>%
  eBayes()
topGenes.donor_type_aab_status.inc_outliers <-
  vfit.donor_type_aab_status.inc_outliers %>%
  topTable(
    coef = which(str_detect(colnames(.), "donor_type_aab_status")),
    number=Inf, sort.by="F")

# isolate genes that have an FDR < .1
p_cut.tmp <- 0.1
topGenes.donor_type_aab_status.inc_outliers$threshold <-
  topGenes.donor_type_aab_status.inc_outliers$adj.P.Val < p_cut.tmp


## Heatmaps of counts by donor_type_aab_status

# use genes with FDR < 0.1
counts.tmp <-
  get_counts_sig_genes(
    counts=vwts.donor_type_aab_status.inc_outliers,
    topGenes=topGenes.donor_type_aab_status.inc_outliers,
    p_cut=0.1, fc_cut=NULL)

plot_gene_heatmap(
  counts.tmp,
  design=master.donor_type_aab_status.inc_outliers, libID_col="library_id",
  order_by_var="donor_type_aab_status",
  color_by_var="donor_type_aab_status",
  my_var_colors=
    pal.donor_type_aab_status[
      levels(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)],
  add_legend=TRUE,
  ColSideLabs=NA,
  labCol=
    str_replace(master.donor_type_aab_status.inc_outliers$batlab_id[
      order(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)], "-[0-9]F$", ""))

plot_gene_heatmap(
  counts.tmp, 
  design=master.donor_type_aab_status.inc_outliers, libID_col="library_id",
  order_by_var="aab_sample_ordinal",
  color_by_var="donor_type_aab_status",
  my_var_colors=
    pal.donor_type_aab_status[
      levels(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)],
  add_legend=TRUE,
  ColSideLabs=NA,
  labCol=
    str_replace(master.donor_type_aab_status.inc_outliers$batlab_id[
      order(master.donor_type_aab_status.inc_outliers$aab_sample_ordinal)], "-[0-9]F$", ""))

plot_gene_heatmap(
  counts.tmp, design=master.donor_type_aab_status.inc_outliers, libID_col="library_id",
  color_by_var="donor_type_aab_status",
  my_var_colors=
    pal.donor_type_aab_status[
      levels(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)],
  col_dendro=TRUE,
  add_legend=TRUE,
  ColSideLabs=NA,
  labCol=
    str_replace(master.donor_type_aab_status.inc_outliers$batlab_id, "-[0-9]F$", ""))


# Write list of significant + FC genes
write_sig_genes(
  topGenes.donor_type_aab_status.inc_outliers,
  file_prefix="donor_type_aab_status.inc_outliers.P0.1",
  method="combined",
  threshold_col="threshold")

# output list of all genes in ensembl format for DAVID enrichment with background list
write_sig_genes(
  topGenes.donor_type_aab_status.inc_outliers,
  file_prefix="donor_type_aab_status.inc_outliers.ensgene",
  method=c("ranked_list"),
  input_type="symbol", output_type="ensgene")


##### Run limma for donor_type_aab_status.inc_outliers.age #####

DesignMat.donor_type_aab_status.inc_outliers.age <-
  model.matrix(~ donor_type_aab_status + age + sex,
               data=master.donor_type_aab_status.inc_outliers)
vwts.donor_type_aab_status.inc_outliers.age <-
  voomWithQualityWeights(
    DGECounts.donor_type_aab_status.inc_outliers.tmp,
    design=DesignMat.donor_type_aab_status.inc_outliers.age, plot=TRUE, span=0.1)
vfit.donor_type_aab_status.inc_outliers.age <-
  lmFit(vwts.donor_type_aab_status.inc_outliers.age) %>%
  eBayes()
topGenes.donor_type_aab_status.inc_outliers.age <-
  vfit.donor_type_aab_status.inc_outliers.age %>%
  topTable(
    coef = which(str_detect(colnames(.), "donor_type_aab_status")),
    number=Inf, sort.by="F")

# isolate genes that have an FDR < .1
p_cut.tmp <- 0.1
topGenes.donor_type_aab_status.inc_outliers.age$threshold <-
  topGenes.donor_type_aab_status.inc_outliers.age$adj.P.Val < p_cut.tmp

# counts of sig genes and FC cut
vwts_sig.donor_type_aab_status.inc_outliers.age <-
  get_counts_sig_genes(
    vwts.donor_type_aab_status.inc_outliers.age,
    topGenes=topGenes.donor_type_aab_status.inc_outliers.age,
    threshold_col="threshold")

# what are the most significant and highest-abs(FC) genes?
head(topGenes.donor_type_aab_status.inc_outliers.age[
  order(topGenes.donor_type_aab_status.inc_outliers.age$P.Value),], 20)
head(topGenes.donor_type_aab_status.inc_outliers.age[
  order(abs(topGenes.donor_type_aab_status.inc_outliers.age$logFC), decreasing=TRUE),], 20)

## Heatmaps of counts by donor_type_aab_status

# use genes with FDR < 0.1
counts.tmp <-
  vwts.donor_type_aab_status.inc_outliers.age$E[
    match(
      rownames(topGenes.donor_type_aab_status.inc_outliers.age[
        topGenes.donor_type_aab_status.inc_outliers.age$adj.P.Val < 0.1,]),
      rownames(vwts.donor_type_aab_status.inc_outliers.age$E)), ]

plot_gene_heatmap(
  counts.tmp, design=master.donor_type_aab_status.inc_outliers, libID_col="library_id",
  order_by_var="donor_type_aab_status",
  color_by_var="donor_type_aab_status",
  my_var_colors=
    pal.donor_type_aab_status[
      levels(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)],
  add_legend=TRUE,
  ColSideLabs=NA,
  labCol=
    str_replace(master.donor_type_aab_status.inc_outliers$batlab_id[
      order(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)], "-[0-9]F$", ""))

plot_gene_heatmap(
  counts.tmp, design=master.donor_type_aab_status.inc_outliers, libID_col="library_id",
  color_by_var="donor_type_aab_status",
  my_var_colors=
    pal.donor_type_aab_status[
      levels(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)],
  col_dendro=TRUE,
  add_legend=TRUE,
  ColSideLabs=NA,
  labCol=
    str_replace(master.donor_type_aab_status.inc_outliers$batlab_id, "-[0-9]F$", ""))

# Write list of significant + FC genes
write_sig_genes(
  topGenes.donor_type_aab_status.inc_outliers.age,
  file_prefix="donor_type_aab_status.inc_outliers.age.P0.1",
  method="combined",
  threshold_col="threshold")


##### Run limma for donor_type_aab_status.inc_outliers.pur_eos_perc #####

DesignMat.donor_type_aab_status.inc_outliers.pur_eos_perc <-
  model.matrix(~ donor_type_aab_status + eosinophil_purity_percent + sex,
               data=master.donor_type_aab_status.inc_outliers)
vwts.donor_type_aab_status.inc_outliers.pur_eos_perc <-
  voomWithQualityWeights(
    DGECounts.donor_type_aab_status.inc_outliers.tmp,
    design=DesignMat.donor_type_aab_status.inc_outliers.pur_eos_perc, plot=TRUE, span=0.1)
vfit.donor_type_aab_status.inc_outliers.pur_eos_perc <-
  lmFit(vwts.donor_type_aab_status.inc_outliers.pur_eos_perc) %>%
  eBayes()
topGenes.donor_type_aab_status.inc_outliers.pur_eos_perc <-
  vfit.donor_type_aab_status.inc_outliers.pur_eos_perc %>%
  topTable(
    coef = which(str_detect(colnames(.), "donor_type_aab_status")),
    number=Inf, sort.by="F")

# isolate genes that have an FDR < 0.1
p_cut.tmp <- 0.1
topGenes.donor_type_aab_status.inc_outliers.pur_eos_perc$threshold <-
  topGenes.donor_type_aab_status.inc_outliers.pur_eos_perc$adj.P.Val < p_cut.tmp


## Heatmaps of counts by donor_type_aab_status

# use genes with FDR < 0.1
counts.tmp <-
  vwts.donor_type_aab_status.inc_outliers.pur_eos_perc$E[
    match(
      rownames(topGenes.donor_type_aab_status.inc_outliers.pur_eos_perc[
        topGenes.donor_type_aab_status.inc_outliers.pur_eos_perc$adj.P.Val < 0.1,]),
      rownames(vwts.donor_type_aab_status.inc_outliers.pur_eos_perc$E)), ]

plot_gene_heatmap(
  counts.tmp, design=master.donor_type_aab_status.inc_outliers, libID_col="library_id",
  order_by_var="donor_type_aab_status",
  color_by_var="donor_type_aab_status",
  my_var_colors=
    pal.donor_type_aab_status[
      levels(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)],
  add_legend=TRUE,
  ColSideLabs=NA,
  labCol=
    str_replace(master.donor_type_aab_status.inc_outliers$batlab_id[
      order(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)], "-[0-9]F$", ""))

plot_gene_heatmap(
  counts.tmp, design=master.donor_type_aab_status.inc_outliers, libID_col="library_id",
  color_by_var="donor_type_aab_status",
  my_var_colors=
    pal.donor_type_aab_status[
      levels(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)],
  col_dendro=TRUE,
  add_legend=TRUE,
  ColSideLabs=NA,
  labCol=
    str_replace(master.donor_type_aab_status.inc_outliers$batlab_id, "-[0-9]F$", ""))

# Write list of significant + FC genes
write_sig_genes(
  topGenes.donor_type_aab_status.inc_outliers.pur_eos_perc,
  file_prefix="donor_type_aab_status.inc_outliers.pur_eos_perc.P0.1",
  method="combined",
  threshold_col="threshold")

# output list of all genes in ensembl format for DAVID enrichment with background list
write_sig_genes(
  topGenes.donor_type_aab_status.inc_outliers.pur_eos_perc,
  file_prefix="donor_type_aab_status.inc_outliers.pur_eos_perc.ensgene",
  method=c("ranked_list"),
  input_type="symbol", output_type="ensgene")


##### Run limma for donor_type_aab_status.inc_outliers.age.pur_eos_perc #####

DesignMat.donor_type_aab_status.inc_outliers.age.pur_eos_perc <-
  model.matrix(~ donor_type_aab_status + age + eosinophil_purity_percent + sex,
               data=master.donor_type_aab_status.inc_outliers)
vwts.donor_type_aab_status.inc_outliers.age.pur_eos_perc <-
  voomWithQualityWeights(
    DGECounts.donor_type_aab_status.inc_outliers.tmp,
    design=DesignMat.donor_type_aab_status.inc_outliers.age.pur_eos_perc, plot=TRUE, span=0.1)
vfit.donor_type_aab_status.inc_outliers.age.pur_eos_perc <-
  lmFit(vwts.donor_type_aab_status.inc_outliers.age.pur_eos_perc) %>%
  eBayes()
topGenes.donor_type_aab_status.inc_outliers.age.pur_eos_perc <-
  vfit.donor_type_aab_status.inc_outliers.age.pur_eos_perc %>%
  topTable(
    coef = which(str_detect(colnames(.), "donor_type_aab_status")),
    number=Inf, sort.by="F")

# isolate genes that have an FDR < 0.1
p_cut.tmp <- 0.1
topGenes.donor_type_aab_status.inc_outliers.age.pur_eos_perc$threshold <-
  topGenes.donor_type_aab_status.inc_outliers.age.pur_eos_perc$adj.P.Val < p_cut.tmp


## Heatmaps of counts by donor_type_aab_status

# use genes with FDR < 0.1
counts.tmp <-
  vwts.donor_type_aab_status.inc_outliers.age.pur_eos_perc$E[
    match(
      rownames(topGenes.donor_type_aab_status.inc_outliers.age.pur_eos_perc[
        topGenes.donor_type_aab_status.inc_outliers.age.pur_eos_perc$adj.P.Val < 0.1,]),
      rownames(vwts.donor_type_aab_status.inc_outliers.age.pur_eos_perc$E)), ]

plot_gene_heatmap(
  counts.tmp, design=master.donor_type_aab_status.inc_outliers, libID_col="library_id",
  order_by_var="donor_type_aab_status",
  color_by_var="donor_type_aab_status",
  my_var_colors=
    pal.donor_type_aab_status[
      levels(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)],
  add_legend=TRUE,
  ColSideLabs=NA,
  labCol=
    str_replace(master.donor_type_aab_status.inc_outliers$batlab_id[
      order(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)], "-[0-9]F$", ""))


# use genes with FDR < 0.1, after removing age and eosinophil purity using removeBatchEffect
counts.tmp <-
  removeBatchEffect(
    vwts.donor_type_aab_status.inc_outliers.age.pur_eos_perc$E,
    covariates=
      master.donor_type_aab_status.inc_outliers[
        , c("age", "eosinophil_purity_percent")],
    design=
      model.matrix(
      ~ donor_type_aab_status + sex,
      data=master.donor_type_aab_status.inc_outliers))[
        match(
          rownames(topGenes.donor_type_aab_status.inc_outliers.age.pur_eos_perc[
            topGenes.donor_type_aab_status.inc_outliers.age.pur_eos_perc$adj.P.Val < 0.1,]),
          rownames(vwts.donor_type_aab_status.inc_outliers.age.pur_eos_perc$E)), ]

plot_gene_heatmap(
  counts.tmp,
  design=master.donor_type_aab_status.inc_outliers,
  libID_col="library_id",
  order_by_var="donor_type_aab_status",
  color_by_var="donor_type_aab_status",
  my_var_colors=
    pal.donor_type_aab_status[
      levels(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)],
  add_legend=TRUE,
  ColSideLabs=NA,
  labCol=
    str_replace(master.donor_type_aab_status.inc_outliers$batlab_id[
      order(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)], "-[0-9]F$", ""))

# Write list of significant + FC genes
write_sig_genes(
  topGenes.donor_type_aab_status.inc_outliers.age.pur_eos_perc,
  file_prefix="donor_type_aab_status.inc_outliers.age.pur_eos_perc.P0.1",
  method="combined",
  threshold_col="threshold")

rm_tmp(ask=FALSE)


##### Venn diagram of DE genes for each group vs HC in donor_type_aab_status.inc_outliers #####

## extract topGenes object for each comparison vs. HC
topGenes.donor_type_aab_status.inc_outliers.each_vs_HC <- list()
for (i in (levels(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)[-1])) {
  topGenes.donor_type_aab_status.inc_outliers.each_vs_HC[[i]] <-
    topTable(vfit.donor_type_aab_status.inc_outliers,
             coef = paste0("donor_type_aab_status", i), number=Inf, sort.by="P")
}

# all DE genes in each with FDR < 0.1
gplots::venn(
  lapply(
    topGenes.donor_type_aab_status.inc_outliers.each_vs_HC,
    function(x) {rownames(x)[x[,"adj.P.Val"] < 0.1]}))

rm_tmp(ask=FALSE)


##### Extract results for each pairwise comparison in donor_type_aab_status.inc_outliers #####

## fit all pairwise contrasts in one object; this specification of contrasts uses the intercept,
# which corresponds to the value for the HC group
# so the other coefficients are only the deviations from that intercept
contrasts.tmp <-
  cbind(c(0,1,0,0,0),
        c(0,0,1,0,0),
        c(0,0,0,1,0),
        c(0,-1,1,0,0),
        c(0,-1,0,1,0),
        c(0,0,-1,1,0))
colnames(contrasts.tmp) <-
  c("at_risk_Aab_neg_vs_HC",
    "at_risk_Aab_pos_vs_HC",
    "T1Dnew_vs_HC",
    "at_risk_Aab_pos_vs_at_risk_Aab_neg",
    "T1Dnew_vs_at_risk_Aab_neg",
    "T1Dnew_vs_at_risk_Aab_pos")
vfit.donor_type_aab_status.inc_outliers.pairwise <- 
  contrasts.fit(vfit.donor_type_aab_status.inc_outliers, contrasts=contrasts.tmp) %>%
  eBayes()

## already generated the topTable object for the full fit including all factor levels
topGenes.donor_type_aab_status.inc_outliers

## extract topGenes object for each pairwise comparison
topGenes.donor_type_aab_status.inc_outliers.pairwise <- list()
for (i in colnames(contrasts.tmp)) {
  topGenes.donor_type_aab_status.inc_outliers.pairwise[[i]] <-
    topTable(vfit.donor_type_aab_status.inc_outliers.pairwise, coef=i, number=Inf, sort.by="P")
}

## use only the significant genes from the F test
# and determine if they're significant in each comparison
vennData.donor_type_aab_status.inc_outliers.pairwise.sig_F_only <-
  data.frame(
    gene=
      rownames(topGenes.donor_type_aab_status.inc_outliers)[
        topGenes.donor_type_aab_status.inc_outliers$adj.P.Val < 0.1],
    at_risk_Aab_neg_vs_HC=NA,
    at_risk_Aab_pos_vs_HC=NA,
    T1Dnew_vs_HC=NA,
    at_risk_Aab_pos_vs_at_risk_Aab_neg=NA,
    T1Dnew_vs_at_risk_Aab_neg=NA,
    T1Dnew_vs_at_risk_Aab_pos=NA)
for (i in names(topGenes.donor_type_aab_status.inc_outliers.pairwise)) {
  vennData.donor_type_aab_status.inc_outliers.pairwise.sig_F_only[,i] <-
    topGenes.donor_type_aab_status.inc_outliers.pairwise[[i]][
      match(vennData.donor_type_aab_status.inc_outliers.pairwise.sig_F_only$gene,
            rownames(topGenes.donor_type_aab_status.inc_outliers.pairwise[[i]])),
      "adj.P.Val"] < 0.1}

# plot a venn diagram
vennDiagram(
  apply(
    as.matrix(
      vennData.donor_type_aab_status.inc_outliers.pairwise.sig_F_only[,-1]),
    2, as.numeric),
  names=
    str_replace_all(
      colnames(vennData.donor_type_aab_status.inc_outliers.pairwise.sig_F_only)[-1],
      "_vs_", " vs. "),
  cex=c(1.2,1.0,0.7))
# 0 genes between any of them


## use all the genes, and determine if they're significant in each comparison
vennData.donor_type_aab_status.inc_outliers.pairwise.all_genes <-
  data.frame(
    gene=
      rownames(topGenes.donor_type_aab_status.inc_outliers),
    at_risk_Aab_neg_vs_HC=NA,
    at_risk_Aab_pos_vs_HC=NA,
    T1Dnew_vs_HC=NA,
    at_risk_Aab_pos_vs_at_risk_Aab_neg=NA,
    T1Dnew_vs_at_risk_Aab_neg=NA,
    T1Dnew_vs_at_risk_Aab_pos=NA)
for (i in names(topGenes.donor_type_aab_status.inc_outliers.pairwise)) {
  vennData.donor_type_aab_status.inc_outliers.pairwise.all_genes[,i] <-
    topGenes.donor_type_aab_status.inc_outliers.pairwise[[i]][
      match(vennData.donor_type_aab_status.inc_outliers.pairwise.all_genes$gene,
            rownames(topGenes.donor_type_aab_status.inc_outliers.pairwise[[i]])),
      "adj.P.Val"] < 0.1}

# plot a venn diagram
vennDiagram(
  apply(
    as.matrix(
      vennData.donor_type_aab_status.inc_outliers.pairwise.all_genes[,-1]),
    2, as.numeric),
  names=
    str_replace_all(
      colnames(vennData.donor_type_aab_status.inc_outliers.pairwise.all_genes)[-1],
      "_vs_", " vs. "),
  cex=c(1.2,1.0,0.7))

# compile differences and output to a file
pairwise_sig_genes.donor_type_aab_status.inc_outliers <-
  data.frame(
    comparison=
      c(names(topGenes.donor_type_aab_status.inc_outliers.pairwise), "total_genes"),
    sig_genes_only=
      c(colSums(vennData.donor_type_aab_status.inc_outliers.pairwise.sig_F_only[,-1]),
        nrow(vennData.donor_type_aab_status.inc_outliers.pairwise.sig_F_only)),
    all_genes=
      c(colSums(vennData.donor_type_aab_status.inc_outliers.pairwise.all_genes[,-1]),
        nrow(vennData.donor_type_aab_status.inc_outliers.pairwise.all_genes)))
write.csv(
  pairwise_sig_genes.donor_type_aab_status.inc_outliers,
  file="pairwise_sig_genes.donor_type_aab_status.inc_outliers.csv",
  row.names=FALSE)

rm_tmp(ask=FALSE)


##### Set up and run limma for donor_type_aab_status #####

condition.tmp <- with(master, !is.na(donor_type_aab_status))
master.donor_type_aab_status <-
  master %>%
  filter(condition.tmp) %>%
  fix_factors()
DGECounts.donor_type_aab_status.tmp <-
  calc_norm_counts(
    counts=counts.final,
    design=master.donor_type_aab_status,
    libID_col="library_id",
    min_cpm = 1, min_libs_perc = 0.15,
    return_DGEcounts=TRUE,
    group=master.donor_type_aab_status$participant_id)
master.donor_type_aab_status <-
  master.donor_type_aab_status[
    match(colnames(DGECounts.donor_type_aab_status.tmp),
          master.donor_type_aab_status[,"library_id"]),]

DesignMat.donor_type_aab_status <-
  model.matrix(
    ~ donor_type_aab_status + sex,
    data=master.donor_type_aab_status)
vwts.donor_type_aab_status <-
  voomWithQualityWeights(
    DGECounts.donor_type_aab_status.tmp,
    design=DesignMat.donor_type_aab_status, plot=TRUE, span=0.1)
vfit.donor_type_aab_status <-
  lmFit(vwts.donor_type_aab_status) %>%
  eBayes()
topGenes.donor_type_aab_status <-
  vfit.donor_type_aab_status %>%
  topTable(
    coef = which(str_detect(colnames(.), "donor_type_aab_status")),
    number=Inf, sort.by="F")

# isolate genes that have an FDR < .1
p_cut.tmp <- 0.1
topGenes.donor_type_aab_status$threshold <-
  topGenes.donor_type_aab_status$adj.P.Val < p_cut.tmp

# counts of sig genes and FC cut
vwts_sig.donor_type_aab_status <-
  get_counts_sig_genes(
    vwts.donor_type_aab_status,
    topGenes=topGenes.donor_type_aab_status,
    threshold_col="threshold")

## Heatmaps of counts by donor_type_aab_status

## use genes with FDR < 0.1
counts.tmp <-
  vwts.donor_type_aab_status$E[
    match(rownames(topGenes.donor_type_aab_status[topGenes.donor_type_aab_status$adj.P.Val < 0.1,]),
          rownames(vwts.donor_type_aab_status$E)), ]

plot_gene_heatmap(
  counts.tmp, design=master.donor_type_aab_status, libID_col="library_id",
  order_by_var="donor_type_aab_status",
  color_by_var="donor_type_aab_status",
  my_var_colors=
    pal.donor_type_aab_status[levels(master.donor_type_aab_status$donor_type_aab_status)],
  add_legend=TRUE,
  ColSideLabs=NA,
  labCol=
    str_replace(master.donor_type_aab_status$batlab_id[
      order(master.donor_type_aab_status$donor_type_aab_status)], "-[0-9]F$", ""))

plot_gene_heatmap(
  counts.tmp, design=master.donor_type_aab_status, libID_col="library_id",
  color_by_var="donor_type_aab_status",
  my_var_colors=
    pal.donor_type_aab_status[levels(master.donor_type_aab_status$donor_type_aab_status)],
  col_dendro=TRUE,
  add_legend=TRUE,
  ColSideLabs=NA,
  labCol=
    str_replace(master.donor_type_aab_status$batlab_id, "-[0-9]F$", ""))

# Write list of significant + FC genes
write_sig_genes(
  topGenes.donor_type_aab_status,
  file_prefix="donor_type_aab_status.P0.1",
  method="combined",
  threshold_col="threshold")

rm_tmp(ask=FALSE)


##### read in MSigDB Hallmark interferon gene sets #####

gene_set.Hallmark_interferon_alpha <-
  read.table(
    paste0(
      "interferon_gene_modules/",
      "MSigDB_HALLMARK_INTERFERON_ALPHA_RESPONSE.txt"))[,1]

gene_set.Hallmark_interferon_gamma <-
  read.table(
    paste0(
      "interferon_gene_modules/",
      "MSigDB_HALLMARK_INTERFERON_GAMMA_RESPONSE.txt"))[,1]


##### read in GO associations and interferon terms #####

GO_BP.all_associations <-
  read.table(
    "interferon_gene_modules/goa_human.gaf",
    header=F, skip=34, sep="\t", quote="",
    col.names=
      c("DB", "DB Object ID", "DB Object Symbol", "Qualifier", "GO ID", "DB:Reference",
        "Evidence Code", "With (or) From", "Aspect", "DB Object Name", "DB Object Synonym",
        "DB Object Type", "Taxon", "Date", "Assigned By", "Annotation Extension",
        "Gene Product Form ID")) %>%
  standardize_dimnames()
for (i in 1:ncol(GO_BP.all_associations))
  GO_BP.all_associations[
    GO_BP.all_associations[,i] %in% "", i] <- NA
GO_BP.all_associations <-
  remove_all_NA_rowcols(GO_BP.all_associations)

## read in GO BP terms with "interferon" and "response" or "signaling pathway"
# excluded those with "biosynthetic process", "secretion", or "production"
GO_interferon_reponse_terms <-
  read.table(
    "GO_interferon_response_BP_terms.txt",
    header=F, sep="\t",
    col.names=c("go_id", "go_term", "go_synonym")) %>%
  remove_all_NA_rowcols()

## check for GO term overlap
GO_interferon_reponse_terms[
  GO_interferon_reponse_terms$go_id %nin% 
    GO_BP.all_associations$go_id,]

GO_BP.all_associations[
  GO_BP.all_associations$go_id %in% GO_interferon_reponse_terms$go_id,
  "db_object_symbol"] %>%
  unique() %>%
  sort()
# 239 genes

## read in GO BP terms with "interferon"
# and "response", "signaling pathway", "biosynthetic process", "secretion", or "production"
# excluded "viral suppression of..."
GO_interferon_terms <-
  read.table(
    "GO_interferon_response_production_secretion_BP_terms.txt",
    header=F, sep="\t",
    col.names=c("go_id", "go_term", "go_synonym")) %>%
  remove_all_NA_rowcols()

## check for GO term overlap
GO_interferon_terms[
  GO_interferon_terms$go_id %nin% 
    GO_BP.all_associations$go_id,]
# lots missing from the human GO biological processes; not sure why

GO_BP.interferon_term_related_genes <-
  GO_BP.all_associations[
    GO_BP.all_associations$go_id %in% GO_interferon_terms$go_id,
    "db_object_symbol"] %>%
  unique() %>%
  sort()
# 410 genes

rm_tmp(ask=FALSE)


##### plot eosinophil_purity_percent distributions #####

ggplot(
  master.inc_outliers,
  mapping=aes(x=eosinophil_purity_percent, color=donor_type_aab_status)) +
  geom_density(size=3) +
  scale_color_manual(values=pal.donor_type_aab_status)

ggplot(
  master.inc_outliers,
  mapping=aes(x=eosinophil_purity_percent, color=donor_type)) +
  geom_density(size=3) +
  scale_color_manual(values=pal.donor_type)


##### check for correlations of eosinophil_purity_percent with eosinophil-specific genes #####

gene_set.eosinophil_genes <-
  c("PRG2", "RNASE3", "RNASE2", "EPX", "CLC", "IL5RA", "CCR3")
gene_set.eosinophil_genes %in% rownames(vwts.all$E)

vwts.all.inc_outliers$E[
  match(gene_set.eosinophil_genes, rownames(vwts.all.inc_outliers$E)),] %>%
  na.omit() %>%
  # t() %>%
  apply(MARGIN=1, FUN=cor.test, y=master.inc_outliers$eosinophil_purity_percent) %>%
  lapply(function(x) x[["p.value"]]) %>%
  unlist() %>%
  p.adjust(method="BH")


##### check for correlations of eosinophil_purity_percent with DE genes from donor_type_aab_status.inc_outliers#####

# get the actual gene names
gene_set.DE.donor_type_aab_status.inc_outliers <-
  get_sig_genes(
    topGenes=topGenes.donor_type_aab_status.inc_outliers,
    method="combined",
    adj_p_cut=0.1, fc_cut=NULL)

# get the counts, and calculate correlations, within the healthy controls
counts.hc.tmp <-
  vwts.donor_type_aab_status.inc_outliers$E[
    match(
      gene_set.DE.donor_type_aab_status.inc_outliers,
      rownames(vwts.donor_type_aab_status.inc_outliers$E)),
    match(
      master.donor_type_aab_status.inc_outliers$library_id[
        master.donor_type_aab_status.inc_outliers$donor_type %in% "HC"],
      colnames(vwts.donor_type_aab_status.inc_outliers$E))]

counts.hc.tmp %>%
  na.omit() %>%
  apply(
    MARGIN=1, FUN=cor.test,
    y=master.donor_type_aab_status.inc_outliers$eosinophil_purity_percent[
      master.donor_type_aab_status.inc_outliers$donor_type %in% "HC"]) %>%
  # lapply(function(x) x[["estimate"]]) %>%
  # unlist()
  lapply(function(x) x[["p.value"]]) %>% # extract the p-value and adjust for multiple comparisons
  unlist() %>%
  p.adjust(method="BH") %>%
  sort()

# get the counts, and calculate correlations, within the healthy controls, or the at_risk patients
counts.at_risk.tmp <-
  vwts.donor_type_aab_status.inc_outliers$E[
    match(
      gene_set.DE.donor_type_aab_status.inc_outliers,
      rownames(vwts.donor_type_aab_status.inc_outliers$E)),
    match(
      master.donor_type_aab_status.inc_outliers$library_id[
        master.donor_type_aab_status.inc_outliers$donor_type %in% "at_risk"],
      colnames(vwts.donor_type_aab_status.inc_outliers$E))]

counts.at_risk.tmp %>%
  na.omit() %>%
  apply(
    MARGIN=1, FUN=cor.test,
    y=master.donor_type_aab_status.inc_outliers$eosinophil_purity_percent[
      master.donor_type_aab_status.inc_outliers$donor_type %in% "at_risk"]) %>%
  # lapply(function(x) x[["estimate"]]) %>%
  # unlist()
  lapply(function(x) x[["p.value"]]) %>% # extract the p-value and adjust for multiple comparisons
  unlist() %>%
  p.adjust(method="BH") %>%
  sort()


## repeat for the interferon genes only
gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers

counts.hc.tmp <-
  vwts.donor_type_aab_status.inc_outliers$E[
    match(
      gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers,
      rownames(vwts.donor_type_aab_status.inc_outliers$E)),
    match(
      master.donor_type_aab_status.inc_outliers$library_id[
        master.donor_type_aab_status.inc_outliers$donor_type %in% "HC"],
      colnames(vwts.donor_type_aab_status.inc_outliers$E))]

counts.hc.tmp %>%
  na.omit() %>%
  apply(
    MARGIN=1, FUN=cor.test,
    y=master.donor_type_aab_status.inc_outliers$eosinophil_purity_percent[
      master.donor_type_aab_status.inc_outliers$donor_type %in% "HC"]) %>%
  lapply(function(x) x[["estimate"]]) %>%
  unlist()
  # lapply(function(x) x[["p.value"]]) %>% # extract the p-value and adjust for multiple comparisons
  # unlist()
  # p.adjust(method="BH") %>%
  # sort()

counts.at_risk.tmp <-
  vwts.donor_type_aab_status.inc_outliers$E[
    match(
      gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers,
      rownames(vwts.donor_type_aab_status.inc_outliers$E)),
    match(
      master.donor_type_aab_status.inc_outliers$library_id[
        master.donor_type_aab_status.inc_outliers$donor_type %in% "at_risk"],
      colnames(vwts.donor_type_aab_status.inc_outliers$E))]

counts.at_risk.tmp %>%
  na.omit() %>%
  apply(
    MARGIN=1, FUN=cor.test,
    y=master.donor_type_aab_status.inc_outliers$eosinophil_purity_percent[
      master.donor_type_aab_status.inc_outliers$donor_type %in% "at_risk"]) %>%
  # lapply(function(x) x[["estimate"]]) %>%
  # unlist()
  lapply(function(x) x[["p.value"]]) %>% # extract the p-value and adjust for multiple comparisons
  unlist() %>%
  p.adjust(method="BH") %>%
  sort()

rm_tmp(ask=FALSE)


##### Set up and run limma for donor type, at_risk_vs_HC.inc_outliers #####

condition.tmp <-
  with(master.inc_outliers, donor_type %in% c("HC", "at_risk"))
master.at_risk_vs_HC.inc_outliers <-
  master.inc_outliers %>%
  filter(condition.tmp) %>%
  fix_factors()
DGECounts.at_risk_vs_HC.inc_outliers.tmp <-
  calc_norm_counts(
    counts=counts.final.inc_outliers,
    design=master.at_risk_vs_HC.inc_outliers, libID_col="library_id",
    min_cpm = 1, min_libs_perc = 0.15,
    return_DGEcounts=TRUE, group=master.at_risk_vs_HC.inc_outliers$patient_id)
master.at_risk_vs_HC.inc_outliers <-
  master.at_risk_vs_HC.inc_outliers[
    match(colnames(DGECounts.at_risk_vs_HC.inc_outliers.tmp),
          master.at_risk_vs_HC.inc_outliers[,"library_id"]),]

DesignMat.at_risk_vs_HC.inc_outliers <-
  model.matrix(
    ~ donor_type + sex,
    data=master.at_risk_vs_HC.inc_outliers)
vwts.at_risk_vs_HC.inc_outliers <-
  voomWithQualityWeights(
    DGECounts.at_risk_vs_HC.inc_outliers.tmp,
    design=DesignMat.at_risk_vs_HC.inc_outliers, plot=TRUE, span=0.1)
vfit.at_risk_vs_HC.inc_outliers <-
  lmFit(vwts.at_risk_vs_HC.inc_outliers) %>%
  eBayes()
topGenes.at_risk_vs_HC.inc_outliers <-
  topTable(vfit.at_risk_vs_HC.inc_outliers, coef = 2, number=Inf, sort.by="P")

# isolate genes that have an absolute fold change > 2^1.5  (absolute value of log2 fold change > 1.5)
# and an FDR < .1
pred.tmp <- master.at_risk_vs_HC.inc_outliers$donor_type
adj_factor.tmp <- if (is.numeric(pred.tmp)) 2*sd(pred.tmp) else 1
fc_cut.tmp <- log2(1.5) / adj_factor.tmp # need to adjust this for a reasonable range
p_cut.tmp <- 0.1
topGenes.at_risk_vs_HC.inc_outliers$threshold <-
  as.factor(abs(topGenes.at_risk_vs_HC.inc_outliers$logFC) > fc_cut.tmp &
              topGenes.at_risk_vs_HC.inc_outliers$adj.P.Val < p_cut.tmp)

## Heatmaps of counts by at_risk_vs_HC
plot_gene_heatmap(
  vwts_sig.at_risk_vs_HC.inc_outliers,
  design=master.at_risk_vs_HC.inc_outliers, libID_col="library_id",
  order_by_var="donor_type", color_by_var="donor_type",
  my_var_colors=
    pal.donor_type[
      levels(master.at_risk_vs_HC.inc_outliers$donor_type)],
  add_legend=TRUE, leg_x=0.92, leg_y=0.5)

# Write list of significant + FC genes
write_sig_genes(
  topGenes.at_risk_vs_HC.inc_outliers,
  "at_risk_vs_HC.inc_outliers",
  method="combined",
  adj_p_cut=p_cut.tmp, fc_cut=fc_cut.tmp, fc_adj_factor=adj_factor.tmp)

# output list of all genes in ensembl format for DAVID enrichment with background list
write_sig_genes(
  topGenes.at_risk_vs_HC.inc_outliers,
  file_prefix="at_risk_vs_HC.inc_outliers.ensgene",
  method=c("ranked_list"),
  input_type="symbol", output_type="ensgene")

rm_tmp(ask=FALSE)


##### Plot enriched biological processes in genes DE in at_risk_vs_HC.inc_outliers #####

# file read in here was generated using the DAVID enrichment tool, using the two files exported above as the target and background gene lists
DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1 <-
  read.table(
    "at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1.DAVID_GO_BP_DIRECT.background_all_genes_expressed.txt",
    sep="\t", header=T, quote="") %>%
  remove_all_NA_rowcols() %>%
  standardize_dimnames() %>%
  plyr::rename(replace=c("x"="percent"))
glimpse(DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1)
DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1$term_simple <-
  DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1$term %>%
  str_replace("GO:[0-9]+~", "")

ggplot(
  DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1[
    DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1$benjamini <= 0.1,],
  mapping=aes(x=reorder(term_simple, -benjamini), y=-log10(benjamini))) +
  scale_x_discrete(
    name=NULL,
    breaks=DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1$term_simple[
      DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1$benjamini <= 0.1],
    labels=
      str_wrap(DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1$term_simple[
        DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1$benjamini <= 0.1],
        width=25)) +
  geom_bar(stat='identity') +
  theme(axis.text.y=element_text(size=14)) +
  scale_y_continuous(
    expand = c(0, 0),
    limits =
      c(0,
        max(-log10(DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1$benjamini)) + 0.5)) +
  coord_flip() +
  labs(y="-log10 (adjusted p-value)")


##### Extract interferon-related genes from DAVID GO BP enrichment #####

# extract interferon-related genes
gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1 <-
  DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1 %>%
  dplyr::filter(str_detect(term, "interferon")) %>%
  dplyr::pull(genes) %>%
  str_split(", ") %>%
  unlist() %>%
  unique() %>%
  sort()

# write out gene list to a file
write.table(
  gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1,
  file=
    "at_risk_vs_HC.inc_outliers.FC1.5_and_P0.1.DAVID_GO_BP_DIRECT.background_all_genes_expressed.interferon_module_genes.txt",
  quote=FALSE, row.names=FALSE, col.names=FALSE)


##### Heatmaps of genes DE in at_risk_vs_HC.inc_outliers, in all patients, with interferon genes only ######

# the gene set with the significantly DE interferon-related genes
gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1

## for donor_type
# heatmap of genes with FDR < 0.1 and |FC| > 1.5, interferon-related genes only
plot_gene_heatmap(
  vwts.all.inc_outliers$E[
    na.omit(
      match(
        gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1,
        rownames(vwts.all.inc_outliers$E))),],
  design=master.inc_outliers, libID_col="library_id",
  color_by_var="donor_type",
  order_by_var="aab_sample_ordinal",
  my_var_colors=
    pal.donor_type[
      levels(master.inc_outliers$donor_type)],
  ColSideLabs=NA,
  col_dendro=TRUE,
  labCol=
    str_replace(master.inc_outliers$batlab_id[
      order(master.inc_outliers$aab_sample_ordinal)], "-[0-9]F$", ""),
  add_legend=TRUE, leg_x=0.92, leg_y=0.5)


## for donor_type, including some lower-count genes so that IFNL1 doesn't get filtered out
# generate a new version of vwts.all.inc_outliers that I can use
vwts.all.inc_outliers.14perc_min_cpm <-
  counts.final.inc_outliers %>%
  calc_norm_counts(
    design=master.inc_outliers, libID_col="library_id",
    min_cpm = 1, min_libs_perc = 0.14,
    return_DGEcounts=TRUE) %>%
  voomWithQualityWeights(plot=TRUE)

# heatmap of genes with FDR < 0.1 and |FC| > 1.5, interferon-related genes only
plot_gene_heatmap(
  vwts.all.inc_outliers.14perc_min_cpm$E[
    na.omit(
      match(
        gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers.genes_FC1.5_and_P0.1,
        rownames(vwts.all.inc_outliers.14perc_min_cpm$E))),],
  design=master.inc_outliers, libID_col="library_id",
  color_by_var="donor_type",
  order_by_var="aab_sample_ordinal",
  my_var_colors=
    pal.donor_type[
      levels(master.inc_outliers$donor_type)],
  ColSideLabs=NA,
  col_dendro=TRUE,
  labCol=
    str_replace(master.inc_outliers$batlab_id[
      order(master.inc_outliers$aab_sample_ordinal)], "-[0-9]F$", ""),
  add_legend=TRUE, leg_x=0.92, leg_y=0.5)


##### check for correlations of eosinophil_purity_percent with DE genes from at_risk_vs_HC.inc_outliers#####

# get the actual gene names
gene_set.DE.at_risk_vs_HC.inc_outliers <-
  get_sig_genes(
    topGenes=topGenes.at_risk_vs_HC.inc_outliers,
    method="combined",
    adj_p_cut=0.1, fc_cut=NULL)

# get the counts, and calculate correlations, within the healthy controls
counts.hc.tmp <-
  vwts.at_risk_vs_HC.inc_outliers$E[
    match(
      gene_set.DE.at_risk_vs_HC.inc_outliers,
      rownames(vwts.at_risk_vs_HC.inc_outliers$E)),
    match(
      master.at_risk_vs_HC.inc_outliers$library_id[
        master.at_risk_vs_HC.inc_outliers$donor_type %in% "HC"],
      colnames(vwts.at_risk_vs_HC.inc_outliers$E))]

counts.hc.tmp %>%
  na.omit() %>%
  apply(
    MARGIN=1, FUN=cor.test,
    y=master.at_risk_vs_HC.inc_outliers$eosinophil_purity_percent[
      master.at_risk_vs_HC.inc_outliers$donor_type %in% "HC"]) %>%
  # lapply(function(x) x[["estimate"]]) %>% 
  # unlist()
  lapply(function(x) x[["p.value"]]) %>% # extract the p-value and adjust for multiple comparisons
  unlist() %>%
  p.adjust(method="BH") %>%
  sort()

# get the counts, and calculate correlations, within the healthy controls, or the at_risk patients
counts.at_risk.tmp <-
  vwts.at_risk_vs_HC.inc_outliers$E[
    match(
      gene_set.DE.at_risk_vs_HC.inc_outliers,
      rownames(vwts.at_risk_vs_HC.inc_outliers$E)),
    match(
      master.at_risk_vs_HC.inc_outliers$library_id[
        master.at_risk_vs_HC.inc_outliers$donor_type %in% "at_risk"],
      colnames(vwts.at_risk_vs_HC.inc_outliers$E))]

counts.at_risk.tmp %>%
  na.omit() %>%
  apply(
    MARGIN=1, FUN=cor.test,
    y=master.at_risk_vs_HC.inc_outliers$eosinophil_purity_percent[
      master.at_risk_vs_HC.inc_outliers$donor_type %in% "at_risk"]) %>%
  # lapply(function(x) x[["estimate"]]) %>% 
  # unlist()
  lapply(function(x) x[["p.value"]]) %>% # extract the p-value and adjust for multiple comparisons
  unlist() %>%
  p.adjust(method="BH") %>%
  sort()


## repeat for the interferon genes only
gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers

counts.hc.tmp <-
  vwts.at_risk_vs_HC.inc_outliers$E[
    match(
      gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers,
      rownames(vwts.at_risk_vs_HC.inc_outliers$E)),
    match(
      master.at_risk_vs_HC.inc_outliers$library_id[
        master.at_risk_vs_HC.inc_outliers$donor_type %in% "HC"],
      colnames(vwts.at_risk_vs_HC.inc_outliers$E))]

counts.hc.tmp %>%
  na.omit() %>%
  apply(
    MARGIN=1, FUN=cor.test,
    y=master.at_risk_vs_HC.inc_outliers$eosinophil_purity_percent[
      master.at_risk_vs_HC.inc_outliers$donor_type %in% "HC"]) %>%
  lapply(function(x) x[["estimate"]]) %>%
  unlist()
# lapply(function(x) x[["p.value"]]) %>% # extract the p-value and adjust for multiple comparisons
# unlist()
# p.adjust(method="BH") %>%
# sort()

counts.at_risk.tmp <-
  vwts.at_risk_vs_HC.inc_outliers$E[
    match(
      gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.inc_outliers,
      rownames(vwts.at_risk_vs_HC.inc_outliers$E)),
    match(
      master.at_risk_vs_HC.inc_outliers$library_id[
        master.at_risk_vs_HC.inc_outliers$donor_type %in% "at_risk"],
      colnames(vwts.at_risk_vs_HC.inc_outliers$E))]

counts.at_risk.tmp %>%
  na.omit() %>%
  apply(
    MARGIN=1, FUN=cor.test,
    y=master.at_risk_vs_HC.inc_outliers$eosinophil_purity_percent[
      master.at_risk_vs_HC.inc_outliers$donor_type %in% "at_risk"]) %>%
  # lapply(function(x) x[["estimate"]]) %>%
  # unlist()
  lapply(function(x) x[["p.value"]]) %>% # extract the p-value and adjust for multiple comparisons
  unlist() %>%
  p.adjust(method="BH") %>%
  sort()

rm_tmp(ask=FALSE)


##### Set up and run limma for donor type, at_risk_vs_HC #####

condition.tmp <-
  with(master, donor_type %in% c("HC", "at_risk"))
master.at_risk_vs_HC <- fix_factors(filter(master, condition.tmp))
DGECounts.at_risk_vs_HC.tmp <-
  calc_norm_counts(
    counts=counts.final,
    design=master.at_risk_vs_HC, libID_col="library_id",
    min_cpm = 1, min_libs_perc = 0.15,
    return_DGEcounts=TRUE, group=master.at_risk_vs_HC$patient_id)
master.at_risk_vs_HC <-
  master.at_risk_vs_HC[
    match(colnames(DGECounts.at_risk_vs_HC.tmp),
          master.at_risk_vs_HC[,"library_id"]),]

DesignMat.at_risk_vs_HC <-
  model.matrix(
    ~ donor_type + sex,
    data=master.at_risk_vs_HC)
vwts.at_risk_vs_HC <-
  voomWithQualityWeights(
    DGECounts.at_risk_vs_HC.tmp,
    design=DesignMat.at_risk_vs_HC, plot=TRUE, span=0.1)
vfit.at_risk_vs_HC <-
  lmFit(vwts.at_risk_vs_HC) %>%
  eBayes()
topGenes.at_risk_vs_HC <-
  topTable(vfit.at_risk_vs_HC, coef = 2, number=Inf, sort.by="P")

# isolate genes that have an absolute fold change > 2^1.5  (absolute value of log2 fold change > 1.5)
# and an FDR < .1
pred.tmp <- master.at_risk_vs_HC$donor_type
adj_factor.tmp <- if (is.numeric(pred.tmp)) 2*sd(pred.tmp) else 1
fc_cut.tmp <- log2(1.5) / adj_factor.tmp # need to adjust this for a reasonable range
p_cut.tmp <- 0.1
topGenes.at_risk_vs_HC$threshold <-
  as.factor(
    abs(topGenes.at_risk_vs_HC$logFC) > fc_cut.tmp &
      topGenes.at_risk_vs_HC$adj.P.Val < p_cut.tmp)

## Heatmaps of counts by at_risk_vs_HC
plot_gene_heatmap(
  vwts_sig.at_risk_vs_HC,
  design=master.at_risk_vs_HC, libID_col="library_id",
  order_by_var="donor_type", color_by_var="donor_type",
  my_var_colors=
    pal.donor_type[
      levels(master.at_risk_vs_HC$donor_type)],
  add_legend=TRUE, leg_x=0.92, leg_y=0.5)

# Write list of significant + FC genes
write_sig_genes(
  topGenes.at_risk_vs_HC, "at_risk_vs_HC",
  method=c("ranked_list", "combined", "directional"),
  adj_p_cut=p_cut.tmp, fc_cut=fc_cut.tmp, fc_adj_factor=adj_factor.tmp)

# output list of all genes in ensembl format for DAVID enrichment with background list
write_sig_genes(
  topGenes.at_risk_vs_HC, "at_risk_vs_HC.ensgene",
  method="ranked_list",
  input_type="symbol", output_type="ensgene")

rm_tmp(ask=FALSE)


##### Extract interferon-related genes from DAVID GO BP enrichment #####

# extract interferon-related genes
gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.genes_FC1.5_and_P0.1 <-
  DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.genes_FC1.5_and_P0.1 %>%
  dplyr::filter(str_detect(term, "interferon")) %>%
  dplyr::pull(genes) %>%
  str_split(", ") %>%
  unlist() %>%
  unique() %>%
  sort()


##### Heatmaps of genes DE in at_risk_vs_HC, in all patients, with interferon genes only ######

# the gene set with the significant interferon-related genes
gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.genes_FC1.5_and_P0.1

## for donor_type
# heatmap of genes with FDR < 0.1 and |FC| > 1.5, interferon-related genes only
plot_gene_heatmap(
  vwts.all$E[
    na.omit(
      match(
        gene_set.interferon_module_genes.DAVID_GO_BP_DIRECT.background_all_genes_expressed.at_risk_vs_HC.genes_FC1.5_and_P0.1,
        rownames(vwts.all$E))),],
  design=master, libID_col="library_id",
  color_by_var="donor_type",
  order_by_var="aab_sample_ordinal",
  my_var_colors=
    pal.donor_type[
      levels(master$donor_type)],
  ColSideLabs=NA,
  col_dendro=TRUE,
  labCol=
    str_replace(master$batlab_id[
      order(master$aab_sample_ordinal)], "-[0-9]F$", ""),
  add_legend=TRUE, leg_x=0.92, leg_y=0.5)


##### summary statistics for transcriptomic data for donor_type_aab_status.inc_outliers #####

table(master.donor_type_aab_status.inc_outliers$donor_type_aab_status)
# matches table from Manuela Battaglia

contrasts.tmp <-
  cbind(c(0,-1,1,0),
        c(0,-1,0,1),
        c(0,0,-1,1))
lm.age_vs_donor_type_aab_status.inc_outliers <-
  lm(age ~ donor_type_aab_status, data=master.donor_type_aab_status.inc_outliers)
summary(lm.age_vs_donor_type_aab_status.inc_outliers)
contrast::contrast(
  lm.age_vs_donor_type_aab_status.inc_outliers,
  list(donor_type_aab_status="at_risk Aab+"),
  list(donor_type_aab_status="T1Dnew"))

# test for variation in age among patient groups using ANOVA
aov.age_vs_donor_type_aab_status.inc_outliers <-
  aov(age ~ donor_type_aab_status, data=master.donor_type_aab_status.inc_outliers)
summary.aov(aov.age_vs_donor_type_aab_status.inc_outliers)

## test for differences in sex by patient group
table(master.donor_type_aab_status[,c("sex", "donor_type_aab_status")]) %>%
  chisq.test()
  # fisher.test()
# p=0.8689 from chi-squared test

rm_tmp(ask=FALSE)


##### determine median gene expression of interferon signatures #####

for (i in names(gene_sets.T1D_neutrophils_interferon_signatures)) {
  master[,paste0("median_", i)] <-
    gene_set_median_count(
      gene_sets.T1D_neutrophils_interferon_signatures[[i]],
      vwts.all)
  
  master.inc_outliers[,paste0("median_", i)] <-
    gene_set_median_count(
      gene_sets.T1D_neutrophils_interferon_signatures[[i]],
      vwts.all.inc_outliers)
}

rm_tmp(ask=FALSE)


##### plot interferon gene signatures by donor_type #####

## Hallmark_interferon_alpha, with outliers
ggplot(
  master.inc_outliers,
  mapping=aes(
    x=donor_type,
    y=median_gene_set.Hallmark_interferon_alpha,
    fill=donor_type)) +
  geom_boxplot(outlier.color="white") +
  ggbeeswarm::geom_beeswarm(size=4, cex=2.5) +
  labs(x="Donor Status", y="Interferon gene set expression") +
  scale_fill_manual(
    values=
      c("gray50",
        average_colors(ggthemes::colorblind_pal()(4)[3:4]),
        ggthemes::colorblind_pal()(4)[2])) +
  guides(fill=FALSE) +
  theme(axis.text.x=element_text(angle=0, size=16, hjust=0.5))

## now test for a difference using Kruskal-Wallis test
kruskal.test(
  median_gene_set.Hallmark_interferon_alpha ~ donor_type,
  data=master.inc_outliers)

# Wilcoxon rank-sum test (Mann-Whitney U test) for pairwise differences
combn.tmp <- combn(levels(master.inc_outliers$donor_type), m=2)
wilcox.pairwise.median_gene_set.Hallmark_interferon_alpha_vs_donor_type.inc_outliers <- list()
for (i in 1:ncol(combn.tmp)) {
  data.tmp <-
    master.inc_outliers[
      master.inc_outliers$donor_type %in% combn.tmp[,i],]
  wilcox.pairwise.median_gene_set.Hallmark_interferon_alpha_vs_donor_type.inc_outliers[[
    paste(combn.tmp[,i], collapse="_")]] <-
    wilcox.test(
      median_gene_set.Hallmark_interferon_alpha ~ donor_type,
      data=data.tmp)
}


## Hallmark_interferon_gamma, with outliers
ggplot(
  master.inc_outliers,
  mapping=aes(
    x=donor_type,
    y=median_gene_set.Hallmark_interferon_gamma,
    fill=donor_type)) +
  geom_boxplot(outlier.color="white") +
  ggbeeswarm::geom_beeswarm(size=4, cex=2.5) +
  labs(x="Donor Status", y="Interferon gene set expression") +
  scale_fill_manual(
    values=
      c("gray50",
        average_colors(ggthemes::colorblind_pal()(4)[3:4]),
        ggthemes::colorblind_pal()(4)[2])) +
  guides(fill=FALSE) +
  theme(axis.text.x=element_text(angle=0, size=16, hjust=0.5))

## now test for a difference using Kruskal-Wallis test
kruskal.test(
  median_gene_set.Hallmark_interferon_gamma ~ donor_type,
  data=master.inc_outliers)

# Wilcoxon rank-sum test (Mann-Whitney U test) for pairwise differences
combn.tmp <- combn(levels(master.inc_outliers$donor_type), m=2)
wilcox.pairwise.median_gene_set.Hallmark_interferon_gamma_vs_donor_type.inc_outliers <- list()
for (i in 1:ncol(combn.tmp)) {
  data.tmp <-
    master.inc_outliers[
      master.inc_outliers$donor_type %in% combn.tmp[,i],]
  wilcox.pairwise.median_gene_set.Hallmark_interferon_gamma_vs_donor_type.inc_outliers[[
    paste(combn.tmp[,i], collapse="_")]] <-
    wilcox.test(
      median_gene_set.Hallmark_interferon_gamma ~ donor_type,
      data=data.tmp)
}


# Hallmark_interferon_alpha, no outliers
ggplot(
  master,
  mapping=aes(
    x=donor_type,
    y=median_gene_set.Hallmark_interferon_alpha,
    fill=donor_type)) +
  geom_boxplot(outlier.color="white") +
  ggbeeswarm::geom_beeswarm(size=4, cex=2.5) +
  labs(x="Donor Status", y="Interferon gene set expression") +
  scale_fill_manual(
    values=
      c("gray50",
        average_colors(ggthemes::colorblind_pal()(4)[3:4]),
        ggthemes::colorblind_pal()(4)[2])) +
  guides(fill=FALSE) +
  theme(axis.text.x=element_text(angle=0, size=16, hjust=0.5))

## Hallmark_interferon_gamma, no outliers
ggplot(
  master,
  mapping=aes(
    x=donor_type,
    y=median_gene_set.Hallmark_interferon_gamma,
    fill=donor_type)) +
  geom_boxplot(outlier.color="white") +
  ggbeeswarm::geom_beeswarm(size=4, cex=2.5) +
  labs(x="Donor Status", y="Interferon gene set expression") +
  scale_fill_manual(
    values=
      c("gray50",
        average_colors(ggthemes::colorblind_pal()(4)[3:4]),
        ggthemes::colorblind_pal()(4)[2])) +
  guides(fill=FALSE) +
  theme(axis.text.x=element_text(angle=0, size=16, hjust=0.5))
dev.off()

rm_tmp(ask=FALSE)

