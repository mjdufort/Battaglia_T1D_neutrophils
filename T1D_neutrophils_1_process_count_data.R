#### scripts for analysis of RNA-seq data for Vecchio, Lo Buono, ..., Battaglia. 2018. Abnormal neutrophil signature in the blood and pancreas of pre-symptomatic and symptomatic type 1 diabetes.
### this file includes scripts for upstream processing of RNA-seq data to remove problematic samples
## it should be used if starting with all samples
## it can be skipped if starting with data filtered to include only valid samples (e.g. data from GEO)

##### set up environment: load packages #####

## load general packages
library(xlsx) 
library(dplyr)
library(stringr)
library(ggplot2)
theme_set(
  theme_bw(20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black", fill=NA, size=1)))
library(ggthemes)

## load analysis-specific packages
library(gplots)
library(ggthemes)
library(VariantAnnotation)

# load packages with custom functions (available at https:://github.com/mjdufort)
library(RNAseQC)
library(countSubsetNorm)
library(miscHelpers)

# load maggritr (load last to ensure default pipe function selection)
library(magrittr)


##### read in count values #####

# read counts files
counts.merged <-
  merge(
    read.csv("counts/Tranche1_combined_counts.csv"),
    read.csv("counts/Tranche2_combined_counts.csv"),
    by="geneName")
rownames(counts.merged) <- counts.merged$geneName
counts.merged <- counts.merged[, colnames(counts.merged) != "geneName"]

# Reduce counts column names to sample identifier (4-digit number)
colnames(counts.merged) <-
  colnames(counts.merged) %>%
  str_extract(pattern="[0-9]+$") %>%
  make.unique(sep="_")

# Keep only protein coding genes with HGNC symbols
counts.merged.tmp <- counts.merged
genesHGNC <- get_HGNC(rownames(counts.merged.tmp), type="protein_coding")
counts.merged.tmp$HGNC.symbol <- genesHGNC

## sum counts.merged for duplicated HGNC symbols
# this also drops rows with HGNC.symbols==NA, which should include any genes that are not protein-coding
counts.merged.aggregated <-
  aggregate(counts.merged.tmp[, -ncol(counts.merged.tmp)],
            by=list(counts.merged.tmp$HGNC.symbol), sum)
rownames(counts.merged.aggregated) <- counts.merged.aggregated$Group.1
counts.merged.aggregated <-
  counts.merged.aggregated[,-which(colnames(counts.merged.aggregated) == "Group.1")]
dim(counts.merged.aggregated)
# 19751 genes, 48 libraries

rm_tmp(ask=FALSE)  # clean up unneeded variables


##### plot saturation curves of read counts #####

saturation.merged.aggregated <-
  estimate_saturation(
    counts.merged.aggregated, method="division", ndepths=30, min_counts=1,
    verbose=TRUE)
plot_saturation_curve(
  saturation.merged.aggregated, plot_lines=TRUE, plot_terminal_points=TRUE)
# looks good! reaching saturation at around 12,000-13,000 genes at about 6,000,000 reads in all libraries
# all libraries have lots of reads


##### read in sample and library prep data #####

sample_annotation.full <-
  read.xlsx(
    file="annotation_RNAseq_neutr_Battaglia.xlsx",
    sheetIndex=1, colIndex=1:60, check.names=FALSE) %>%
  remove_all_NA_rowcols()
glimpse(sample_annotation.full)

# drop the last two rows (which are not actually data rows)
sample_annotation.full <-
  sample_annotation.full[-((nrow(sample_annotation.full)-1):nrow(sample_annotation.full)),]
glimpse(sample_annotation.full)

# standardize column names
colnames(sample_annotation.full) <-
  colnames(sample_annotation.full) %>%
  str_replace_all("%", "_perc") %>%
  str_replace_all("#", "_abs") %>%
  str_replace_all("\\(U\\/L\\)", "") %>%
  str_replace_all("10x[36] ul", "") %>%
  str_replace_all("\\(mg\\/d?L\\)", "") %>%
  str_replace_all("(ng\\/)?ul", "") %>%
  str_replace_all("dl|fL|pg", "") %>%
  standardize_names()
glimpse(sample_annotation.full)

# exclude sample 4169, which is problematic
sample_annotation.full <-
  sample_annotation.full[
    sample_annotation.full$library_id != 4169,]

# standardize variable classes and levels
sample_annotation.full$donor_type[
  sample_annotation.full$donor_type == "FDR"] <-
  "at_risk"
sample_annotation.full$donor_type <-
  factor(sample_annotation.full$donor_type,
         levels=c("HC", "at_risk", "T1Dnew"))

sample_annotation.full$patient_id <-
  sample_annotation.full$patient_id %>%
  str_trim()
sample_annotation.full$library_id <-
  as.character(sample_annotation.full$library_id)

## check for duplicates in columns, to understand the data better
sapply(sample_annotation.full, function(x) any(duplicated(x, incomparables=NA)))
sample_annotation.full[
  duplicated(sample_annotation.full$dob) |
    duplicated(sample_annotation.full$dob, fromLast=TRUE),]
# essentially all of them have duplicates, except library_id
# patients are sometimes included more than once with
# different preps from the same visit (HC167), or different visit (preT1D071)

## create a duplicate of patient_id without the "-#F" stuff
sample_annotation.full$patient_id <-
  sample_annotation.full$patient_id %>%
  str_replace_all("\\-[0-9]F", "")

rm_tmp(ask=FALSE)

### read in additional library prep data, and merge it with existing sample_annotation
sample_annotation.updated.tmp <-
  read.xlsx(
    "annotation_RNAseq_neutr_Battaglia.xlsx",
    sheetName="ON SAMPLE",
    rowIndex=1:51) %>%
  remove_all_NA_rowcols() %>%
  standardize_dimnames() %>%
  dplyr::select(library_id, t1d_aab_on_sample, neu_purity, eosinophils) %>%
  plyr::rename(
    replace=c("t1d_aab_on_sample"="aab_sample_pret1d_only",
              "neu_purity"="neutrophil_purity_percent",
              "eosinophils"="eosinophil_purity_percent"))

sample_annotation.full <-
  sample_annotation.full %>%
  merge(sample_annotation.updated.tmp, by="library_id", all.x=TRUE)

# fill in HC and T1D for sample_annotation.full$aab_sample
sample_annotation.full$aab_sample_ordinal <-
  sample_annotation.full$aab_sample_pret1d_only
sample_annotation.full$aab_sample_ordinal[
  is.na(sample_annotation.full$aab_sample_ordinal) &
    sample_annotation.full$donor_type %in% "HC"] <- -1
sample_annotation.full$aab_sample_ordinal[
  is.na(sample_annotation.full$aab_sample_ordinal) &
    sample_annotation.full$donor_type %in% "T1Dnew"] <- 5
sample_annotation.full$aab_sample <-  
  c("HC", "AutoAb Neg", "AutoAb 1", "AutoAb 2", "AutoAb 3", "AutoAb 4", "T1D")[
    sample_annotation.full$aab_sample_ordinal + 2] %>%
  factor(
    levels=c("HC", "AutoAb Neg", "AutoAb 1", "AutoAb 2", "AutoAb 3", "AutoAb 4", "T1D"))

# create a new grouping with HC, at_risk_AbNeg, at_risk_AbPos, T1Dnew
sample_annotation.full$donor_type_aab_status <-
  sample_annotation.full$donor_type %>%
  as.character()
sample_annotation.full$donor_type_aab_status[
  with(sample_annotation.full,
       donor_type_aab_status %in% "at_risk" &
         aab_sample_ordinal %in% 0)] <-
  "at_risk Aab-"
sample_annotation.full$donor_type_aab_status[
  with(sample_annotation.full,
       donor_type_aab_status %in% "at_risk" &
         aab_sample_ordinal %in% 1:4)] <-
  "at_risk Aab+"
sample_annotation.full$donor_type_aab_status <-
  sample_annotation.full$donor_type_aab_status %>%
  factor(
    levels=c("HC", "at_risk Aab-", "at_risk Aab+", "T1Dnew"))

rm_tmp(ask=FALSE)


##### create color palettes for use in plots and such #####

pal.donor_type_aab_status <-
  colorblind_pal()(4)[c(1,3,4,2)] %>%
  setNames(levels(sample_annotation.full$donor_type_aab_status))
pal.donor_type <-
  c(colorblind_pal()(4)[1],
    average_colors(colorblind_pal()(4)[3:4]),
    colorblind_pal()(4)[2]) %>%
  setNames(levels(sample_annotation.full$donor_type))


##### read in library metrics #####

# read the raw data, merge them, fix names
metrics.merged <-
  rbind(
    read.csv("Tranche1_combined_metrics.csv"),
    read.csv("Tranche2_combined_metrics.csv")) %>%
  standardize_dimnames()
# glimpse(metrics.merged)

# remove "sample" from library IDs
metrics.merged$libid <-
  metrics.merged$libid %>%
  str_extract(pattern="[0-9]+")

# convert %s to decimal values
for (i in 2:ncol(metrics.merged)) {
  if (!is.numeric(metrics.merged[,i])) {
    if (sum(str_detect(metrics.merged[,i], "%")) > 0) {
      metrics.merged[,i] <- as.numeric(str_replace(metrics.merged[,i], "%", "")) / 100
    } else metrics.merged[,i] <- as.numeric(metrics.merged[,i])
  }
}

rm_tmp(ask=FALSE)


##### remove libraries with poor quality metrics #####

#Plot total counts
metrics.merged <- arrange(metrics.merged, fastq_total_reads)

plot_read_counts(
  metrics.merged,
  n_lowcount=30,
  id_col="libid")
# 1 "low"-count library, with only 18 million reads; keep them all!

## plot total reads, mapped_reads_w_dups, median_cv_coverage all against each other
plot_metrics(
  metrics.merged, metrics.libID_col="libid",
  design=sample_annotation.full, design.libID_col="library_id",
  point_size=3,
  plotdims=c(11,9))
# 3481 has a relatively low % alignment (77%), but probably OK to include

# problematic libraries using standard thresholds: 3481
lib.bad.tmp <-
  metrics.merged$libid[
    with(metrics.merged,
         (fastq_total_reads < 5e6) | 
           (mapped_reads_w_dups < 0.8) |
           (median_cv_coverage > 1.0))]

# inspect metrics and sample data for problematic libraries
metrics.merged[metrics.merged$libid %in% lib.bad.tmp,]
sample_annotation.full[sample_annotation.full$library_id %in% lib.bad.tmp,]

# check sample percent duplication
metrics.merged[
  order(metrics.merged$percent_duplication, decreasing=F),
  c("libid", "percent_duplication")]
# pretty high, but that's probably due to the high read counts

## Make quality control cuts based on slightly modified metrics
metrics.merged.qc <-
  metrics.merged[
    with(metrics.merged,
         fastq_total_reads > 5e6 & mapped_reads_w_dups > 0.75 & median_cv_coverage < 1.0),]
nrow(metrics.merged.qc)  # 48 libraries (QC removed 0, or 0% of libraries)

# Remove from counts data libraries that fail QC cuts, or aren't in the annotation object
counts.merged.aggregated.qc <-
  counts.merged.aggregated[
    , colnames(counts.merged.aggregated) %in% metrics.merged.qc$libid]
counts.merged.aggregated.qc <-
  counts.merged.aggregated.qc[,order(colnames(counts.merged.aggregated.qc))]
counts.merged.aggregated.qc <-
  counts.merged.aggregated.qc[
    ,colnames(counts.merged.aggregated.qc) %in% sample_annotation.full$library_id]
# dropped one library with no annotation data

# Remove from sample annotation data libraries that fail QC cuts
sample_annotation.full.qc <-
  sample_annotation.full[
    sample_annotation.full$library_id %in%
      colnames(counts.merged.aggregated.qc),]
sample_annotation.full.qc <-
  sample_annotation.full.qc[order(sample_annotation.full.qc$library_id),]

rm_tmp(ask=FALSE)


##### check sample annotated sex against inferred sex from RNAseq data #####

## calculate log ratio of reads mapping to X and to Y chromosome, for each library
logXY.tmp <- logXYratio(counts.merged.aggregated.qc, gene_ID="symbol")

# examine distribution
ggplot(data=data.frame(logXYratio=logXY.tmp), mapping=aes(x=logXYratio) )+
  geom_histogram(colour="black", fill="gray60")
sort(logXY.tmp) # obvious break between 6 and 10

# examine list of libraries by assigned sex and logXYratio
cbind(sample_annotation.full.qc$sex,
      logXY.tmp[match(sample_annotation.full.qc$library_id, names(logXY.tmp))])
# looks like libraries 4174 is supposed to be male, but actually female
# looks like libraries 4175 is supposed to be female, but actually male
# library 4167 is supposed to be the same individual as 4175, not 4174

# assign inferred sex based on total number of reads mapping to X and Y chromosomes, for each library
sample_annotation.full.qc$sex_by_rna <-
  ifelse(logXY.tmp[sample_annotation.full.qc$library_id] > 8, "F", "M")
sample_annotation.full.qc[,c("library_id","sex","sex_by_rna")]

# check that inferred sex is same for all libraries for each patient
table(sample_annotation.full.qc[,c("patient_id", "sex_by_rna")])
which(rowSums(table(sample_annotation.full.qc[,c("patient_id", "sex_by_rna")])==0)!=1)
# 1 patient has 2 libraries inferred to be male and female

# plot them
data.tmp <- data.frame(
  sex=sample_annotation.full.qc$sex[
    match(names(logXY.tmp), sample_annotation.full.qc$library_id)],
  logXYratio=logXY.tmp)

# plot logXYratio by annotated sex
ggplot(data=data.tmp, mapping=aes(x=logXYratio, fill=sex)) +
  geom_histogram(position="dodge") +
  scale_fill_manual(values=c("red", "blue"))
# note apparent swap of 2 libraries

# check that inferred sex matches sex from database
all(sample_annotation.full.qc$sex == sample_annotation.full.qc$sex_by_rna)
# no; identify the ones that look wrong
sample_annotation.full.qc[
  which(sample_annotation.full.qc$sex != sample_annotation.full.qc$sex_by_rna),]

# check UTY expression for those libraries with mis-matched sex
counts.merged.aggregated.qc[
  "UTY",
  sample_annotation.full.qc$library_id[
    which(sample_annotation.full.qc$sex !=
            sample_annotation.full.qc$sex_by_rna)],]

# extract library names to include in tests
libs.test_sex.tmp <- colnames(counts.merged.aggregated)

counts.test_sex.tmp <-
  data.frame(
    libid=libs.test_sex.tmp,
    sex=factor(sample_annotation.full$sex[
      match(libs.test_sex.tmp, sample_annotation.full$library_id)]),
    UTY=log(unlist(counts.merged[match("ENSG00000229807", rownames(counts.merged)), libs.test_sex.tmp])+0.5),
    XIST=log(unlist(counts.merged[match("ENSG00000183878", rownames(counts.merged)), libs.test_sex.tmp])+0.5))

counts.test_sex.tmp <-
  counts.test_sex.tmp %>%
  group_by(sex) %>%
  mutate(UTY_outlier = ifelse(is_outlier(UTY), libid, NA)) %>%
  mutate(XIST_outlier = ifelse(is_outlier(XIST), libid, NA))

## scatterplot of points on UTY and XIST to identify sex (with points colored by sex)
ggplot(data=counts.test_sex.tmp, aes(x=UTY, y=XIST, colour=sex)) +
  geom_point(position=position_jitter(w=0.02,h=0.02), size=3) +
  geom_text(aes(label=UTY_outlier), colour="black", vjust=-1.3) +
  geom_text(aes(label=XIST_outlier), colour="black", vjust=-1.3) +
  scale_color_manual(values=c("red", "blue")) +
  labs(x="log10 UTY expression", y="log10 XIST expression")
# 4174 and 4175 end up with the wrong group, as they do with logXYratio

rm_tmp(ask=FALSE)  # remove unneeded variables


##### check identity of patients using SNP variants called from RNAseq reads #####
## used samtools mpileup on BAM alignment files against a BED file with common variants to get BCF files
## then used bcftools call to convert bcf files to vcf files

## import data
variant_data <-
  readVcf(
    "T1D_neutrophils_4167_4174_4175_mpileup.vcf",
    genome="hg38")

hist(qual(variant_data), breaks=50) # plot histogram of quality scores

variant_data.summary <- snpSummary(variant_data) # generate summary of variants

variant_data.snpmat <- genotypeToSnpMatrix(variant_data) # convert to SNP matrix

# filter SNP matrix to include high-quality variants, substitutions only, variable among samples
variant_data.snpmat.filtered <-
  variant_data.snpmat$genotypes[
    ,(!variant_data.snpmat$map$ignore) & # remove NAs
      (qual(variant_data) > 500) & # restrict to quality score > 500
      (isSubstitution(variant_data))] # restrict to substitutions
class(variant_data.snpmat.filtered) <- "matrix"
variant_data.snpmat.filtered.variable <-
  variant_data.snpmat.filtered[
    ,apply(variant_data.snpmat.filtered, 2, function(x) length(unique(x)) > 1)] # remove variants with same genotype in all samples
dist(variant_data.snpmat.filtered, method="manhattan")
n_samples.tmp <- nrow(variant_data.snpmat.filtered)
variant_data.similarity_mat <-
  matrix(NA, nrow=n_samples.tmp, ncol=n_samples.tmp,
         dimnames=list(rownames(variant_data.snpmat.filtered.variable),
                       rownames(variant_data.snpmat.filtered.variable)))
for (i in 1:n_samples.tmp)
  for (j in 1:n_samples.tmp) {
    variant_data.similarity_mat[i,j] <-
      sum(variant_data.snpmat.filtered.variable[i,]==
            variant_data.snpmat.filtered.variable[j,]) /
      ncol(variant_data.snpmat.filtered.variable)
  }
    
variant_data.similarity_mat

# 0.814 of variable variants shared by 4167 and 4174
# 0.063 of variable variants shared by 4167 and 4175
# 0.068 of variable variants shared by 4174 and 4175

# plot a venn diagram of the shared variants
venndata.variant_similarity <-
  list(
    "4167"=paste(seq_len(nrow(variant_data.snpmat.filtered.variable)), variant_data.snpmat.filtered.variable["Sample_4167.bam",], sep="_"),
    "4174"=paste(seq_len(nrow(variant_data.snpmat.filtered.variable)), variant_data.snpmat.filtered.variable["Sample_4174.bam",], sep="_"),
    "4175"=paste(seq_len(nrow(variant_data.snpmat.filtered.variable)), variant_data.snpmat.filtered.variable["Sample_4175.bam",], sep="_"))
venndata.variant_similarity <-
  lapply(venndata.variant_similarity,
         function(x) paste(1:length(venndata.variant_similarity[[1]]), x, sep="_"))
  
venn(data=venndata.variant_similarity)

# feel fairly confident that 4167 and 4174 are the same person
# that means the annotation for sample 4175 belongs with sample 4174
# but what about sample 4175? is it just a one-for-one swap with sample 4174? hard to tell

## swap annotation data for 4174/4175, and drop 4175
## drop 4174 later, because we only want one sample from each patient, and this one is suspect

rm_tmp(ask=FALSE)


##### create a corrected version of the sample annotation object, and rerun the steps above #####

sample_annotation.full.corrected.qc <- sample_annotation.full.qc

# remove sample 4174 from annotation, and rename 4175 to 4174
sample_annotation.full.corrected.qc <- sample_annotation.full.corrected.qc[
  sample_annotation.full.corrected.qc$library_id != "4174",]
sample_annotation.full.corrected.qc$library_id[
  sample_annotation.full.corrected.qc$library_id == "4175"] <- "4174"

# also drop samples 4186 and 4187 from annotation, as they were duplicates purified using a different method
sample_annotation.full.corrected.qc <- sample_annotation.full.corrected.qc[
  sample_annotation.full.corrected.qc$library_id %nin% c("4186", "4187"),]

# remove the dropped libraries from my counts and metrics objects
counts.merged.aggregated.corrected.qc <-
  counts.merged.aggregated.qc[
    ,match(sample_annotation.full.corrected.qc$library_id,
           colnames(counts.merged.aggregated.qc))]
metrics.merged.corrected.qc <-
  metrics.merged.qc[
    match(sample_annotation.full.corrected.qc$library_id, metrics.merged.qc$libid),]


##### repeat basic checks for problematic data #####

## calculate total number of reads mapping to X and to Y chromosome, for each library
logXY.tmp <- logXYratio(counts.merged.aggregated.corrected.qc, gene_ID="symbol")
hist(logXY.tmp, breaks=40, main="log-ratio of X reads to Y reads")
sort(logXY.tmp) # obvious break between 5 and 9

# examine list of libraries by assigned sex and logXYratio
cbind(sample_annotation.full.corrected.qc$sex,
      logXY.tmp[match(sample_annotation.full.corrected.qc$library_id, names(logXY.tmp))])
# pattern is now consistent

# assign sex based on total number of reads mapping to X and Y chromosomes, for each library
sample_annotation.full.corrected.qc$sex_by_rna <-
  ifelse(logXY.tmp[sample_annotation.full.corrected.qc$library_id] > 7, "F", "M")
sample_annotation.full.corrected.qc[,c("library_id","sex","sex_by_rna")]

# check that inferred sex is same for all libraries for each patient
table(sample_annotation.full.corrected.qc[,c("patient_id", "sex_by_rna")])
which(rowSums(table(sample_annotation.full.corrected.qc[,c("patient_id", "sex_by_rna")])==0)!=1)
# now they all match

# plot them
data.tmp <-
  data.frame(
    sex=sample_annotation.full.corrected.qc$sex[
      match(names(logXY.tmp), sample_annotation.full.corrected.qc$library_id)],
    logXYratio=logXY.tmp)

ggplot(data=data.tmp, mapping=aes(x=logXYratio, fill=sex)) +
  geom_histogram(position="dodge") +
  scale_fill_manual(values=c("red", "blue"))

# check that inferred sex matches sex from database
all(sample_annotation.full.corrected.qc$sex ==
      sample_annotation.full.corrected.qc$sex_by_rna)
# looks good now!

rm_tmp(ask=FALSE)  # remove unneeded variables


##### run PCA on (normalized) full data set with duplicates included #####
# run PCA here because we dropped one questionable library above

# Remove libraries not in sample annotation object,
# filter out genes that have a count of at least one in < 15% of libraries,
# and normalize counts across libraries
counts.merged.aggregated.corrected.qc.normalized <-
  calc_norm_counts(
    counts=counts.merged.aggregated.corrected.qc,
    min_cpm=1,
    design=sample_annotation.full.corrected.qc,
    libID_col="library_id")
# glimpse(counts.merged.aggregated.corrected.qc.normalized)

## run PCA on filtered, normalized counts

# run PCA
pcaAll.full.corrected.qc <-
  calc_PCAs(counts.merged.aggregated.corrected.qc.normalized, log2_transform=TRUE)

# scree plot, the number of informative PCs = elbow
plot(pcaAll.full.corrected.qc, type="l")
# 2-3 informative PCs

# attach sample info to PC scores for ease of plotting
scores_sample_annotation_pcaAll.full.corrected.qc <-
  merge(sample_annotation.full.corrected.qc,
        as.data.frame(pcaAll.full.corrected.qc$x), by.x = "library_id", by.y="row.names")

## plots of PCs1 1-3, with no color
plot_PCAs(
  scores_sample_annotation_pcaAll.full.corrected.qc,
  pvars.labs=pcaAll.full.corrected.qc$pvars.labs, PCs = 1:3,
  plotdims=c(7.3,6))

## plot PC1 and PC2, with text for the duplicated samples
ggplot(scores_sample_annotation_pcaAll.full.corrected.qc,
       mapping=aes(x=PC1, y=PC2)) +
  geom_point(size=3) +
  geom_text(data=scores_sample_annotation_pcaAll.full.corrected.qc[
    str_detect(scores_sample_annotation_pcaAll.full.corrected.qc$patient_id, "preT1D071"),],
    mapping=aes(label=patient_id), nudge_y=-5)


##### drop duplicate sample and re-normalize #####

lib.dup.tmp <- c("4174")

sample_annotation.full.corrected.qc.ex_dups <-
  sample_annotation.full.corrected.qc[
    sample_annotation.full.corrected.qc$library_id %nin% lib.dup.tmp,]
metrics.merged.corrected.qc.ex_dups <-
  metrics.merged.corrected.qc[
    metrics.merged.corrected.qc$libid %nin% lib.dup.tmp,]
counts.merged.aggregated.corrected.qc.ex_dups <-
  counts.merged.aggregated.corrected.qc[
    ,colnames(counts.merged.aggregated.corrected.qc) %nin% lib.dup.tmp]


##### run PCA on (normalized) full data set to check for bias and batch effects, with outliers included #####
# run PCA here because we dropped one questionable library

# Remove libraries not in sample annotation object,
# filter out genes that have a count of at least one in < 15% of libraries,
# and normalize counts across libraries
counts.merged.aggregated.corrected.qc.ex_dups.normalized <-
  calc_norm_counts(
    counts=counts.merged.aggregated.corrected.qc.ex_dups,
    min_cpm=1,
    design=sample_annotation.full.corrected.qc.ex_dups,
    libID_col="library_id")
# glimpse(counts.merged.aggregated.corrected.qc.ex_dups.normalized)

# run PCA
# only need pcaAll for the scores, the sds, and the plot
pcaAll.full.corrected.qc.ex_dups <-
  calc_PCAs(counts.merged.aggregated.corrected.qc.ex_dups.normalized,
            log2_transform=TRUE)

# scree plot, the number of informative PCs = elbow
plot(pcaAll.full.corrected.qc.ex_dups, type="l")
# 2-4 informative PCs

# attach sample info to PC scores for ease of plotting
scores_sample_annotation_pcaAll.full.corrected.qc.ex_dups <-
  merge(sample_annotation.full.corrected.qc.ex_dups,
        as.data.frame(pcaAll.full.corrected.qc.ex_dups$x),
        by.x = "library_id", by.y="row.names")

### Color plots of PCs by different variables

## plots of PCs1 1-3, with no color
plot_PCAs(
  scores_sample_annotation_pcaAll.full.corrected.qc.ex_dups,
  pvars.labs=pcaAll.full.corrected.qc.ex_dups$pvars.labs, PCs = 1:3,
  plotdims=c(7.3,6))

## plots of PC1s 1-3, colored by donor_type (T1D, at_risk, HC)
plot_PCAs(
  scores_sample_annotation_pcaAll.full.corrected.qc.ex_dups,
  pvars.labs=pcaAll.full.corrected.qc.ex_dups$pvars.labs, PCs = 1:3,
  color_by_var="donor_type",
  my_cols=colorblind_pal()(3),
  plotdims=c(9,6))
# HC look very different from some T1D and at_risk, and a little different from the rest of the T1D and at_risk

## plots of PC1s 1-3, colored by age
plot_PCAs(
  scores_sample_annotation_pcaAll.full.corrected.qc.ex_dups,
  pvars.labs=pcaAll.full.corrected.qc.ex_dups$pvars.labs, PCs = 1:3,
  color_by_var="age_at_sample_coll",
  color_var_lab="Age at\ncollection",
  plotdims=c(9,6))
# possibly some variation with age on PC1, but not strong


## quantify correlations of annotation variables with PCs, and output a heatmap of those correlations
sample_annotation.metrics.full.corrected.qc.ex_dups <-
  merge(sample_annotation.full.corrected.qc.ex_dups,
        metrics.merged.corrected.qc.ex_dups,
        by.x="library_id", by.y="libid")

cor_all_clinvars_vs_pcaAll.counts.shared.qc.filtered_genes <-
  calc_PCcors(
    pcaAll.full.corrected.qc.ex_dups,
    sample_annotation.metrics.full.corrected.qc.ex_dups,
    "id_col"="library_id")

plot_PCcor_heatmap(
  cor_all_clinvars_vs_pcaAll.counts.shared.qc.filtered_genes,
  plotdims=c(12,8))
# nothing really jumps out as being strongly correlated
# other than visit_id variables and others that capture duplicate visits
# PC2 appears to be weakly correlated with lots of things related to study or metrics or something
# donor type appears to segregate most on PC1 and PC4

# plot PCAs colored by donor type, with PCs where there is the greatest separation
plot_PCAs(
  scores_sample_annotation_pcaAll.full.corrected.qc.ex_dups,
  pvars.labs=pcaAll.full.corrected.qc.ex_dups$pvars.labs, PCs = c(1,4),
  color_by_var="donor_type",
  my_cols=pal.donor_type,
  color_var_lab="Donor\nType",
  plotdims=c(8.35,6))

plot_PCAs(
  scores_sample_annotation_pcaAll.full.corrected.qc.ex_dups,
  pvars.labs=pcaAll.full.corrected.qc.ex_dups$pvars.labs, PCs = c(1,4),
  color_by_var="donor_type_aab_status",
  my_cols=pal.donor_type_aab_status,
  color_var_lab="Donor\nType",
  plotdims=c(9,6))

rm_tmp(ask=FALSE)


##### drop outliers, re-normalize, and re-run PCA #####

## plot individual variables vs. PC scores for interesting variables
ggplot(data=scores_sample_annotation_pcaAll.full.corrected.qc.ex_dups,
       mapping=aes(x=PC1, y=PC2)) +
  geom_point() +
  geom_text(mapping=aes(label=library_id), nudge_y=-3)

## some points that have really high values for PC1, but not sure why
# 3474, 4177, 4183
# 3474 is a preT1D subject with anomalous data
# 4177 is a preT1D subject with anomalous data
# 4183 is a T1D subject with anomalous data; drop it, as subject had a chronic inflammatory condition
# there are NOT the samples with worse scores in some of the metrics (3469, 3477, 3480)
outlier_libs.tmp <- c("3474", "4177", "4183")

sample_annotation.full.corrected.qc.ex_dups.ex_outliers <-
  sample_annotation.full.corrected.qc.ex_dups[
    sample_annotation.full.corrected.qc.ex_dups$library_id %nin% outlier_libs.tmp,]
metrics.merged.corrected.qc.ex_dups.ex_outliers <-
  metrics.merged.corrected.qc.ex_dups[
    metrics.merged.corrected.qc.ex_dups$libid %nin% outlier_libs.tmp,]
counts.merged.aggregated.corrected.qc.ex_dups.ex_outliers <-
  counts.merged.aggregated.corrected.qc.ex_dups[
    ,colnames(counts.merged.aggregated.corrected.qc.ex_dups) %nin% outlier_libs.tmp]

# Remove libraries not in sample annotation object,
# filter out genes that have a count of at least one in < 15% of libraries,
# and normalize counts across libraries
counts.merged.aggregated.corrected.qc.ex_dups.ex_outliers.normalized <-
  calc_norm_counts(
    counts=counts.merged.aggregated.corrected.qc.ex_dups.ex_outliers,
    min_cpm=1,
    design=sample_annotation.full.corrected.qc.ex_dups.ex_outliers,
    libID_col="library_id")
# glimpse(counts.merged.aggregated.corrected.qc.ex_dups.ex_outliers.normalized)

# run PCA
# only need pcaAll for the scores, the sds, and the plot
pcaAll.full.corrected.qc.ex_dups.ex_outliers <-
  calc_PCAs(counts.merged.aggregated.corrected.qc.ex_dups.ex_outliers.normalized,
            log2_transform=TRUE)

# scree plot, the number of informative PCs = elbow
plot(pcaAll.full.corrected.qc.ex_dups.ex_outliers, type="l")
# ??? informative PCs

# attach sample info to PC scores for ease of plotting
scores_sample_annotation_pcaAll.full.corrected.qc.ex_dups.ex_outliers <-
  merge(sample_annotation.full.corrected.qc.ex_dups.ex_outliers,
        as.data.frame(pcaAll.full.corrected.qc.ex_dups.ex_outliers$x),
        by.x = "library_id", by.y="row.names")

## plots of PCs1 1-3, with no color
plot_PCAs(
  scores_sample_annotation_pcaAll.full.corrected.qc.ex_dups.ex_outliers,
  pvars.labs=pcaAll.full.corrected.qc.ex_dups.ex_outliers$pvars.labs, PCs = 1:3,
  plotdims=c(7.3,6))

## plots of PC1s 1-3, colored by donor_type (T1D, at_risk, HC)
plot_PCAs(
  scores_sample_annotation_pcaAll.full.corrected.qc.ex_dups.ex_outliers,
  pvars.labs=pcaAll.full.corrected.qc.ex_dups.ex_outliers$pvars.labs, PCs = 1:3,
  color_by_var="donor_type",
  my_cols=colorblind_pal()(3),
  plotdims=c(9,6))
# HC look very different from some T1D and at_risk, and a little different from the rest of the T1D and at_risk

## quantify correlations of annotation variables with PCs, and output a heatmap
sample_annotation.metrics.full.corrected.qc.ex_dups.ex_outliers <-
  merge(sample_annotation.full.corrected.qc.ex_dups.ex_outliers,
        metrics.merged.corrected.qc.ex_dups.ex_outliers,
        by.x="library_id", by.y="libid")

cor_all_clinvars_vs_pcaAll.counts.shared.qc.filtered_genes.ex_outliers <-
  calc_PCcors(
    pcaAll.full.corrected.qc.ex_dups.ex_outliers,
    sample_annotation.metrics.full.corrected.qc.ex_dups.ex_outliers,
    "id_col"="library_id")

plot_PCcor_heatmap(
  cor_all_clinvars_vs_pcaAll.counts.shared.qc.filtered_genes.ex_outliers,
  plotdims=c(12,8))
# nothing really jumps out as being strongly correlated
# other than visit_id variables and others that capture duplicate visits
# PC2 appears to be weakly correlated with lots of things related to study or metrics or something
# donor type appears to segregate most on PC1 and PC5

# plot PCAs colored by donor type, using PCs where there is the greatest separation
plot_PCAs(
  scores_sample_annotation_pcaAll.full.corrected.qc.ex_dups.ex_outliers,
  pvars.labs=pcaAll.full.corrected.qc.ex_dups.ex_outliers$pvars.labs, PCs = c(1,5),
  color_by_var="donor_type",
  my_cols=pal.donor_type,
  color_var_lab="Donor\nType",
  plotdims=c(8.35,6))

plot_PCAs(
  scores_sample_annotation_pcaAll.full.corrected.qc.ex_dups.ex_outliers,
  pvars.labs=pcaAll.full.corrected.qc.ex_dups.ex_outliers$pvars.labs, PCs = c(1,5),
  color_by_var="donor_type_aab_status",
  my_cols=pal.donor_type_aab_status,
  color_var_lab="Donor\nType",
  plotdims=c(9,6))

rm_tmp(ask=FALSE)


##### plot various clinical variables vs. age, to understand age-related variation in disease #####

## examine age distribution by donor type
ggplot(data=sample_annotation.full.corrected.qc.ex_dups.ex_outliers,
       mapping=aes(x=donor_type, y=age_at_sample_coll)) +
  scale_x_discrete(limits=c("HC", "at_risk", "T1Dnew")) +
  geom_boxplot(outlier.color="white") +
  geom_jitter(width=0.1) +
  labs(x="Donor type", y="Age")

ggplot(data=sample_annotation.full.corrected.qc.ex_dups.ex_outliers,
       mapping=aes(x=donor_type_aab_status, y=age_at_sample_coll)) +
  geom_boxplot(fill="gray70", outlier.color="white") +
  geom_jitter(size=2, width=0.1) +
  labs(x="Donor type", y="Age")
 
## plot density of age by donor type
ggplot(data=sample_annotation.full.corrected.qc.ex_dups.ex_outliers) +
  geom_density(mapping=aes(x=age_at_sample_coll, group=donor_type, color=donor_type), size=1,
               position="identity") +
  scale_color_manual(values=colorblind_pal()(3)) +
  labs(x="Age")


##### check for variation among patient groups in neutrophil purity % #####

# ANOVA of neutrophil_purity_percent by donor_type_aab_status
lm.neutrophil_purity_percent_vs_donor_type_aab_status.inc_outliers <-
  lm(neutrophil_purity_percent ~ donor_type_aab_status,
     data=sample_annotation.final.inc_outliers)
aov(lm.neutrophil_purity_percent_vs_donor_type_aab_status.inc_outliers) %>%
  summary()

# plot of neutrophil_purity_percent by donor_type_aab_status
ggplot(
  sample_annotation.final.inc_outliers,
  mapping=
    aes(x=donor_type_aab_status,
        y=neutrophil_purity_percent)) +
  geom_violin() +
  ggbeeswarm::geom_beeswarm(size=3, cex=2) 


##### write out counts file for samples included in downstream analyses #####

counts.merged.for_GEO <-
  counts.merged[
    , match(sample_annotation.final.inc_outliers$library_id,
            colnames(counts.merged))]
colnames(counts.merged.for_GEO) <-
  paste0("sample_", colnames(counts.merged.for_GEO))
  
write.table(
  counts.merged.for_GEO,
  file="data_output/raw_counts_T1D_neutrophils.txt",
  sep="\t", row.names=TRUE, quote=FALSE)


##### save Rdata objects for downstream use #####

## these objects have had bad libraries and outliers removed
## the data are not normalized, and include all treatment groups

## generate simply-named combined objects, with all libraries passing qc and outlier cuts
# generate versions with and without outliers for ease of running analyses both ways

# the regular versions exclude the three outliers (3474, 4177, 4183; or preT1D082, preT1D005, T1D317)
counts.final <- counts.merged.aggregated.corrected.qc.ex_dups.ex_outliers
metrics.final <- metrics.merged.corrected.qc.ex_dups.ex_outliers
sample_annotation.final <- sample_annotation.full.corrected.qc.ex_dups.ex_outliers

# the inc_outliers versions include the three outliers (3474, 4177, 4183; or preT1D082, preT1D005, T1D317)
counts.final.inc_outliers <- counts.merged.aggregated.corrected.qc.ex_dups
metrics.final.inc_outliers <- metrics.merged.corrected.qc.ex_dups
sample_annotation.final.inc_outliers <- sample_annotation.full.corrected.qc.ex_dups

## save objects for downstream analyses
save(file="T1D_neutrophils_data_for_analysis.Rdata",
     list=c(ls_grep("\\.final"), "pcaAll.full.corrected.qc.ex_dups.ex_outliers"))
