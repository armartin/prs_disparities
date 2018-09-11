library(plyr)
library(dplyr)
library(tidyr)
library(doParallel)

cl = makeCluster(8, 'PSOCK')
registerDoParallel(cl)

afr=FALSE

if(afr) {
  setwd('/Users/alicia/daly_lab/manuscripts/knowles_ashley_response/analysis/ukbb_afr/')
} else {
  setwd('/Users/alicia/daly_lab/manuscripts/knowles_ashley_response/analysis/')
}

covariates = c("isFemale", "age", "age_squared", "age_isFemale", "age_squared_isFemale", paste0("PC", seq(20)))

# Bootstrap 95% CI around r^2
bootstrap_CIs <- function(my_subset) {
  sample_subset <- sample_n(my_subset, nrow(my_subset), replace=TRUE)
  model1 <- lm(as.formula(sprintf("true_pheno ~ PRS + %s", paste(covariates, collapse = " + "))), data=sample_subset)
  model0 <- lm(as.formula(sprintf("true_pheno ~ %s", paste(covariates, collapse = " + "))), data=sample_subset)
  sample_r2 = summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  return(sample_r2)
}

# Get R^2 and p-value for nested model
compute_r2 <- function(dataset, phenotype, p) {
  print(paste(phenotype, p))
  my_subset <- subset(dataset, p_threshold==p & pheno_name==phenotype)
  model1 <- lm(as.formula(sprintf("true_pheno ~ PRS + %s", paste(covariates, collapse = " + "))), data=my_subset)
  model0 <- lm(as.formula(sprintf("true_pheno ~ %s", paste(covariates, collapse = " + "))), data=my_subset)
  r2 <- summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  pval <- anova(model1, model0)$'Pr(>F'[2]
  date(); bs_means <- replicate(1000, bootstrap_CIs(my_subset)); date()
  delta <- bs_means - r2
  d = quantile(delta, c(0.025, 0.975))
  ci = r2 - c(d[1], d[2])
  r2_vec <- c(r2, rev(ci), pval, my_subset$pheno_name[1], my_subset$p_threshold[1])
  names(r2_vec) <- c('r2', '2.5%', '97.5%', 'p', 'pheno', 'p_threshold')
  return(r2_vec)
}

# Read in and wrangle PRS + phenotype + covariate data
read_phenos <- function(pheno_name, afr) {
  if(afr) {
    my_pheno <- read.table(paste0('ukbb/UKB_', pheno_name, '_afr_PRS.txt.bgz'), header=T)  
    afr_subset <- read.table('../ukb31063_afr_subset_samples.tsv', header=T)
    afr_subset <- subset(afr_subset, afr_select)
    my_pheno <- subset(my_pheno, s %in% afr_subset$s)
  } else {
    my_pheno <- read.table(paste0('ukbb/UKB_', pheno_name, '_PRS.txt.bgz'), header=T)
  }
  my_pheno <- my_pheno %>%
    select(s:PC20,pheno_name,s1:s5)
  colnames(my_pheno)[27] <- 'true_pheno'
  my_pheno <- my_pheno %>%
    gather('p_threshold', 'PRS', num_range('s', 1:10))
  my_pheno$pheno_name <- pheno_name
  # my_pheno$pheno <- pheno_name
  return(my_pheno)
}

phenos <- c('height', 'bmi', 'sbp', 'dbp', 'wbc', 'monocyte', 'neutrophil', 'eosinophil', 'basophil', 'lymphocyte',
            'rbc', 'mch', 'mcv', 'mchc', 'hb', 'ht', 'plt')


date(); all_phenos <- ldply(phenos, function(x) { ldply(afr, function(y) {read_phenos(x, y) } ) }); date()

p_thresholds <- paste0('s', 1:5)

# run this for all phenos, takes ~3 minutes per phenotype
date(); r2_combined <- ldply(phenos, function(x) { ldply(p_thresholds, function(y) { compute_r2(all_phenos, x, y) } ) }, 
                             .parallel=T, .paropts=list(.export = c("compute_r2", "bootstrap_CIs", "p_thresholds", "all_phenos", "covariates"),
                                                        .packages = c('dplyr', 'plyr'))); date()

write.table(r2_combined, 'ukbb_from_ukbb.r2.txt', quote = F, row.names = F, sep = '\t')

stopCluster(cl)
