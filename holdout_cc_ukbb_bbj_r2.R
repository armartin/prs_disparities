library(plyr)
library(dplyr)
library(tidyr)
library(doParallel)
library(rcompanion)

cl = makeCluster(8, 'PSOCK')
registerDoParallel(cl)

covariates = c("isFemale", "age", "age_squared", "age_isFemale", "age_squared_isFemale", paste0("PC", seq(20)))

prevalence = list(crc = 0.02, t2d = 0.05, glaucoma = 0.02, afib = 0.02, ra = 0.01)

# reference: https://github.com/armartin/pgc_scz_asia/blob/master/eur_eas_prs.R
h2l_R2N <- function(k, r2n, p) {
  # k = population prevalence
  # r2n = Nagelkerke's attributable to genomic profile risk score
  # p = proportion of cases in the case-control sample
  # calculates proportion of variance explained on the liability scale
  # from ABC at http://www.complextraitgenomics.com/software/
  # Lee SH, Goddard ME, Wray NR, Visscher PM. (2012) A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012 Apr;36(3):214-24.
  x <- qnorm(1 - k)
  z <- dnorm(x)
  i <- z / k
  cc <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
  theta <- i * ((p - k)/(1 - k)) * (i * ((p - k) / ( 1 - k)) - x)
  e <- 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
  h2l_R2N <- cc * e * r2n / (1 + cc * e * theta * r2n)
}

# Bootstrap 95% CI around r^2
bootstrap_CIs <- function(my_subset, phenotype) {
  sample_subset <- sample_n(my_subset, nrow(my_subset), replace=TRUE)
  model1 <- glm(as.formula(sprintf("true_pheno ~ PRS + %s", paste(covariates, collapse = " + "))), data=sample_subset, family = 'binomial')
  model0 <- glm(as.formula(sprintf("true_pheno ~ %s", paste(covariates, collapse = " + "))), data=sample_subset, family = 'binomial')
  my_r2 <- nagelkerke(model1, null=model0)
  r2n <- my_r2$Pseudo.R.squared.for.model.vs.null[3]
  r2l <- h2l_R2N(k=prevalence[[phenotype]], r2n=r2n, p=sum(my_subset$true_pheno==1)/sum(my_subset$true_pheno %in% c(0,1)))
  return(data.frame(r2n = r2n, r2l = r2l))
}

# Get R^2 and p-value for nested model
compute_r2 <- function(dataset, phenotype, p) {
  print(paste(phenotype, p))
  my_subset <- subset(dataset, p_threshold==p & pheno_name==phenotype)
  model1 <- glm(as.formula(sprintf("true_pheno ~ PRS + %s", paste(covariates, collapse = " + "))), data=my_subset, family = 'binomial')
  model0 <- glm(as.formula(sprintf("true_pheno ~ %s", paste(covariates, collapse = " + "))), data=my_subset, family = 'binomial')
  my_r2 <- nagelkerke(model1, null=model0)
  r2n <- my_r2$Pseudo.R.squared.for.model.vs.null[3]
  r2l <- h2l_R2N(k=prevalence[[phenotype]], r2n=r2n, p=sum(my_subset$true_pheno==1)/sum(my_subset$true_pheno %in% c(0,1)))
  pval <- my_r2$Likelihood.ratio.test[4]
  bs_means <- ldply(1:1000, function(i){bootstrap_CIs(my_subset, phenotype)})
  delta_r2n <- bs_means$r2n - r2n
  d_r2n = quantile(delta_r2n, c(0.025, 0.975))
  ci_r2n = r2n - c(d_r2n[1], d_r2n[2])
  delta_r2l <- bs_means$r2l - r2l
  d_r2l = quantile(delta_r2l, c(0.025, 0.975))
  ci_r2l = r2l - c(d_r2l[1], d_r2l[2])
  r2_vec <- c(r2n, rev(ci_r2n), r2l, rev(ci_r2l), pval, my_subset$pheno_name[1], my_subset$p_threshold[1])
  names(r2_vec) <- c('r2n', 'r2n_2.5%', 'r2n_97.5%', 'r2l', 'r2l_2.5%', 'r2l_97.5%', 'p', 'pheno', 'p_threshold')
  return(r2_vec)
}


# Read in and wrangle PRS + phenotype + covariate data
read_phenos <- function(pheno_name) {
  my_pheno <- read.table(paste0('bbj/UKB_', pheno_name, '_PRS.txt.bgz'), header=T)
  my_pheno <- my_pheno %>%
    select(s:PC20,pheno_name,s1:s5)
  colnames(my_pheno)[27] <- 'true_pheno'
  my_pheno <- my_pheno %>%
    gather('p_threshold', 'PRS', num_range('s', 1:10))
  my_pheno$pheno_name <- pheno_name
  # my_pheno$pheno <- pheno_name
  return(my_pheno)
}

phenos <- c('crc', 't2d', 'glaucoma', 'afib', 'ra')


date(); all_phenos <- ldply(phenos, function(x) { read_phenos(x) }); date()

p_thresholds <- paste0('s', 1:5)

# run this for all phenos, takes ~20 minutes for all phenotypes
date(); r2_combined <- ldply(phenos, function(x) { ldply(p_thresholds, function(y) { compute_r2(all_phenos, x, y) } ) },
                             .parallel=T,
                             .paropts=list(.export = c("h2l_R2N", "all_phenos", "p_thresholds", "covariates", "prevalence", "compute_r2", "bootstrap_CIs"),
                             .packages = c('dplyr', 'plyr', 'rcompanion'))); date()

write.table(r2_combined, 'ukbb_from_bbj_cc.r2.txt', quote = F, row.names = F, sep = '\t')

stopCluster(cl)

