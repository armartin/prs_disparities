library(plyr)
library(tidyverse)
library(data.table)
library(rcompanion)
library(RColorBrewer)

setwd('/Users/alicia/daly_lab/manuscripts/knowles_ashley_response/analysis/ukbb_pops')
# get population labels from here
date(); ukbb_pca <- read.table('/Users/alicia/daly_lab/skin_pigmentation/pca/ukbb_pca_pops_rf.txt.gz', header=T) %>%
  select(s, pop); date() 
# get PCs and unrelated indicator from here
date(); pca_unrel <- read.table('../ukbb_covariates_all.txt.bgz', header=T) %>%
  filter(used_in_pca_calculation=='true') %>%
  left_join(ukbb_pca, by='s'); date()

true_phenos <- read.table('../UKB_phenos_ALL17.txt.bgz', header=T) %>%
  gather('phenotype', 'true_pheno', height:plt)

covariates = c("isFemale", "age", "age_squared", "age_isFemale", "age_squared_isFemale", paste0("PC", seq(20)))

# Bootstrap 95% CI around r^2
bootstrap_CIs <- function(my_subset) {
  sample_subset <- sample_n(my_subset, nrow(my_subset), replace=TRUE)
  model1 <- lm(as.formula(sprintf("true_pheno ~ PRS + %s", paste(covariates, collapse = " + "))), data=sample_subset)
  model0 <- lm(as.formula(sprintf("true_pheno ~ %s", paste(covariates, collapse = " + "))), data=sample_subset)
  sample_r2 = summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  return(sample_r2)
}

compute_r2 <- function(my_subset) {
  model1 <- lm(as.formula(sprintf("true_pheno ~ PRS + %s", paste(covariates, collapse = " + "))), data=my_subset)
  model0 <- lm(as.formula(sprintf("true_pheno ~ %s", paste(covariates, collapse = " + "))), data=my_subset)
  my_r2 <- nagelkerke(model1, null=model0)
  r2 <- summary(model1)$adj.r.squared - summary(model0)$adj.r.squared
  date(); bs_means <- replicate(100, bootstrap_CIs(my_subset)); date()
  delta <- bs_means - r2
  d <- quantile(delta, c(0.025, 0.975))
  ci <- r2 - c(d[1], d[2])
  df <- data.frame(r2=r2,
                   ci95_low=ci[2],
                   ci95_high=ci[1],
                   p=my_r2$Likelihood.ratio.test[4],
                   pheno=my_subset$phenotype[1])
  return(df)
}

pheno_r2 <- function(pheno) {
  pheno_holdout_inds <- read.table(paste0('../ukbb/UKB_', pheno, '_PRS.txt.bgz'), header=T) %>% select(s)
  prs_inds <- pca_unrel %>%
    filter((pop=='EUR' & s %in% pheno_holdout_inds$s) | pop!='EUR')
  prs_pheno <- fread(input = paste0('zcat < UKB_', pheno, '_PRS.txt.bgz')) %>%
    select(s, s1:s10)
  
  prs_pheno_pop <- prs_pheno %>% 
    right_join(prs_inds, by='s') %>%
    left_join(true_phenos %>% subset(phenotype==pheno) %>% 
                select(eid, true_pheno, phenotype), by=c('s'='eid')) %>%
    subset(!is.na(true_pheno)) %>%
    gather('p_threshold', 'PRS', num_range('s', 1:10))
  
  r2_pheno <- prs_pheno_pop %>%
    group_by(p_threshold, pop, phenotype) %>%
    do(compute_r2(.))
  return(r2_pheno)
}

phenos <- c('height', 'bmi', 'sbp', 'dbp', 'wbc', 'monocyte', 'neutrophil', 'eosinophil', 'basophil', 'lymphocyte',
            'rbc', 'mch', 'mcv', 'mchc', 'hb', 'ht', 'plt')

date(); r2_results <- ldply(phenos, pheno_r2); date()

write.table(r2_results, 'ukbb_cross_pop_r2.txt', sep='\t', quote=F, row.names=F)
r2_results <- read.table('ukbb_cross_pop_r2.txt', header=T)
# r2_results <- read.table('ukbb_cross_pop_r2.txt', header=T)

r2_summarize <- r2_results %>%
  group_by(pop, phenotype) %>%
  dplyr::mutate(best_pred = max(r2)) %>%
  group_by(phenotype) %>%
  dplyr::mutate(best_eur_pred=max(best_pred * (pop=='EUR'))) %>%
  group_by(pop, phenotype) %>%
  dplyr::mutate(rel_eur = best_pred / best_eur_pred) %>%
  subset(r2==best_pred) 

r2_summarize_subset <- r2_summarize %>%
  subset(pop!='oth')
pheno_order <- r2_summarize_subset %>%
  subset(pop=='EUR') %>%
  arrange(r2)
r2_summarize_subset$phenotype <- factor(r2_summarize_subset$phenotype, levels=pheno_order$phenotype)

color_vec <- c(brewer.pal(8, 'Set1'))
label_vec <- c('European', 'East Asian', 'South Asian', 'African', 'American')
names(color_vec) <- names(label_vec) <- c('EUR', 'EAS', 'SAS', 'AFR', 'AMR')

pd <- position_dodge(0.5) # move them .05 to the left and right
p_r2_pop <- ggplot(r2_summarize_subset, aes(x=phenotype, y=r2, color=pop)) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=ci95_low, ymax=ci95_high), width=0, position=pd) +
  scale_color_manual(values=color_vec, name='Population', labels=label_vec) +
  labs(x='Phenotype', y=expression(R^2)) +
  theme_classic() +
  theme(axis.text = element_text(color='black', size=14),
        text = element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('r2_pheno_pop.pdf', p_r2_pop, width=8, height=6)

pop_r2_summarize <- r2_summarize %>%
  ungroup %>%
  group_by(pop) %>%
  dplyr::mutate(mean_rel_eur=mean(rel_eur), sem_rel_eur=sd(rel_eur)/sqrt(17)) %>%
  subset(phenotype=='height'& pop!='oth') %>%
  arrange(desc(mean_rel_eur))

pop_r2_summarize$pop <- factor(pop_r2_summarize$pop, levels=pop_r2_summarize$pop)
 
r2_sub <- r2_summarize %>%
  subset(pop!='oth') %>%
  arrange(desc(rel_eur))
r2_sub$pop <- factor(r2_sub$pop, levels=c('EUR', 'AMR', 'SAS', 'EAS', 'AFR'))

p_rel <- ggplot(r2_sub, aes(x=pop, y=rel_eur, fill=pop, color=pop)) +
  geom_violin() +
  geom_point(data=pop_r2_summarize, aes(x=pop, y=mean_rel_eur), color='black') +
  geom_errorbar(data=pop_r2_summarize, aes(ymin=mean_rel_eur - sem_rel_eur, ymax=mean_rel_eur + sem_rel_eur),color='black', width=0.1) +
  scale_fill_manual(values=color_vec, name='Population', labels=label_vec) +
  scale_color_manual(values=color_vec, name='Population', labels=label_vec) +
  ylim(c(0,max(r2_sub$rel_eur))) +
  scale_x_discrete(labels=label_vec) +
  labs(x='Population', y='Prediction accuracy\n(relative to Europeans)') +
  theme_classic() +
  guides(fill=F, color=F) +
  theme(axis.text = element_text(color='black', size=16),
        text = element_text(size=18),
        axis.text.x = element_text(angle = 45, hjust = 1))

p_rel <- ggplot(pop_r2_summarize, aes(x=pop, y=mean_rel_eur, fill=pop, color=pop)) +
  #geom_bar(stat='identity') +
  geom_point() +
  geom_errorbar(aes(ymin=mean_rel_eur - sem_rel_eur, ymax=mean_rel_eur + sem_rel_eur), width=0.2) +
  scale_color_manual(values=color_vec, name='Population', labels=label_vec) +
  scale_x_discrete(labels=label_vec) +
  labs(x='Population', y='Prediction accuracy\n(relative to Europeans)') +
  theme_classic() +
  guides(fill=F) +
  theme(axis.text = element_text(color='black', size=16),
        text = element_text(size=18),
        axis.text.x = element_text(angle = 45, hjust = 1))

save(p_rel, file='p_rel.RData')
ggsave('prs_r2_rel_eur2.pdf', p_rel, width=8, height=5)
