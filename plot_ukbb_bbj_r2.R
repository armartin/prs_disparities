library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(plyr)
library(tidyr)
library(dplyr)

####################################################
### Read in data to make plots

library(ggplot2)

read_r2 <- function(filename, target, sumstats) {
  # nomenclature is [target cohort]_[sumstats cohort]
  target_sumstats <- read.table(filename, header=T)
  colnames(target_sumstats) <- c('r2', 'CI95_low', 'CI95_high', 'p', 'pheno', 'p_threshold')
  target_sumstats$pheno <- revalue(target_sumstats$pheno, c("basophil"="Basophil", "bmi"="BMI", 'dbp'='DBP', 'eosinophil'='Eosinophil', 
                                                            'hb'='HB', 'height'='Height', 'ht'='HT', 'lymphocyte'='Lymphocyte', 
                                                            'mch'='MCH', 'mchc'='MCHC','mcv'='MCV', 'monocyte'='Monocyte', 
                                                            'neutrophil'='Neutrophil', 'plt'='Platelet', 'rbc'='RBC', 'sbp'='SBP', 'wbc'='WBC'))
  target_sumstats$p_threshold <- factor(target_sumstats$p_threshold, levels = paste0('s', 1:10))
  target_sumstats$p_threshold <- revalue(target_sumstats$p_threshold, c('s1'='5e-8', 's2'='1e-6', 's3'='1e-4', 's4'='1e-3', 
                                                                        's5'='1e-2', 's6'='.05', 's7'='.1', 's8'='.2', 
                                                                        's9'='.5', 's10'='1'))
  target_sumstats$target <- target
  target_sumstats$sumstats <- sumstats
  return(target_sumstats)
}

bbj_bbj <- read_r2('bbj_from_bbj.r2.txt', 'BBJ', 'BBJ')
bbj_ukbb <- read_r2('bbj_from_ukbb.r2.txt', 'BBJ', 'UKBB')
ukbb_ukbb <- read_r2('ukbb_from_ukbb.r2.txt', 'UKBB', 'UKBB')
ukbb_bbj <- read_r2('ukbb_from_bbj.r2.txt', 'UKBB', 'BBJ')

####################################################
### Make supplementary plots for all p-value thresholds

bbj <- rbind(bbj_bbj, bbj_ukbb)
ukbb <- rbind(ukbb_ukbb, ukbb_bbj)

color_vec <- brewer.pal(3, 'Set1')
names(color_vec) <- c('UKBB', 'BBJ', '')

plot_all_r2 <- function(dataset, data_name) {
  pd <- position_dodge(0.2) # move them .05 to the left and right
  ggplot(dataset, aes(x=p_threshold, y=r2, color=sumstats)) +
    geom_point(position=pd) +
    #ylim(-0.01,0.17) +
    geom_errorbar(aes(ymin=CI95_low, ymax=CI95_high), width=.1, position=pd) +
    facet_wrap(~pheno) +
    scale_color_manual(values=color_vec, name='Prediction\n(sumstats)') +
    labs(x='P-value threshold', y=bquote(R^2 ~'(in'~ .(data_name) * ')')) +
           #bquote(R^2 ~ (in) ~.(data_name) ~ ())) +
           #expression(R^2~'(in ', data_name, ')')) +
    theme_classic() +
    theme(strip.background = element_rect(fill = "lightgrey"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text = element_text(color='black'),
          text = element_text(size=14))
}


p_bbj <- plot_all_r2(bbj, 'BBJ')
p_ukbb <- plot_all_r2(ukbb, 'UKBB')

p_bbj_ukbb <- plot_grid(p_bbj, p_ukbb, labels=c('A', 'B'), ncol=1)
save_plot('bbj_ukbb_r2_all_p.pdf', p_bbj_ukbb, nrow=2, base_height=6, base_width=10)

####################################################
### Make main figures

all_cohorts <- rbind(bbj_bbj, bbj_ukbb, ukbb_ukbb, ukbb_bbj)
write.table(all_cohorts, 'all_cohorts.txt', quote=F, row.names=F, sep='\t')

all_filt <- all_cohorts %>%
  group_by(target, sumstats, pheno) %>%
  arrange(desc(r2)) %>%
  slice(1)

bbj_filt <- subset(all_filt, target=='BBJ') %>% arrange(r2)
ukbb_filt <- subset(all_filt, target=='UKBB') %>% arrange(r2)
bbj_bbj_filt <- subset(bbj_filt, sumstats=='BBJ')
ukbb_ukbb_filt <- subset(ukbb_filt, sumstats=='UKBB')
bbj_filt$pheno <- factor(bbj_filt$pheno, levels=as.character(bbj_bbj_filt$pheno))
ukbb_filt$pheno <- factor(ukbb_filt$pheno, levels=as.character(ukbb_ukbb_filt$pheno))

plot_top_r2 <- function(dataset, data_name, legend=FALSE) {
  pd <- position_dodge(0.2) # move them .05 to the left and right
  p <- ggplot(dataset, aes(x=pheno, y=r2, color=sumstats)) +
    geom_point(position=pd) +
    geom_errorbar(aes(ymin=CI95_low, ymax=CI95_high), width=.1, position=pd) +
    scale_color_manual(values=color_vec, name='GWAS\ncohort\n(sumstats)') +
    labs(x='Phenotype', y=bquote(R^2 ~'(in'~ .(data_name) * ')')) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          text = element_text(size=16, color='black'))
  if(!legend) {
    p <- p + guides(color=F)
  }
  return(p)
}

p_bbj_top <- plot_top_r2(bbj_filt, 'BBJ')
p_ukbb_top <- plot_top_r2(ukbb_filt, 'UKBB', T)

p_bbj_ukbb_top <- plot_grid(p_bbj_top, p_ukbb_top, labels=c('A', 'B'), rel_widths = c(1, 1.2))
save_plot('bbj_ukbb_r2_top.pdf', p_bbj_ukbb_top, base_height=6, base_width=12)

