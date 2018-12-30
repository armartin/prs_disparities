library(tidyverse)
library(RColorBrewer)
library(cowplot)

setwd('/Users/alicia/daly_lab/manuscripts/knowles_ashley_response/analysis')
####################################################
### Read in data to make plots

read_r2 <- function(filename, target, sumstats, cc=FALSE) {
  # nomenclature is [target cohort]_[sumstats cohort]
  target_sumstats <- read.table(filename, header=T)
  if (cc) {
    # observed: 1-3, liability: 4-6
    target_sumstats <- target_sumstats[,c(4:ncol(target_sumstats))]
  }
  colnames(target_sumstats) <- c('r2', 'CI95_low', 'CI95_high', 'p', 'pheno', 'p_threshold')
  if (cc) {
    target_sumstats$pheno <- revalue(target_sumstats$pheno, c("crc" = "CRC", "glaucoma" = "Glaucoma", "afib" = "AFib", "t2d" = "T2D", "ra" = "RA"))
  } else {
    target_sumstats$pheno <- revalue(target_sumstats$pheno, c("basophil"="Basophil", "bmi"="BMI", 'dbp'='DBP', 'eosinophil'='Eosinophil',
                                                              'hb'='HB', 'height'='Height', 'ht'='HT', 'lymphocyte'='Lymphocyte',
                                                              'mch'='MCH', 'mchc'='MCHC','mcv'='MCV', 'monocyte'='Monocyte',
                                                              'neutrophil'='Neutrophil', 'plt'='Platelet', 'rbc'='RBC', 'sbp'='SBP', 'wbc'='WBC'))
  }

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
afr_ukbb <- read_r2('afr_from_ukbb.r2.txt', 'UKBB AFR', 'UKBB')
afr_bbj <- read_r2('afr_from_bbj.r2.txt', 'UKBB AFR', 'BBJ')

bbj_bbj_cc <- read_r2('bbj_from_bbj_cc.r2.txt', 'BBJ', 'BBJ', cc=TRUE)
bbj_ukbb_cc <- read_r2('bbj_from_ukbb_cc.r2.txt', 'BBJ', 'UKBB', cc=TRUE)
ukbb_ukbb_cc <- read_r2('ukbb_from_ukbb_cc.r2.txt', 'UKBB', 'UKBB', cc=TRUE)
ukbb_bbj_cc <- read_r2('ukbb_from_bbj_cc.r2.txt', 'UKBB', 'BBJ', cc=TRUE)

####################################################
### Make supplementary plots for all p-value thresholds

bbj <- rbind(bbj_bbj, bbj_ukbb)
ukbb <- rbind(ukbb_ukbb, ukbb_bbj)
afr <- rbind(afr_ukbb, afr_bbj)

bbj_cc = rbind(bbj_bbj_cc, bbj_ukbb_cc)
ukbb_cc = rbind(ukbb_ukbb_cc, ukbb_bbj_cc)

color_vec <- brewer.pal(3, 'Set1')
names(color_vec) <- c('UKBB', 'BBJ', 'UKBB AFR')

plot_all_r2 <- function(dataset, data_name, cc=FALSE) {
  pd <- position_dodge(0.2) # move them .05 to the left and right
  if (cc) {
    ylab = bquote(R[Liability]^2 ~'(in'~ .(data_name) * ')')
  } else {
    ylab = bquote(R^2 ~'(in'~ .(data_name) * ')')
  }
  ggplot(dataset, aes(x=p_threshold, y=r2, color=sumstats)) +
    geom_point(position=pd) +
    #ylim(-0.01,0.17) +
    geom_errorbar(aes(ymin=CI95_low, ymax=CI95_high), width=.1, position=pd) +
    facet_wrap(~pheno) +
    scale_color_manual(values=color_vec, name='Prediction\n(sumstats)') +
    labs(x='P-value threshold', y=ylab) +
    theme_classic() +
    theme(strip.background = element_rect(fill = "lightgrey"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),#, size=8),
          axis.text = element_text(color='black'),
          text = element_text(size=14))
}


p_bbj <- plot_all_r2(bbj, 'BBJ')
p_ukbb <- plot_all_r2(ukbb, 'UKBB European descent')
p_afr <- plot_all_r2(afr, 'UKBB African descent')

p_bbj_cc <- plot_all_r2(bbj_cc, 'BBJ', cc=TRUE)
p_ukbb_cc <- plot_all_r2(ukbb_cc, 'UKBB European descent', cc=TRUE)

save_plot('bbj_r2_all_p.pdf', p_bbj, base_width=8, base_height=5)
save_plot('ukbb_r2_all_p.pdf', p_ukbb, base_width=8, base_height=5)
save_plot('afr_r2_all_p.pdf', p_afr, base_width=8, base_height=5)

save_plot('bbj_cc_r2_all_p.pdf', p_bbj_cc, base_width=8, base_height=5)
save_plot('ukbb_cc_r2_all_p.pdf', p_ukbb_cc, base_width=8, base_height=5)

p_bbj_ukbb <- plot_grid(p_bbj, p_ukbb, p_afr, labels=c('A', 'B', 'C'), ncol=2)
save_plot('bbj_ukbb_r2_all_p.pdf', p_bbj_ukbb, nrow=2, base_height=4, base_width=15)

p_bbj_ukbb <- plot_grid(p_bbj_cc, p_ukbb_cc, labels=c('A', 'B'), ncol=2)
save_plot('bbj_ukbb_cc_r2_all_p.pdf', p_bbj_ukbb, nrow=2, base_height=4, base_width=15)


####################################################
### Make main figures

all_cohorts <- rbind(bbj_bbj, bbj_ukbb, ukbb_ukbb, ukbb_bbj, afr_ukbb, afr_bbj)
write.table(all_cohorts, 'all_cohorts.txt', quote=F, row.names=F, sep='\t')

all_filt <- all_cohorts %>%
  group_by(target, sumstats, pheno) %>%
  arrange(desc(r2)) %>%
  slice(1)

bbj_filt <- subset(all_filt, target=='BBJ') %>% arrange(r2)
ukbb_filt <- subset(all_filt, target=='UKBB') %>% arrange(r2)
afr_filt <- subset(all_filt, target=='UKBB AFR') %>% arrange(r2)
bbj_bbj_filt <- subset(bbj_filt, sumstats=='BBJ')
ukbb_ukbb_filt <- subset(ukbb_filt, sumstats=='UKBB')
afr_ukbb_filt <- subset(afr_filt, sumstats=='UKBB')
bbj_filt$pheno <- factor(bbj_filt$pheno, levels=as.character(bbj_bbj_filt$pheno))
ukbb_filt$pheno <- factor(ukbb_filt$pheno, levels=as.character(ukbb_ukbb_filt$pheno))
afr_filt$pheno <- factor(afr_filt$pheno, levels=as.character(afr_ukbb_filt$pheno))

all_filt$pheno <- factor(all_filt$pheno, levels=as.character(ukbb_ukbb_filt$pheno))
ukbb_sumstats <- subset(all_filt, sumstats=='UKBB')
all_filt$pheno <- factor(all_filt$pheno, levels=as.character(bbj_bbj_filt$pheno))
bbj_sumstats <- subset(all_filt, sumstats=='BBJ')

# disease endpoints
all_cohorts_cc <- rbind(bbj_bbj_cc, bbj_ukbb_cc, ukbb_ukbb_cc, ukbb_bbj_cc)
write.table(all_cohorts_cc, 'all_cohorts_cc.txt', quote=F, row.names=F, sep='\t')

all_filt_cc <- all_cohorts_cc %>%
  group_by(target, sumstats, pheno) %>%
  arrange(desc(r2)) %>%
  slice(1)

bbj_filt_cc <- subset(all_filt_cc, target=='BBJ') %>% arrange(r2)
ukbb_filt_cc <- subset(all_filt_cc, target=='UKBB') %>% arrange(r2)
bbj_bbj_filt_cc <- subset(bbj_filt_cc, sumstats=='BBJ')
ukbb_ukbb_filt_cc <- subset(ukbb_filt_cc, sumstats=='UKBB')
bbj_filt_cc$pheno <- factor(bbj_filt_cc$pheno, levels=as.character(bbj_bbj_filt_cc$pheno))
ukbb_filt_cc$pheno <- factor(ukbb_filt_cc$pheno, levels=as.character(ukbb_ukbb_filt_cc$pheno))

all_filt_cc$pheno <- factor(all_filt_cc$pheno, levels=as.character(ukbb_ukbb_filt_cc$pheno))
ukbb_sumstats_cc <- subset(all_filt_cc, sumstats=='UKBB')
all_filt_cc$pheno <- factor(all_filt_cc$pheno, levels=as.character(bbj_bbj_filt_cc$pheno))
bbj_sumstats_cc <- subset(all_filt_cc, sumstats=='BBJ')

pd <- position_dodge(0.2) # move them .05 to the left and right
p1 <- ggplot(ukbb_sumstats, aes(x=pheno, y=r2, color=target)) +
  geom_point(position=pd) +
  #facet_grid(~sumstats) +
  geom_errorbar(aes(ymin=CI95_low, ymax=CI95_high), width=.1, position=pd) +
  scale_color_manual(values=color_vec, name='Target cohort') +
  labs(x='Phenotype', y=bquote(R^2), title='UKBB sumstats') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=16, color='black'))

p2 <- ggplot(bbj_sumstats, aes(x=pheno, y=r2, color=target)) +
  geom_point(position=pd) +
  #facet_grid(~sumstats) +
  geom_errorbar(aes(ymin=CI95_low, ymax=CI95_high), width=.1, position=pd) +
  scale_color_manual(values=color_vec, name='Target cohort') +
  labs(x='Phenotype', y=bquote(R^2), title='BBJ sumstats') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=16, color='black'))

p_bbj_ukbb_test <- plot_grid(p1, p2, labels=c('A', 'B'))#, hjust=c(-0.5,-0.5,-0.5,1))

ggsave('bbj_ukbb_r2_test.pdf', p_bbj_ukbb_test, width=12, height=6)

plot_top_r2 <- function(dataset, data_name, legend=FALSE, cc=FALSE) {
  pd <- position_dodge(0.2) # move them .05 to the left and right
  if (cc) {
    ylab = bquote(R[Liability]^2 ~'(in'~ .(data_name) * ')')
  } else {
    ylab = bquote(R^2 ~'(in'~ .(data_name) * ')')
  }
  p <- ggplot(dataset, aes(x=pheno, y=r2, color=sumstats)) +
    geom_point(position=pd) +
    geom_errorbar(aes(ymin=CI95_low, ymax=CI95_high), width=.1, position=pd) +
    scale_color_manual(values=color_vec, name='GWAS\ncohort\n(sumstats)') +
    labs(x='Phenotype', y=ylab) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          text = element_text(size=16, color='black'))
  if(!legend) {
    p <- p + guides(color=F)
  }
  return(p)
}

p1 <- plot_top_r2(bbj_filt, 'BBJ', legend=TRUE)
p_bbj_top <- plot_top_r2(bbj_filt, 'BBJ')
p_ukbb_top <- plot_top_r2(ukbb_filt, 'UKBB')
p_afr_top <- plot_top_r2(afr_filt, 'UKBB African descent')

legend <- get_legend(p1)

p_bbj_ukbb_top <- plot_grid(p_ukbb_top, p_bbj_top, p_afr_top, legend, labels=c('A', 'B', 'C', ''))#, hjust=c(-0.5,-0.5,-0.5,1))
save_plot('bbj_ukbb_r2_top.pdf', p_bbj_ukbb_top, base_height=10, base_width=10)


# Plot h2 -----------------------------------------------------------------

h2 <- read.table('baselineLD.txt', header=T, sep='\t') %>%
  mutate(se_low = h2 - SE, se_high = h2 + SE)
trait_order <- (h2 %>% filter(population=='UKBB'&model=='S-LDSC (baseline-LD)') %>% arrange(desc(h2)))$Trait
h2$Trait <- factor(h2$Trait, levels=trait_order)
pd <- position_dodge(0.2) # move them .05 to the left and right

p_h2 <- ggplot(h2, aes(x=Trait, y=h2, color=population)) +
  geom_point(position=pd) +
  facet_wrap(~model) +
  geom_errorbar(aes(ymin=se_low, ymax=se_high), width=.1, position=pd) +
  scale_color_manual(values=color_vec, name='Population') +
  scale_shape(name='LDSC model') +
  labs(y=expression(h^2)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size=16, color='black'))

save_plot('h2_bbj_ukbb.pdf', p_h2, base_height=5,base_width=10)

p_bbj_cc_top <- plot_top_r2(bbj_filt_cc, 'BBJ', cc=TRUE)
p_ukbb_cc_top <- plot_top_r2(ukbb_filt_cc, 'UKBB', cc=TRUE)

legend <- get_legend(p1)

p_bbj_ukbb_top <- plot_grid(p_bbj_top, p_ukbb_top, p_afr_top, p_bbj_cc_top, p_ukbb_cc_top, legend, labels=c('A', 'B', 'C', 'D', 'E', ''))#, hjust=c(-0.5,-0.5,-0.5,1))
save_plot('bbj_ukbb_r2_top.pdf', p_bbj_ukbb_top, base_height=10, base_width=15)
