library(ggplot2)
library(dplyr)
library(RColorBrewer)

setwd('/Users/alicia/daly_lab/manuscripts/knowles_ashley_response/analysis')

read_pops <- function(filename, popname) {
  pop <- read.table(filename, header=T)
  pop$pop <- popname
  return(pop)
}

eur <- read_pops('EUR_gwas_catalog_v1.0.2-associations_e93_r2018-08-14.frq', 'EUR')
eas <- read_pops('EAS_gwas_catalog_v1.0.2-associations_e93_r2018-08-14.frq', 'EAS')
afr <- read_pops('AFR_gwas_catalog_v1.0.2-associations_e93_r2018-08-14.frq', 'AFR')

all_pops <- rbind(eur, eas, afr)
all_pops$pop <- factor(all_pops$pop, levels=c('EAS', 'AFR', 'EUR'))

brewer_vec <- brewer.pal(3, 'Set1')
names(brewer_vec) <- c('EUR', 'EAS', 'AFR')

p1 <- ggplot(all_pops, aes(x=MAF, fill=pop, color=pop)) +
  geom_density(alpha=0.7) +
  scale_fill_manual(values=brewer_vec, name='Population') +
  scale_color_manual(values=brewer_vec, name='Population') +
  labs(x='Minor allele frequency', y='Density', title='GWAS catalog sites') +
  guides(color=F, fill=guide_legend(override.aes = list(alpha=1))) +
  theme_classic() +
  theme(text = element_text(size=16),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        axis.text = element_text(color='black'),
        plot.title = element_text(hjust = 0.5))

ggsave(filename='gwas_sfs.pdf', p1, width=5, height=4)
