library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)

setwd('/Users/alicia/daly_lab/manuscripts/knowles_ashley_response/analysis')

prs <- read.table('prs_r2_anc.txt2', header=T, sep='\t')
prs_target <- subset(prs, Vary.study.vs.target.population=='Target')
prs_target$trait_ref <- paste0(prs_target$Phenotype, '\n(', prs_target$Ref_short, ')')
prs_target$x=1:nrow(prs_target)

color_vec <- brewer.pal(3, 'Set1')
names(color_vec) <- c('European', 'East Asian', 'African American')

p1 <- ggplot(prs_target, aes(x=x, y=R2, fill=Target.population)) +
  #geom_bar(stat='identity', position='dodge') +
  geom_col(position = "dodge") +
  facet_grid(~trait_ref, scales = 'free_x', space='free_x', switch='x') +
  scale_fill_manual(values = color_vec, name='Target\npopulation') +
  labs(x='Trait/study', y='Variance explained')+
  #geom_text(aes(label=round(Relative.to.European, 2)), position=position_dodge(width=0.9), vjust=0.25, hjust=-0.3, size=4, angle=90) +
  theme_classic() +
  theme(text = element_text(size=14),
        axis.text = element_text(color='black'),
        strip.text = element_text(angle = 90),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        legend.position = 'bottom')

ggsave('multiethnic_prs_empirical.pdf', p1)

rel_prs <- subset(prs_target, Target.population!='European')
rel_prs$trait_ref_cohort <- paste0(rel_prs$Phenotype, '\n(', rel_prs$Ref_short, ', ', rel_prs$Target.cohort, ')')
rel_prs$trait_ref2 <- paste0(rel_prs$Phenotype, ' (', rel_prs$Ref_short, ')')
rel_prs$shape='Study'
rel_prs <- rel_prs %>%
  group_by(Target.population) %>%
  arrange(Relative.to.European)
total <- rel_prs %>%
  group_by(Target.population) %>%
  summarize_if(is.numeric, mean, na.rm=T)
total$trait_ref_cohort <- total$trait_ref2 <- c('Mean (African American)', 'Mean (East Asian)')
total$shape='Total'
rel_prs <- bind_rows(rel_prs, total)
rel_prs$trait_ref_cohort <- factor(rel_prs$trait_ref_cohort, levels=rev(c(
  subset(rel_prs, Target.population=='East Asian')$trait_ref_cohort,
  subset(rel_prs, Target.population=='African American')$trait_ref_cohort)))

shape_vec <- c(20, 18)
names(shape_vec) <- c('Study', 'Total')
size_vec <- c(4, 8)
names(size_vec) <- c('Study', 'Total')
p2 <- ggplot(rel_prs, aes(y=trait_ref_cohort, x=Relative.to.European, color=Target.population, shape=shape, size=shape)) +
  geom_point() +
  xlim(0,1) +
  scale_color_manual(values=color_vec, name='Target population') +
  scale_shape_manual(values=shape_vec) +
  scale_size_manual(values=size_vec) +
  geom_vline(xintercept=1, color="gray20",linetype="longdash") +
  geom_hline(yintercept=9.5, color="gray60",linetype="dotted") +
  scale_y_discrete(labels=rel_prs$trait_ref2[order(rel_prs$trait_ref_cohort)]) +
  theme_bw() +
  labs(x='Proportion variance explained\nrelative to Europeans', y='Study cohort') +
  guides(shape=F, size=F, color=F) +
  theme(text = element_text(size=14),
        axis.text = element_text(color='black'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank())

p3 <- plot_grid(p1, p2, labels=c('A', 'B'))
save_plot(filename='multiethnic_prs_empirical2.pdf', p3, base_width=12, base_height=5)
