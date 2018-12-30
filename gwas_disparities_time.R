library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(zoo)
library(cowplot)

setwd('/Users/alicia/daly_lab/manuscripts/knowles_ashley_response/analysis')

gwas_header <- colnames(read.delim('gwas_catalog-ancestry_r2018-07-17.tsv', header=T, sep='\t', row.names=NULL))
gwas <- read.delim('gwas_catalog-ancestry_r2018-07-17.tsv', header=F, sep='\t', skip=1)[,1:12]
colnames(gwas) <- gwas_header[2:13]
gwas$DATE <- as.Date(gwas$DATE)

# Order GWAS catalog by date and plot cumulative number of individuals
gwas <- gwas %>% arrange(DATE)
gwas$cum_num <- gwas$NUMBER.OF.INDIVDUALS
for(i in 2:nrow(gwas)) {
  gwas$cum_num[i] <- sum(gwas$cum_num[i-1], gwas$cum_num[i], na.rm=T)
}

# Group individuals by broad ancestral categories
gwas_easy <- gwas[!grepl(', ', gwas$BROAD.ANCESTRAL.CATEGORY) & gwas$BROAD.ANCESTRAL.CATEGORY!='NR',]
mideast <- gwas[gwas$BROAD.ANCESTRAL.CATEGORY=='Greater Middle Eastern (Middle Eastern, North African or Persian)',]
africa <- gwas[gwas$BROAD.ANCESTRAL.CATEGORY %in% c('Sub-Saharan African, African American or Afro-Caribbean',
                                                    'Sub-Saharan African, African unspecified'),]
africa$BROAD.ANCESTRAL.CATEGORY <- 'African unspecified'
asia <- gwas[gwas$BROAD.ANCESTRAL.CATEGORY %in% c('East Asian, Asian unspecified', 'South Asian, East Asian ',
                                                  'South Asian, South East Asian', 'South Asian, South East Asian, East Asian',
                                                  'South East Asian, East Asian', 'South East Asian, South Asian, East Asian'),]
asia$BROAD.ANCESTRAL.CATEGORY <- 'Asian unspecified'
nr <- gwas[gwas$BROAD.ANCESTRAL.CATEGORY=='NR',]
multiple <- gwas[grepl(', ', gwas$BROAD.ANCESTRAL.CATEGORY) & gwas$BROAD.ANCESTRAL.CATEGORY!='NR' &
                   gwas$BROAD.ANCESTRAL.CATEGORY!='Greater Middle Eastern (Middle Eastern, North African or Persian)' &
                   !gwas$BROAD.ANCESTRAL.CATEGORY %in% c('Sub-Saharan African, African American or Afro-Caribbean',
                                                         'Sub-Saharan African, African unspecified') &
                   !gwas$BROAD.ANCESTRAL.CATEGORY %in% c('East Asian, Asian unspecified', 'South Asian, East Asian ',
                                                         'South Asian, South East Asian', 'South Asian, South East Asian, East Asian',
                                                         'South East Asian, East Asian', 'South East Asian, South Asian, East Asian') &
                   gwas$BROAD.ANCESTRAL.CATEGORY!='NR',]
multiple$BROAD.ANCESTRAL.CATEGORY <- 'Multiple'

gwas_simplified <- bind_rows(gwas_easy, mideast, africa, asia, nr, multiple)

anc_categories <- data.frame(ancestry=sort(unique(gwas_simplified$BROAD.ANCESTRAL.CATEGORY)), 
                             category=c(rep('Non-EURASN', 3), rep('ASN', 3), 'EUR', rep('Non-EURASN', 2), 'Multiple', 'Non-EURASN', 'NR', rep('Non-EURASN', 3), rep('ASN', 2), 'Non-EURASN'),
                             category2=c('Oceanic', rep('African', 2), rep('South Asian/Other Asian', 2), 'East Asian', 'European', 'Greater Middle Eastern', 'Hispanic/Latino', 'Multiple', 'Hispanic/Latino', 'Not Reported', 'Oceanic', rep('Other', 2), rep('South Asian/Other Asian', 2), 'African'))
anc_merge <- merge(gwas_simplified, anc_categories, by.x='BROAD.ANCESTRAL.CATEGORY', 'ancestry')
anc_merge$category2 <- factor(anc_merge$category2, levels=(c('European', 'East Asian', 'South Asian/Other Asian', 'African', 'Hispanic/Latino',
                                                                'Greater Middle Eastern', 'Oceanic', 'Other', 'Multiple', 'Not Reported'))) #rev
anc_merge <- subset(anc_merge, category2 != 'Not Reported')
gwas_pop_date_agg <- anc_merge %>%
  select(STUDY.ACCCESSION, PUBMEDID, FIRST.AUTHOR, DATE, STAGE, NUMBER.OF.INDIVDUALS, BROAD.ANCESTRAL.CATEGORY, category2) %>%
  arrange(DATE) %>%
  subset(!is.na(NUMBER.OF.INDIVDUALS)) %>%
  group_by(category2) %>%
  mutate(date_total = cumsum(NUMBER.OF.INDIVDUALS)) %>%
  group_by(DATE, category2) %>%
  slice(which.max(date_total))

# Set colors for population plot
color_vec <- c(brewer.pal(4, 'Set1'), 'grey')
color_vec <- c(color_vec, brewer.pal(3, 'Reds')[2:3], color_vec[2], brewer.pal(5, 'Greens')[2:5], color_vec[4], 'grey')
names(color_vec) <- (c('ASN', 'EUR', 'Non-EURASN', 'Multiple', 'NR', 'East Asian', 'Other Asian', 'European', 'African', 
                      'Hispanic or Latin American', 'MidNatOce', 'Other', 'Multiple', 'NR'))

color_vec <- c(brewer.pal(8, 'Set1'), brewer.pal(3, 'Greys')[2:3])
labels <- levels(anc_merge$category2)
names(color_vec) <- c(labels[1:2], labels[4:3], labels[5:length(labels)])

#anc_merge <- anc_merge %>% arrange(DATE)
my_vals <- gwas_pop_date_agg %>% arrange(DATE) %>% expand(DATE, category2) %>% distinct()
my_vals2 <- merge(my_vals, gwas_pop_date_agg, all.x=T, by=c('DATE', 'category2'))
my_vals2$date_total[2:10] <- 0
my_vals2 <- my_vals2 %>%
  subset(category2!='Not Reported') %>%
  group_by(category2) %>%
  mutate(fill_gap = na.locf(date_total, fromLast=F, na.rm=F)) 
# na.locf is filling from the wrong direction. later dates first
my_vals3 <- my_vals2 %>%
  subset(category2 != 'Not Reported') %>%
  group_by(DATE) %>%
  mutate(pop_frac=fill_gap/sum(fill_gap))

p2 <- ggplot(my_vals2, aes(x=DATE, y=fill_gap/1e6, fill=category2, color=category2)) +
#p2 <- ggplot(my_vals3, aes(x=DATE, y=pop_frac, fill=category2, color=category2)) +
  geom_area(position='stack') +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_fill_manual(values=color_vec, name='Population') +
  scale_color_manual(values=color_vec, name='Population') +
  labs(x='', y='Individuals in GWAS (millions)') +
  #labs(x='Date', y='Fraction of individuals in GWAS') +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90),
        text = element_text(size=16),
        legend.position = c(0.02, 1),
        legend.justification = c(0, 1),
        legend.text = element_text(size=14),
        legend.background = element_rect(fill = "transparent", colour = NA))
ggsave('gwas_time.pdf', p2)

p3 <- ggplot(my_vals3, aes(x=DATE, y=pop_frac, fill=category2, color=category2)) +
  geom_area(position='stack') +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y", position='top') +
  scale_fill_manual(values=color_vec, name='Population') +
  scale_color_manual(values=color_vec, name='Population') +
  labs(x='', y='Fraction') +
  theme_classic() +
  guides(fill=F, color=F) +
  theme(axis.text.x = element_blank(),#element_text(angle=45, hjust=1),
        text = element_text(size=16),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))


# http://www.worldometers.info/world-population/#region
populationAncestries <- data.frame(
  pop=c("East Asian", "South Asian/Other Asian", "European", "Greater Middle Eastern", "African", "Hispanic/Latino", "Oceanic", "DiverseOther", "AfricanEuropeanAdmixed"),
  world=c(1932, 2085, 1145, 410,  1022, 529, 38, NA, NA),
  kgenomes=c(523, 494, 514, 0, 691, 355, 0, NA, NA),
  esp=c(0, 0, 4298, 0, 2217,0,0, NA, NA),
  exac=c(4327, 8256, 36677, 0, 5203, 5789, 0,454, NA),  # Numbers from Monkol on 3/4/15
  pcgcGwas=c(5219, 0, 116766, 0, 0, 0, 0, NA, NA),   # PGC Published 4 (SCZ, BIP, MDD, ADHD)
  ptsd=c(106, 0, 8393, 0, NA, 829, 0, 1295, 9845),
  pop_labs = c('eas', 'sas', 'nfe', 'mde', 'afr', 'amr', 'oce', 'other', 'aa')
)

populationAncestries <- populationAncestries %>%
  mutate(total_world=cumsum(world),
         min_world=total_world - world)
populationAncestries$pop <- factor(populationAncestries$pop, levels=rev(c('Oceanic', 'Greater Middle Eastern', 'Hispanic/Latino', 'African', 'South Asian/Other Asian', 'East Asian', 'European')))
p_global <- ggplot(populationAncestries, aes(x='Present', y=world/1000, fill=pop)) +
  geom_bar(stat='identity') +
  theme_classic() +
  scale_fill_manual(values=color_vec) +
  labs(y='Global population (billions)', x='') +
  guides(fill=F) +
  scale_y_continuous(position='right') +
  theme(axis.text.x = element_text(angle=90),#axis.text.x=element_blank(),
        text = element_text(size=16))

p_global2 <- p_global +
  scale_x_discrete(position='top') +
  labs(y='', x='') +
  theme(axis.text.x = element_blank(),
        axis.text.y= element_text(color='white'),
        axis.ticks.y= element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA))

#legend <- get_legend(p2)
#p3 <- p2 + guides(fill=F, color=F)
p_gwas_global <- plot_grid(p2, p_global, align = "h", rel_widths = c(0.85,0.15))
p_gwas_global2 <- plot_grid(p3, p_global2, align = "h", rel_widths = c(0.85,0.15))
#p_gwas_global <- plot_grid(p3, p_global, legend, align = "v", rel_widths = c(0.65,0.1, 0.25), nrow=1)

p_agg <- ggdraw() +
  draw_plot(p_gwas_global, 0, 0.12, 1, 0.88) +
  draw_plot(p_gwas_global2, 0, 0, 1, 0.24)

save_plot('gwas_time_global3.pdf', p_agg, base_width=7, base_height=5)
#save_plot('gwas_time_global3_wide.pdf', p_agg, base_width=10, base_height=5)
ggsave('gwas_time_global2.pdf', p_gwas_global, width=10)

# Old version with pie charts ---------------------------------------------

# p_global <- ggplot(populationAncestries) +
#   geom_rect(aes(fill=pop, color=pop, ymax=total_world, ymin=min_world, xmax=3, xmin=0)) +
#   xlim(c(0,3)) +
#   coord_polar(theta='y') +
#   geom_text(aes(x=2, y=((total_world+min_world)/2), label=pop)) +
#   labs(x='', y='') +
#   #scale_fill_manual(values=color_vec) +
#   #scale_color_manual(values=color_vec) +
#   theme(aspect.ratio=1) +
#   theme_classic() +
#   guides(fill=F, color=F) +
#   theme(axis.text=element_blank(),
#         axis.line=element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA),
#         axis.ticks.x=element_blank(),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill='transparent', color=NA))
# ggsave('pie_global.pdf', p_global)

# p_gwas_global <- ggdraw() +
#   draw_plot(p2, 0, 0, 1, 1) +
#   draw_plot(p_global, 0.05, 0.5, 0.5, 0.5)

dates <- c('2006-01-01', '2008-01-01', '2010-01-01', '2012-01-01', '2014-01-01', '2016-01-01', '2018-01-01')
pop_date_total <- c()
for(i in 1:length(dates)) {
  date_subset <- subset(gwas_simplified, DATE < dates[i])
  pop_total <- date_subset %>%
    group_by(BROAD.ANCESTRAL.CATEGORY) %>%
    summarize(total=sum(NUMBER.OF.INDIVDUALS, na.rm=T))
  pop_total$date <- dates[i]
  pop_date_total <- rbind(pop_date_total, pop_total)
  print(dim(date_subset))
}


by_2018 <- subset(pop_date_total, date=='2018-01-01')

pie_date <- function(my_date) {
  test <- subset(pop_date_total, date==my_date)
  anc_merge <- merge(test, anc_categories, by.x='BROAD.ANCESTRAL.CATEGORY', 'ancestry') %>%
    arrange(category, category2)
  
  anc_merge$ymin=0
  anc_merge$ymax=anc_merge$total
  for(i in 2:nrow(anc_merge)) {
    anc_merge$ymin[i] = anc_merge$ymin[i-1] + anc_merge$total[i-1]
    anc_merge$ymax[i] = anc_merge$ymin[i] + anc_merge$total[i]
  }
  
  # anc_merge <- anc_merge %>%
  #   gather('super_pop', 'spec_pop', 4:5)
  
  p_year <- ggplot(anc_merge) +
    geom_rect(aes(fill=category, color=category, ymax=ymax, ymin=ymin, xmax=3, xmin=0)) +
    geom_rect(aes(fill=category2, color=category2, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
    xlim(c(0,4)) +
    coord_polar(theta='y') +
    scale_fill_manual(values=color_vec) +
    scale_color_manual(values=color_vec) +
    theme(aspect.ratio=1) +
    theme_classic() +
    guides(fill=F, color=F) +
    theme(axis.text=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill='transparent', color=NA))
  return(p_year)
}
for(i in dates) {
  print(i)
  year <- format(as.Date(i, format="%Y-%m-%d"),"%Y")
  p_year <- pie_date(i)
  ggsave(paste0('pie_', year, '.pdf'), p_year, bg = "transparent", height=4, width=4)
}

date_total <- pop_date_total %>%
  group_by(date) %>%
  summarize(total=sum(total))

# make points/line continuous, but put pie chart at these dates only
p2 <- ggplot(gwas, aes(x=DATE, y=cum_num)) +
  geom_line() +
  theme_classic() +
  ylim(0,2e8) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  #scale_x_continuous(breaks = dates, origin='2005-03-10') +
  labs(x='Date', y='Total individuals in GWAS') +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(size=16))
ggsave('gwas_time.pdf')

date_total$date <- as.Date(date_total$date)
p1 <- ggplot(date_total, aes(x=date, y=total, size=total)) +
  geom_point() +
  labs(x='Date (cumulative)', y='Total # of individuals') +
  theme_classic() +
  scale_x_date(date_labels = '%Y') +
  ylim(0,1.75e8) +
  #xlim('2006-01-01', '2019-01-01') +
  guides(size=F) +
  scale_size_continuous(range=c(5,25)) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(size=16))
ggsave('gwas_size.pdf', p1)
