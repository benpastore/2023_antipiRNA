library(ggpubr)
library(tidyverse)
library(ggplot2)
library(lemon)
library(scales)
library(RColorBrewer)
library(ggrepel)
library(ggpmisc)
library(ggseqlogo)
library(rstatix)
library(ggdist)
library(viridis)
library(ggbeeswarm)

source("/fs/ess/PCON0160/ben/bin/mighty.R")
setwd("/fs/ess/PCON0160/ben/projects/antisense_piRNA_v2")

#################################################
##### Figure 1B ##### 

conditions1 = read.csv("./results_07102023/samples/replicates.csv")
conditions2 = read.csv("./results_07_14_2023/samples/replicates.csv")
conditions = rbind(conditions1, conditions2)
conditions$condition = gsub("par1_prg1SG", "parn1_prg1SG", conditions$condition)

## sense counts
pi1 = read.delim("results_07102023/master_tables/analysis.aligned.v0.m1.count_total_norm.tsv", check.names = F)
pi2 = read.delim("./results_07_14_2023/master_tables/analysis.aligned.v0.m1.count_total_norm.tsv")
pi = pi1 %>% full_join(pi2, by = c("gene", "biotype", "feature", "class"))

pi_grouped = pi %>% 
  filter(feature == "piRNA") %>% 
  pivot_longer(!c(gene, biotype, feature, class), names_to = 'sample', values_to = 'count') %>% 
  group_by(sample, feature) %>% 
  summarise(count = sum(count, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(conditions, by = c("sample" = 'simple_name')) %>% 
  group_by(feature, condition) %>% 
  summarise(M = mean(count), S = sd(count))

pi_grouped %>% filter(condition %in% c("N2_YD", "parn1_YD", "TW29") ) %>% arrange(desc(condition))

p1 = ggplot(pi_grouped %>% filter(condition %in% c("N2_YD", "parn1_YD", "TW29") ) %>% arrange(desc(condition)), 
            aes(x = condition, y = M)) + 
  geom_bar(stat = 'identity', fill = "black", color = NA, width = 0.5) + 
  geom_errorbar(aes(x = condition, ymin = M - S, ymax = M + S), width = 0.3, color = 'red') +
  my_theme() + 
  scale_y_continuous(limits = c(0,60000), breaks = seq(0,60000,by = 10000))
p1

pi %>% 
  filter(feature == "piRNA") %>% 
  pivot_longer(!c(gene, biotype, feature, class), names_to = 'sample', values_to = 'count') %>% 
  group_by(sample, feature) %>% 
  summarise(count = sum(count, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(conditions, by = c("sample" = 'simple_name')) %>% 
  filter(condition %in%  c("parn1_YD", "N2_YD")) %>% 
  mutate(condition = as.character(condition)) %>% 
  t_test(count ~ condition, var.equal = T) %>% 
  adjust_pvalue(method = 'none') %>% 
  add_significance() %>% 
  filter(group1 == "parn1_YD" | group2 == "parn1_YD")

## antisense counts
anti1 = read.delim("./antisense_piRNA_counts.tsv", check.names = F)

wide = anti1 %>% 
  select(-feature) %>% 
  filter(grepl("YD", sample) & !grepl("Oxi", sample)) %>% 
  pivot_wider(names_from = sample, values_from = count) 

wide[is.na(wide)] <- 0

wide %>% filter(parn_YD1 > 0 | parn_YD2 > 0)

write.table(wide, 'anti_rpm_publication.tsv', sep = "\t", col.names = T, row.names = F, quote = F)

anti1 = anti1 %>% 
  group_by(gene, sample) %>% 
  summarise(count = sum(count)) %>% 
  ungroup()

anti2 = read.delim("results_07_14_2023//bed/antipiRNA_bed_master_filtered.tsv", check.names = F)
anti2 %>% as_tibble()
anti2 = anti2 %>% 
  group_by(gene, sample) %>% 
  summarise(count = sum(count)) %>% 
  ungroup()

anti = rbind(anti1, anti2)

anti_grouped = anti %>% 
  group_by(sample) %>% 
  summarise(count = sum(count, na.rm = T)) %>% 
  ungroup() %>% 
  full_join(conditions, by = c("sample" = 'simple_name')) %>% 
  group_by(condition) %>% 
  summarise(M = mean(count), S = sd(count)) %>% 
  replace_na(list(M = 0))

anti_grouped %>% filter(condition %in% c("N2_YD", "parn1_YD", "TW29") ) %>% arrange(desc(condition))


p2 = ggplot(anti_grouped %>% filter(condition %in% c("N2_YD", "parn1_YD", "TW29") ) %>% arrange(desc(condition)), 
            aes(x = condition, y = M)) + 
  geom_bar(stat = 'identity', fill = "black", color = NA, width = 0.5) + 
  geom_errorbar(aes(x = condition, ymin = M - S, ymax = M + S), width = 0.3, color = 'red') +
  my_theme() + 
  scale_y_continuous(limits = c(0,600), breaks = seq(0,600,by = 100))

p2

p = ggarrange(p1, p2)
ggsave(p, filename = "./piRNA_anti_level_parn1prg1.pdf", dpi = 300, height = 6, width = 12, device = cairo_pdf)
p


#################################################
##### Figure 1D and E ##### 

# antipiRNA level from parn-1 oxi
anti_dat = read.delim("./antisense_piRNA_counts.tsv")
sense_dat =  read.delim("./results_07102023/master_tables/analysis.aligned.v0.m1.count_total_norm.tsv")
sample1 = read.csv("results_07102023/samples/replicates.csv", colClasses = c("character", "character"))
sample2 = read.csv("results_07_14_2023/samples/replicates.csv", colClasses = c("character", "character"))
sample3 = read.csv("results_SNAP_ligation_07_24_2023/samples/replicates.csv", colClasses = c("character", "character"))
samples = rbind(sample1, sample2, sample3)

sense_piRNA = sense_dat %>% 
  pivot_longer(!c(gene, biotype, feature, class), names_to = 'sample', values_to = 'count') %>% 
  group_by(sample, feature) %>% 
  summarise(count = sum(count)) %>% 
  filter(feature == "piRNA") %>% 
  ungroup() %>% 
  left_join(samples, c("sample" = "simple_name")) %>% 
  group_by(condition) %>% 
  summarise(M = mean(count), S = sd(count)) %>% 
  mutate(feature = "piRNA")

sense_dat %>% 
  pivot_longer(!c(gene, biotype, feature, class), names_to = 'sample', values_to = 'count') %>% 
  group_by(sample, feature) %>% 
  summarise(count = sum(count)) %>% 
  filter(feature == "piRNA") %>% 
  ungroup() %>% 
  left_join(samples, c("sample" = "simple_name")) %>% 
  filter(condition %in% c("parn1_YD", "parn1_oxi") ) %>% 
  group_by(feature) %>% 
  t_test(count ~ condition) %>% 
  add_significance()

anti_piRNA = anti_dat %>% 
  group_by(sample, feature) %>% 
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  left_join(samples, c("sample" = "simple_name")) %>% 
  group_by(condition) %>% 
  summarise(M = mean(count), S = sd(count)) %>% 
  mutate(feature = "antipiRNA")

anti_dat %>% 
  group_by(sample, feature) %>% 
  summarise(count = sum(count)) %>% 
  ungroup() %>% 
  left_join(samples, c("sample" = "simple_name")) %>% 
  filter(condition %in% c("parn1_YD", "parn1_oxi") ) %>% 
  group_by(feature) %>% 
  t_test(count ~ condition, var.equal = T) %>% 
  add_significance()

dat_all = rbind(sense_piRNA, anti_piRNA)

p = ggplot(dat_all %>% filter(condition %in% c("parn1_YD", "parn1_oxi") ) %>% arrange(desc(condition)), 
           aes(x = condition, y = M)) + 
  geom_bar(stat = 'identity', fill = "black", color = NA, width = 0.5) + 
  geom_errorbar(aes(x = condition, ymin = M - S, ymax = M + S), width = 0.3, color = 'red') +
  my_theme() + 
  facet_wrap(~feature, scales = "free")
p
ggsave(p, filename = 'antilevel.pdf', dpi = 300, height = 6, width = 6, device = cairo_pdf)

#################################################
##### Figure 1F, 1G, S1A, S1B ##### 

#1F/G anti-piRNA length distribution
data = read.delim("./length_distribution_unfiltered_piRNA.tsv")

sample1 = read.csv("results_07102023/samples/replicates.csv", colClasses = c("character", "character"))
sample2 = read.csv("results_07_14_2023/samples/replicates.csv", colClasses = c("character", "character"))
sample3 = read.csv("results_SNAP_ligation_07_24_2023/samples/replicates.csv", colClasses = c("character", "character"))
samples = rbind(sample1, sample2, sample3)

data_grouped = data %>% 
  left_join(samples, by = c("sample" = "simple_name")) %>% 
  ungroup()

conditions = c("N2_input", "N2_prg1ip", "TW28_input", "TW28_prg1ip")

cols = c("A" = "#d55e00", "G" = "#009273", "C" = "#e69f00", "T" = "#56b4e9")

p1 = ggplot(data = data_grouped %>% filter(condition %in% c("TW28_input", "TW28_prg1ip")), 
            aes(x = length, y = count, fill = first)) + 
  geom_bar(stat = 'identity', color = 'grey70', size = 0.2) + 
  my_theme() + 
  theme(aspect.ratio = 1) + 
  facet_wrap(~condition) + 
  scale_fill_manual(values = cols) + 
  scale_x_continuous(limits = c(14, 26), breaks = seq(15,26, by = 2)) + 
  scale_y_continuous(limits = c(0,4000), breaks = seq(0,4000, by = 1000))
p1

p2 = ggplot(data = data_grouped %>% filter(condition %in% c("N2_input", "N2_prg1ip")), 
            aes(x = length, y = count, fill = first)) + 
  geom_bar(stat = 'identity', color = 'grey70', size = 0.2) + 
  my_theme() + 
  theme(aspect.ratio = 1) + 
  facet_wrap(~condition) + 
  scale_fill_manual(values = cols) + 
  scale_x_continuous(limits = c(14, 26), breaks = seq(15,26, by = 2)) + 
  scale_y_continuous(limits = c(0,2000), breaks = seq(0,2000, by = 500))
p2

p = ggarrange(p2, p1, nrow = 2, common.legend = T)
ggsave(p, filename = "anti_lendist.pdf", dpi = 300, height = 7, width = 7, device = cairo_pdf)

#S1A/B sense piRNA length distribution
data = read.delim("./results_08_30_2023/piRNA_length_distribution.tsv")
levels(factor(data$sample))

cols = c("A" = "#d55e00", "G" = "#009273", "C" = "#e69f00", "T" = "#56b4e9")

data_grouped = data %>% 
  group_by(sample, firstnt, seqlen) %>% 
  summarise(rpm = sum(rpm))

p1 = ggplot(data = data_grouped %>% filter(sample %in% c("TW28prgip_raw", "TW28input_raw")), 
            aes(x = seqlen, y = rpm, fill = firstnt)) + 
  geom_bar(stat = 'identity', color = 'grey70', size = 0.2) + 
  my_theme() + 
  theme(aspect.ratio = 1) + 
  scale_fill_manual(values = cols) + 
  scale_x_continuous(limits = c(14,33), breaks = seq(15,33,by = 3)) +
  scale_y_continuous(limits = c(0,200000), breaks = seq(0,200000, by = 50000)) +
  facet_wrap(~sample) 
p1

p2 = ggplot(data = data_grouped %>% filter(sample %in% c("N2input_raw", "N2prgip_raw")), 
            aes(x = seqlen, y = rpm, fill = firstnt)) + 
  geom_bar(stat = 'identity', color = 'grey70', size = 0.2) + 
  my_theme() + 
  theme(aspect.ratio = 1) + 
  scale_fill_manual(values = cols) + 
  scale_x_continuous(limits = c(14,33), breaks = seq(15,33,by = 3)) +
  scale_y_continuous(limits = c(0,9e5), breaks = seq(0,9e5, by = 200000)) +
  facet_wrap(~sample) 
p2

p = ggarrange(p2, p1, nrow = 2, common.legend = T)
ggsave(p, filename = "prg1IP_lendist_piRNA.pdf", dpi = 300, height = 7, width = 7, device = cairo_pdf)


#################################################
##### Figure 3B, S3A, S3B

# parn1 rde3 tailing
tailor_piRNA = read.delim("/fs/ess/PCON0160/ben/projects/antisense_piRNA_v2/results_normalize_miRNA_07_31_2023/tailor_piRNA.tsv", check.names = F)
perfect = read.delim("/fs/ess/PCON0160/ben/projects/antisense_piRNA_v2/results_normalize_miRNA_07_31_2023/master_tables/analysis.aligned.v0.m1.count_total_norm.tsv", 
                     check.names = F) %>% filter(feature == "piRNA")
condition = read.delim("/fs/ess/PCON0160/ben/projects/antisense_piRNA_v2/results_normalize_miRNA_07_31_2023/samples/replicates.csv", sep = ",", col.names = c("sample", "condition"))

tailor_piRNA %>% 
  filter(condition == "N2_YD" | condition == "parn1_rde3") %>% 
  separate(seq, c("seq", "tail"), sep = ":") 

conditions = c("N2_YD", "parn1_YD", "parn1_rde3")

p = plot_tails(tailor_piRNA, perfect %>% filter(feature == "piRNA"), condition, conditions) + 
  theme(aspect.ratio = 1.3) + 
  scale_y_continuous(limits = c(0, 15), breaks = seq(0,15,by=1.5)) + 
  labs(x = "", y = "Tailing Frequency %")
p
ggsave(p, filename = 'tailing_parn1rde3disl2.pdf', dpi = 300, height = 6, width = 6, device = cairo_pdf)

res = calculate_tail_frequency(tailor_piRNA %>% filter(condition %in% conditions), 
                               perfect %>% filter(feature == "piRNA"),
                               condition, group = 'total')

colSums(res %>% filter(condition == "parn1_YD") %>% select(freq))

res %>% filter(condition == "parn1_YD")

# find -1 nucleotide of GU tail before initial G
res = calculate_tail_frequency(tailor_piRNA %>% filter(condition %in% conditions), 
                               perfect %>% filter(feature == "piRNA"),
                               condition, group = "gene_seq")

minus_bias = res %>% 
  select(seq, tail, tail_group, condition) %>% 
  distinct() %>% 
  mutate(tail_length = nchar(tail)) %>% 
  mutate(seq_minus_tail = substr(seq, 1, nchar(seq)-tail_length)) %>% 
  mutate(last_nt = substr(seq_minus_tail, nchar(seq_minus_tail), nchar(seq_minus_tail))) %>% 
  mutate(second_last_nt = substr(seq_minus_tail, nchar(seq_minus_tail)-1, nchar(seq_minus_tail)-1)) 

minus_bias$tail_group_f = factor(minus_bias$tail_group, levels = c("A", "C", "T", "G", "GU_repeat", "UG_repeat", "Other"))

minus_1_bias = minus_bias %>% 
  group_by(last_nt, tail_group_f, condition) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  group_by(tail_group_f, condition) %>% 
  mutate(perc = n / sum(n))

minus_2_bias = minus_bias %>% 
  group_by(second_last_nt, tail_group_f, condition) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  group_by(tail_group_f, condition) %>% 
  mutate(perc = n / sum(n))

cols = c("A" = "#d55e00",
         "G" = "#009273",
         "C" = "#f0e442",
         "T" = "#56b4e9")

p1 = ggplot(data = minus_1_bias %>% filter(condition == "parn1_YD"), aes(x = tail_group_f, y = perc)) + 
  geom_bar(stat = 'identity', aes(fill = last_nt)) + 
  scale_fill_manual(values = cols) + 
  facet_wrap(~condition, ncol = 1) + 
  my_theme() + 
  theme(aspect.ratio = 0.8)
p1

p2 = ggplot(data = minus_2_bias %>% filter(condition == "parn1_YD"), aes(x = tail_group_f, y = perc)) + 
  geom_bar(stat = 'identity', aes(fill = second_last_nt)) + 
  scale_fill_manual(values = cols) + 
  facet_wrap(~condition, ncol = 1) + 
  my_theme() + 
  theme(aspect.ratio = 0.8)
p2

p = ggarrange(p1, p2, nrow = 1, common.legend = T)
p

ggsave(p, filename = 'minus_1_nt_bias_tailed_parn1.pdf', dpi = 300, height = 10, width = 10, device = cairo_pdf)

#################################################
##### Figure 3D and S3C

# 22G RNA level S3C
conditions1 = read.csv("./results_07102023/samples/replicates.csv")
conditions2 = read.csv("./results_07_14_2023/samples/replicates.csv")
samples = rbind(conditions1, conditions2)

d = read.delim("./results_normalize_miRNA_07_31_2023/master_tables/analysis.aligned.v0.m1.count_total_norm.tsv", check.names = F)
conditions = c("parn1_YD", "parn1_rde3", "parn1_ego1RNAi", "parn1_rrf1", "parn1_ekl1RNAi", "parn1_drh3", "parn1_drh3RNAi", "parn1_oma1RNAi")

d_grouped = d %>% 
  filter(!grepl("ego-1|drh-3|ekl-1|oma-1", gene)) %>% 
  pivot_longer(!c(gene, biotype, feature, class), names_to = 'sample', values_to = 'count') %>% 
  group_by(sample, feature) %>% 
  summarise(count = sum(count, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(samples, by = c("sample" = 'simple_name')) %>% 
  group_by(feature, condition) %>% 
  summarise(M = mean(count), S = sd(count)) %>% 
  filter(condition %in% conditions) %>% 
  filter(feature == "siRNA")

d_grouped

d_grouped = d_grouped %>% filter(condition %in% conditions ) %>% filter(feature == "siRNA") %>% arrange(desc(condition))
d_grouped = d_grouped %>% filter(condition %in% conditions)
d_grouped$condition_f = factor(d_grouped$condition, levels = conditions)
p = ggplot(d_grouped, 
           aes(x = condition_f, y = M)) + 
  geom_bar(stat = 'identity', fill = "black", color = NA, width = 0.5) + 
  geom_errorbar(aes(x = condition, ymin = M - S, ymax = M + S), width = 0.3, color = 'red') +
  my_theme() + 
  theme(aspect.ratio = 0.6, axis.text.x = element_text(size = 5))
p
ggsave(p, filename = "22G_level.pdf", dpi = 300, height = 7, width = 12, device = cairo_pdf)

d %>% 
  filter(!grepl("ego-1|drh-3|ekl-1|oma-1", gene)) %>% 
  pivot_longer(!c(gene, biotype, feature, class), names_to = 'sample', values_to = 'count') %>% 
  group_by(sample, feature) %>% 
  summarise(count = sum(count, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(samples, by = c("sample" = 'simple_name')) %>% 
  filter(feature == "siRNA") %>% 
  filter(condition %in% conditions) %>% 
  mutate(condition = as.character(condition)) %>% 
  t_test(count ~ condition, var.equal = T) %>% 
  adjust_pvalue(method = 'none') %>% 
  add_significance() %>% 
  filter(group1 == "parn1_YD" | group2 == "parn1_YD")

## 3D antipiRNA level from parn-1 mutants normalized to miRNAs
dat = read.delim("./antipiRNA_bed_master_filtered_miRNA_norm.counts.tsv")

sample1 = read.csv("results_07102023/samples/replicates.csv", colClasses = c("character", "character"))
sample2 = read.csv("results_07_14_2023/samples/replicates.csv", colClasses = c("character", "character"))
sample3 = read.csv("results_SNAP_ligation_07_24_2023/samples/replicates.csv", colClasses = c("character", "character"))
samples = rbind(sample1, sample2, sample3)

dat_grouped = dat %>% 
  group_by(sample) %>% 
  summarise(count = sum(count)) %>% 
  left_join(samples, by = c("sample" = "simple_name")) %>% 
  group_by(condition) %>% 
  summarise(M = mean(count), S = sd(count))

conditions = c("parn1_YD", "parn1_rde3", "parn1_ego1RNAi", "parn1_rrf1", "parn1_ekl1RNAi", "parn1_drh3", "parn1_drh3RNAi", "parn1_oma1RNAi")
sub = dat_grouped %>% filter(condition %in% conditions)
sub$condition_f = factor(sub$condition, levels = conditions)
p = ggplot(sub, aes(x = condition_f, y = M)) + 
  geom_bar(stat = 'identity', fill = "black", color = NA, width = 0.5) + 
  geom_errorbar(aes(x = condition, ymin = M - S, ymax = M + S), width = 0.3, color = 'red') +
  my_theme() + 
  theme(aspect.ratio = 0.6, axis.text.x = element_text(size = 5))
p
ggsave(p, filename = "antisense_piRNA_levels.pdf", dpi = 300, height = 7, width = 12, device = cairo_pdf)

dat_grouped

dat %>% 
  group_by(sample) %>% 
  summarise(count = sum(count)) %>% 
  left_join(samples, by = c("sample" = "simple_name")) %>% 
  filter(condition %in% conditions ) %>% 
  t_test(count ~ condition, var.equal = T) %>% 
  adjust_pvalue(method = 'none') %>% 
  add_significance() %>% 
  filter(group1 == "parn1_YD" | group2 == "parn1_YD")

anti_piRNA = dat %>% 
  group_by(sample) %>% 
  summarise(count = sum(count))

#################################################
##### Figure 4B, and S4A
dist_dat = read.delim("antisense_to_sense_distances.tsv")

dist_dat$type = as.character(dist_dat$type)

sample1 = read.csv("results_07102023/samples/replicates.csv", colClasses = c("character", "character"))
sample2 = read.csv("results_07_14_2023/samples/replicates.csv", colClasses = c("character", "character"))
sample3 = read.csv("results_SNAP_ligation_07_24_2023/samples/replicates.csv", colClasses = c("character", "character"))
samples = rbind(sample1, sample2, sample3)
samples$condition = gsub("par1_prg1SG", "parn1_prg1SG", samples$condition)

d_grouped = dist_dat %>% 
  left_join(samples, by = c("sample" = "simple_name")) %>% 
  group_by(condition, type, distance) %>% 
  summarise(M = mean(count), S = sd(count)) %>% 
  ungroup() %>% 
  group_by(condition, type) %>% 
  mutate(percent = (M / sum(M)) * 100) %>% 
  ungroup()

conditions = c("parn1_YD", "parn1_prg1DA", "parn1_prg1SG", "parn1_ego1RNAi")
p1 = ggplot(data = d_grouped %>% filter(type == "5p3p") %>% filter(condition %in% conditions), aes(x = distance, y = percent)) + 
  geom_line(aes(color = condition), lwd = 0.8) + 
  geom_point(aes(color = condition), fill = "white", size=2, shape=21, stroke = 1) + 
  scale_x_continuous(limits = c(0, 15), breaks = seq(1,15,by=2)) +
  scale_y_continuous(limits = c(0,80), breaks = seq(0,80,by=10)) + 
  my_theme() + 
  theme(aspect.ratio = 1)

p2 = ggplot(data = d_grouped %>% filter(type == "5p5p") %>% filter(condition %in% conditions), aes(x = distance, y = percent)) + 
  geom_line(aes(color = condition), lwd = 0.8) + 
  geom_point(aes(color = condition), fill = "white", size=2, shape=21, stroke = 1) + 
  scale_x_continuous(limits = c(18, 36), breaks = seq(18,36,by=3)) +
  scale_y_continuous(limits = c(0,80), breaks = seq(0,80,by=10)) + 
  my_theme() + 
  theme(aspect.ratio = 1)

p = ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = T)
p
ggsave(p, filename = "parn1SG_distance.pdf", dpi = 300, height = 8, width = 8, device = cairo_pdf)

# overlay 5'5' with sense piRNA length distribution
len_data = read.delim("./piRNA_length_distribution.tsv")

len_data_grouped = len_data %>% 
  left_join(samples, by = c("sample" = "simple_name")) %>% 
  group_by(sample, condition, seqlen) %>% 
  summarise(rpm = sum(rpm)) %>% 
  ungroup() %>% 
  group_by(condition) %>% 
  mutate(percent  = 100*(rpm/sum(rpm))) %>% 
  ungroup() %>% 
  group_by(condition, seqlen) %>% 
  summarise(M = mean(percent), S = sd(percent)) %>% 
  ungroup() %>% 
  filter(condition == "parn1_YD") %>% 
  dplyr::rename(distance = seqlen) %>% 
  mutate(percent = M)

len_data_grouped['type'] = 'piRNA_lengthdist'

lendat_dist = rbind(d_grouped, len_data_grouped)

p = ggplot(data = lendat_dist %>% filter(type == "5p5p" | type == 'piRNA_lengthdist') %>% filter(condition %in% c("parn1_YD")), aes(x = distance, y = percent)) + 
  geom_line(aes(color = type), lwd = 0.8) + 
  geom_point(aes(color = type), fill = "white", size=2, shape=21, stroke = 0.5) + 
  scale_x_continuous(limits = c(18, 35), breaks = seq(20,35,by=5)) +
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,by=10)) + 
  my_theme() + 
  theme(aspect.ratio = 1)
p
ggsave(p, filename = '5p5pvslengthdist.pdf', dpi = 300, height = 6, width = 6, device = cairo_pdf)

#################################################
##### Figure S4B

# piRNA vs. anti-piRNA level from parn-1 vs. parn-1 prg-1
conditions1 = read.csv("./results_07102023/samples/replicates.csv")
conditions2 = read.csv("./results_07_14_2023/samples/replicates.csv")
conditions = rbind(conditions1, conditions2)
conditions$condition = gsub("par1_prg1SG", "parn1_prg1SG", conditions$condition)

pi1 = read.delim("./results_normalize_miRNA_07_31_2023/master_tables/analysis.aligned.v0.m1.count_all_miRNA_norm.tsv", check.names = F)

pi2 = read.delim("./results_07_14_2023/master_tables/analysis.aligned.v0.m1.count_total_norm.tsv")

pi = pi1 %>% full_join(pi2, by = c("gene", "biotype", "feature", "class"))

pi_grouped = pi %>% 
  filter(feature == "piRNA") %>% 
  pivot_longer(!c(gene, biotype, feature, class), names_to = 'sample', values_to = 'count') %>% 
  group_by(sample, feature) %>% 
  summarise(count = sum(count, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(conditions, by = c("sample" = 'simple_name')) %>% 
  group_by(feature, condition) %>% 
  summarise(M = mean(count), S = sd(count))

p1 = ggplot(pi_grouped %>% filter(condition %in% c("parn1_YD", "parn1_prg1DA", "parn1_prg1SG") ) %>% arrange(desc(condition)), 
            aes(x = condition, y = M)) + 
  geom_bar(stat = 'identity', fill = "black", color = NA, width = 0.5) + 
  geom_errorbar(aes(x = condition, ymin = M - S, ymax = M + S), width = 0.3, color = 'red') +
  my_theme()
p1
ggsave(p1, filename = 'parn1_prg1DASG.pdf', dpi = 300, height = 6, width = 6, device = cairo_pdf)


#################################################
##### Figure 5A

anti_dat = read.delim("./antisense_piRNA_counts.tsv")

anti_dat = anti_dat %>% 
  filter(grepl("parn_YD", sample)) %>% 
  group_by(gene) %>% 
  summarise(anti_rpm = mean(count))


pi_dat = read.delim("./results_07102023/master_tables/analysis.aligned.v0.m1.count_total_norm.average.tsv") 

pi_dat = pi_dat %>% 
  filter(feature == "piRNA") %>% 
  select(gene, parn1_YD) %>% 
  dplyr::rename(sense_rpm = parn1_YD)

dat_merge = anti_dat %>% 
  left_join(pi_dat, by = 'gene') %>%
  filter(sense_rpm > 0)

dat_merge
cor.test(log2(dat_merge$anti_rpm), log2(dat_merge$sense_rpm), method = 'pearson')

p = xy_scat(dat_merge, 'sense_rpm', 'anti_rpm', log2_transform = T, cor_method = 'pearson', axmin = 2^-5, axmax = 2^7, fold_change = 2)
p
ggsave(p, filename = 'anti_sense_cor.png', dpi = 300, height = 5, width = 5)

#################################################
##### Figure 5B

pick_control = function(values, df, column){
  
  compare_df = data.frame(count = values)
  
  med = median(values)
  iqr = IQR(values)
  quants = quantile(values)
  
  print(iqr)
  print(quants)
  
  x = data.frame(count = values)
  
  N = nrow(x %>% filter(count >= med & count <= quants[4]))
  
  df['count'] = log2(df[paste0(column)])
  
  for (i in seq(1,1000,by = 1)) {
    
    N = nrow(compare_df %>% filter(count >= med & count <= quants[4]))
    sample1 = df %>% filter(count >= med & count <= quants[4]) %>% sample_n(N, replace = T)
    
    N = nrow(compare_df %>% filter(count <= med & count >= quants[2]))
    sample2 = df %>% filter(count <= med & count >= quants[2]) %>% sample_n(N, replace = T)
    
    N = nrow(compare_df %>% filter(count >= quants[4] & count <= quants[5]))
    sample3 = df %>% filter(count >= quants[4] & count <= quants[5]) %>% sample_n(N, replace = T)
    
    N = nrow(compare_df %>% filter(count <= quants[2] & count >= quants[1]))
    sample4 = df %>% filter(count <= quants[2] & count >= quants[1]) %>% sample_n(N, replace = T)
    
    N = nrow(compare_df %>% filter(count >= quants[5]))
    sample5 = df %>% filter(count >= quants[5]) %>% sample_n(N, replace = T)
    
    top_outlier = 1.5*iqr + quants[4]
    N_top_outleir = nrow(compare_df %>% filter(count >= top_outlier))
    upper_outlier = df %>% filter(count >= top_outlier) %>% sample_n(N_top_outleir, replace = T)
    
    bottom_outlier = quants[1] - 1.5*iqr
    N_bottom_outleir = nrow(compare_df %>% filter(count <= bottom_outlier))
    bottom_outlier = df %>% filter(count >= top_outlier) %>% sample_n(N_bottom_outleir, replace = T)
    
    sampling = rbind(sample1, sample2, sample3, sample4, upper_outlier, bottom_outlier)
    
    pv = t.test(sampling$count, compare_df$count, var.equal = T, alternative = 'two.sided')$p.value
    #print(pv)
    if (pv >= 0.05) {
      
      sampling = sampling %>% select(-count)
      sampling$sampleN = paste0("sample_",i)
      
      if (i == 1){
        res = sampling
      } else {
        res = rbind(sampling, res)
      }
    } else {
      i = i - 1
    }
  } 
  
  return(res)
  
}
clash = read.delim("/fs/ess/PAS1473/deep_sequencing/prg1CLASH/masterTable.piRNACLASH.Apr.r2.bed.dG15.counts.tsv", header = F, col.names = c("locus_id", "count"))

clash_grouped = clash %>% 
  group_by(locus_id) %>% 
  summarise(clash_count = sum(count, na.rm = T))

antipiRNAs = read.delim("./antipiRNA_bed_master_filtered.tsv")
parn1 = antipiRNAs %>% filter(grepl("parn_YD", sample)) %>% select(gene) %>% distinct()

head(parn1)
nrow(parn1)

piRNA_level = read.delim("./results_07102023/dge/N2_YD_vs_parn1_YD.total_norm.tsv")
piRNA_level = piRNA_level %>% mutate(antisense = ifelse(gene %in% parn1$gene, T, F)) %>% filter(feature == "piRNA")

anti = piRNA_level %>% filter(antisense == T) %>% filter(N2_YD > 0) %>% separate(gene, c("gene_name", "seq_id", "locus_id"), sep = ",") %>% filter(locus_id %in% clash_grouped$locus_id)
non = piRNA_level %>% filter(antisense == F) %>% filter(N2_YD > 0) %>% separate(gene, c("gene_name", "seq_id", "locus_id"), sep = ",") %>% filter(locus_id %in% clash_grouped$locus_id)

nrow(anti)
length(log2(anti$N2_YD))

boxplot.stats(log2(anti$N2_YD))

control = pick_control(log2(anti$N2_YD), non, 'N2_YD')

anti$sampleN = "antisense"

control_clash_counts = control %>% 
  filter(feature == "piRNA") %>% 
  # separate(gene, c("gene_name", "seq_id", "locus_id"), sep = ",") %>% 
  left_join(clash_grouped, by = 'locus_id') %>% 
  drop_na() %>% 
  group_by(sampleN) %>% 
  summarise(clash_count = median(clash_count)) %>% 
  dplyr::rename(sample = sampleN) %>% 
  mutate(antisense = F)

anti_clash_counts = anti %>% 
  filter(feature == "piRNA") %>% 
  #separate(gene, c("gene_name", "seq_id", "locus_id"), sep = ",") %>% 
  left_join(clash_grouped, by = 'locus_id') %>% 
  select(locus_id, clash_count) %>% 
  dplyr::rename(sample = locus_id) %>% 
  mutate(antisense = T)

master = rbind(control_clash_counts, anti_clash_counts)

plot_boxplot(master, 
             samples_col = "antisense", 
             counts_col = "clash_count", distribution = T, pvals = T, ylog2 = T, dots = T)

p2 = ggplot(master, aes(x = antisense, y = clash_count)) + 
  stat_boxplot(geom = "errorbar", width = 0.3) + 
  geom_boxplot() + 
  my_theme() + 
  scale_y_continuous(trans = 'log2') + 
  stat_compare_means(method = 't.test')

p2
ggsave(p2, filename = 'clash.pdf', dpi = 300, height = 6, width = 6, device = cairo_pdf)











