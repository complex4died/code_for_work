# Step0 Check -------------------------------------------------------------

# install.packages("tidyverse")
packages_needed <- c("tidyverse", "ggplot2", "dplyr", "tidyr", "data.table", "readxl")
new_packages <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]
if(length(new_packages)) 
  install.packages(new_packages)




# Step1 -------------------------------------------------------------------
rm(list = ls())
setwd("/Users/rhonin025/Desktop/Work/test/Reports/罗宁生物扩增子测序报告/02.OTUAnalysis/Core_files/")

library(tidyverse)
library(data.table)
# otu = data.table::fread("otu_table_v1.csv")
otu = data.table::fread("otu_table.txt")
df_sample = data.table::fread("sample_info_v1.csv", header = TRUE) %>% 
  column_to_rownames("V1")
# df3 = fread("tax_table_v1.csv", header = T)

df_otu_id = otu %>% 
  dplyr::select(1, ncol(otu)) %>% 
  separate(taxonomy,
           into = c("Kingdom", "Phylum", "Class", "Order", 
                   "Family", "Genus", "Species"),
           sep = ";\\s*",  
           remove = FALSE) 
table(duplicated(df_otu_id$`#OTU_ID`))

otu = otu %>% 
  select(-ncol(otu)) %>% 
  column_to_rownames("#OTU_ID")


# table(df_otu_id$Kingdom)
# table(df_otu_id$Phylum)
# table(df_otu_id$Class)
# table(df_otu_id$Order)
# table(df_otu_id$Family)
# table(df_otu_id$Genus)
# table(df_otu_id$Species)


df_tax_clean = df_otu_id %>%
  mutate(across(Kingdom:Species, ~ ifelse(str_detect(., "^[a-z]__$|^[a-z]__\\s*$"), NA_character_, .))) %>%
  select(-taxonomy)


# # kingdom
# table(is.na(df_tax_clean$Kingdom))
# otu_kingdom = otu %>%
#   rownames_to_column("ASV_ID") %>%
#   left_join(df_tax_clean %>% select(`#OTU_ID`, Kingdom), by = c("ASV_ID" = "#OTU_ID")) %>% 
#   filter(!is.na(Kingdom))





# # phylum
# table(is.na(df_tax_clean$Phylum))
# otu_phylum = otu %>% 
#   rownames_to_column("ASV_ID") %>%
#   left_join(df_tax_clean %>% select(`#OTU_ID`, Phylum), by = c("ASV_ID" = "#OTU_ID")) %>% 
#   filter(!is.na(Phylum))
# 
# 
# 
# # class
# table(is.na(df_tax_clean$Class))
# otu_class = otu %>% 
#   rownames_to_column("ASV_ID") %>%
#   left_join(df_tax_clean %>% select(`#OTU_ID`, Phylum), by = c("ASV_ID" = "#OTU_ID")) %>% 
#   filter(!is.na(Phylum))
# 
# 
# # order
# table(is.na(df_tax_clean$Order))
# 
# 
# 
# 
# 
# # family
# table(is.na(df_tax_clean$Family))
# 
# 
# 
# 
# 
# # genus
# table(is.na(df_tax_clean$Genus))




otu_merged = otu %>%
  rownames_to_column("ASV_ID") %>%
  pivot_longer(cols = -ASV_ID, names_to = "Sample", values_to = "Counts") %>%
  left_join(df_tax_clean, by = c("ASV_ID" = "#OTU_ID"))


# df-Kingdom
kingdom_table = otu_merged %>%
  mutate(Kingdom = replace_na(Kingdom, "k__Unclassified")) %>%
  group_by(Sample, Kingdom) %>%
  summarise(Total_Counts = sum(Counts), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Total_Counts, values_fill = 0) %>%
  column_to_rownames("Kingdom")
kingdom_rel = kingdom_table %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x)))

tax_col = setdiff(colnames(kingdom_rel), df_sample$SampleID)
kingdom_rel = kingdom_rel %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))
kingdom_table = kingdom_table %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))



# df-Phylum
phylum_table = otu_merged %>%
  mutate(Phylum = replace_na(Phylum, "p__Unclassified")) %>%
  group_by(Sample, Phylum) %>%
  summarise(Total_Counts = sum(Counts), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Total_Counts, values_fill = 0) %>%
  column_to_rownames("Phylum")
phylum_rel = phylum_table %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x)))

tax_col = setdiff(colnames(phylum_rel), df_sample$SampleID)
phylum_rel = phylum_rel %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))
phylum_table = phylum_table %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))



# df-Class
class_table = otu_merged %>%
  mutate(Class = replace_na(Class, "c__Unclassified")) %>%
  group_by(Sample, Class) %>%
  summarise(Total_Counts = sum(Counts), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Total_Counts, values_fill = 0) %>%
  column_to_rownames("Class")
class_rel = class_table %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x)))

tax_col = setdiff(colnames(class_rel), df_sample$SampleID)
class_rel = class_rel %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))
class_table = class_table %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))


# df-Order
order_table = otu_merged %>%
  mutate(Order = replace_na(Order, "o__Unclassified")) %>%
  group_by(Sample, Order) %>%
  summarise(Total_Counts = sum(Counts), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Total_Counts, values_fill = 0) %>%
  column_to_rownames("Order")
order_rel = order_table %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x)))

tax_col = setdiff(colnames(order_rel), df_sample$SampleID)
order_rel = order_rel %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))
order_table = order_table %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))




# df-Family
family_table = otu_merged %>%
  mutate(Family = replace_na(Family, "f__Unclassified")) %>%
  group_by(Sample, Family) %>%
  summarise(Total_Counts = sum(Counts), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Total_Counts, values_fill = 0) %>%
  column_to_rownames("Family")
family_rel = family_table %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x)))

tax_col = setdiff(colnames(family_rel), df_sample$SampleID)
family_rel = family_rel %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))
family_table = family_table %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))


# df-Genus
genus_table = otu_merged %>%
  mutate(Genus = replace_na(Genus, "g__Unclassified")) %>%
  group_by(Sample, Genus) %>%
  summarise(Total_Counts = sum(Counts), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Total_Counts, values_fill = 0) %>%
  column_to_rownames("Genus")
genus_rel = genus_table %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x)))

tax_col = setdiff(colnames(genus_rel), df_sample$SampleID)
genus_rel = genus_rel %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))
genus_table = genus_table %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))



# df-Species
species_table = otu_merged %>%
  mutate(Species = replace_na(Species, "s__Unclassified")) %>%
  group_by(Sample, Species) %>%
  summarise(Total_Counts = sum(Counts), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Total_Counts, values_fill = 0) %>%
  column_to_rownames("Species")
species_rel = species_table %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x)))

tax_col = setdiff(colnames(species_rel), df_sample$SampleID)
species_rel = species_rel %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))
species_table = species_table %>%
  select(all_of(tax_col), all_of(df_sample$SampleID))

save(
  df_sample, df_otu_id, df_tax_clean,
  kingdom_table, kingdom_rel,
  phylum_table, phylum_rel, 
  class_table, class_rel, 
  order_table, order_rel, 
  family_table, family_rel, 
  genus_table, genus_rel, 
  species_table, species_rel,
  file = "step1_output.rdata"
  )





# tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# all_tax_tables <- tax_levels %>%
#   set_names() %>% 
#   map(function(level) {
#     otu_merged %>%
#       group_by(Sample, !!sym(level)) %>% 
#       summarise(Total_Counts = sum(Counts), .groups = "drop") %>%
#       mutate(!!sym(level) := replace_na(!!sym(level), paste0(tolower(substr(level, 1, 1)), "__Unclassified"))) %>%
#       pivot_wider(names_from = Sample, values_from = Total_Counts, values_fill = 0)
#   })
# 
# 
# # df_kingdom = all_tax_tables$Kingdom
# # df_phylum = all_tax_tables$Phylum
# # df_class = all_tax_tables$Class
# # df_order = all_tax_tables$Order
# # df_family =  all_tax_tables$Family
# # df_genus = all_tax_tables$Genus
# # df_species = all_tax_tables$Species
# 
# 
# all_rel_tables <- all_tax_tables %>%
#   map(~ .x %>% mutate(across(where(is.numeric), ~ . / sum(.))))
# 
# 
# # list2env(all_tax_tables, envir = .GlobalEnv)
# names(all_rel_tables) %>%
#   walk(~ assign(paste0("df_", tolower(.x)), all_rel_tables[[.x]], envir = .GlobalEnv))
# 



# Step2 bar ------------------------------------------------------------------
rm(list = ls())
setwd("/Users/rhonin025/Desktop/Work/test/Reports/罗宁生物扩增子测序报告/02.OTUAnalysis/Core_files/")


# kingdom
load("step1_output.rdata")
rm(list = setdiff(ls(), c("kingdom_rel", "df_sample")))

library(tidyverse)
df_long_kingdom = kingdom_rel %>%
  rownames_to_column(var = "Kingdom") %>% 
  pivot_longer(cols = -Kingdom,          
               names_to = "SampleID",  
               values_to = "Relative Abundance")  %>%
  left_join(df_sample, by = "SampleID")


library(ggplot2)
library(ggsci)
library(ggrepel)
library(ggthemes)

# Group/SampleID顺序设置
df_long_kingdom$Group = factor(df_long_kingdom$Group, 
                               levels = c("CK", "GW", "CY", "GW_CY"))

sample_order = colnames(kingdom_rel) 
df_long_kingdom$SampleID = factor(df_long_kingdom$SampleID, 
                                  levels = sample_order)


# 颜色设置
n_kingdom = length(unique(df_long_kingdom$Kingdom))
my_cols = ggsci::pal_lancet(alpha = 0.6)(n_kingdom)

# plot
plot_kingdom = ggplot(data = df_long_kingdom, aes(x = SampleID, 
                                   y = `Relative Abundance`, 
                                   fill = Kingdom)) + 
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_manual(values = my_cols) + 
  facet_wrap(~Group, scales = "free_x", nrow = 1) + 
  labs(x = " ", 
       y = "Relative Abundance",
       title = " ") + 
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.text = element_text(size = 12, face = "italic")
  ) 
plot_kingdom
# ggsave(plot = plot_kingdom, 
#        filename = "plot_kingdom.tiff", 
#        device = "tiff", 
#        dpi = 600, 
#        height = 9, 
#        width = 9)





# Phylum
load("step1_output.rdata")
rm(list = setdiff(ls(), c("phylum_rel", "df_sample")))

library(tidyverse)
df_long_phylum = phylum_rel %>%
  rownames_to_column(var = "Phylum") %>% 
  pivot_longer(cols = -Phylum,          
               names_to = "SampleID",  
               values_to = "Relative Abundance")  %>%
  left_join(df_sample, by = "SampleID")

## Top 计算
top_phylum = df_long_phylum %>%
  group_by(Phylum) %>%
  summarise(mean_rel = mean(`Relative Abundance`)) %>%
  arrange(desc(mean_rel))

top_phylum_id = top_phylum %>%
  slice(1:20) %>% 
  pull(Phylum)

df_plot_phylum = df_long_phylum %>%
  mutate(Phylum = if_else(Phylum %in% top_phylum_id, Phylum, "Others")) %>%
  group_by(SampleID, Phylum, Group) %>%
  summarise(`Relative Abundance` = sum(`Relative Abundance`), .groups = "drop")

df_plot_phylum$Phylum = factor(df_plot_phylum$Phylum, 
                               levels = rev(c(top_phylum_id, "Others")))

# cors
my_cols = ggsci::pal_d3("category20")(length(top_phylum_id))
names(my_cols) = top_phylum_id
others_col = c("Others" = "grey60")
final_cols = c(my_cols, others_col)


# 顺序
df_plot_phylum$Group = factor(df_plot_phylum$Group, 
                               levels = c("CK", "GW", "CY", "GW_CY"))

sample_order = colnames(phylum_rel) 
df_plot_phylum$SampleID = factor(df_plot_phylum$SampleID, 
                                  levels = sample_order)

# plot
plot_phylum = ggplot(data = df_plot_phylum, aes(x = SampleID, 
                                   y = `Relative Abundance`, 
                                   fill = Phylum)) + 
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_manual(values = final_cols) + 
  facet_wrap(~Group, scales = "free_x", nrow = 1) + 
  # scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = " ", 
       y = "Relative Abundance", 
       title = " ") + 
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.text = element_text(size = 12, face = "italic")
  ) + 
  guides(fill = guide_legend(ncol = 1))
plot_phylum

# ggsave(plot = plot_phylum,
#        filename = "plot_phylum.tiff",
#        device = "tiff",
#        dpi = 600,
#        height = 9,
#        width = 9)
save(top_phylum_id, 
     file = "top_phylum_id.rdata")





# Class
load("step1_output.rdata")
rm(list = setdiff(ls(), c("class_rel", "df_sample")))

library(tidyverse)
df_long_class = class_rel %>%
  rownames_to_column(var = "Class") %>% 
  pivot_longer(cols = -Class,          
               names_to = "SampleID",  
               values_to = "Relative Abundance")  %>%
  left_join(df_sample, by = "SampleID")

## Top 计算
top_class = df_long_class %>%
  group_by(Class) %>%
  summarise(mean_rel = mean(`Relative Abundance`)) %>%
  arrange(desc(mean_rel))

top_class_id = top_class %>%
  slice(1:20) %>% 
  pull(Class)

df_plot_class = df_long_class %>%
  mutate(Class = if_else(Class %in% top_class_id, Class, "Others")) %>%
  group_by(SampleID, Class, Group) %>%
  summarise(`Relative Abundance` = sum(`Relative Abundance`), .groups = "drop")

df_plot_class$Class = factor(df_plot_class$Class, 
                             levels = rev(c(top_class_id, "Others")))


# cors
my_cols = ggsci::pal_d3("category20")(length(top_class_id))
names(my_cols) = top_class_id
others_col = c("Others" = "grey60")
final_cols = c(my_cols, others_col)


# 顺序
df_plot_class$Group = factor(df_plot_class$Group, 
                              levels = c("CK", "GW", "CY", "GW_CY"))

sample_order = colnames(class_rel) 
df_plot_class$SampleID = factor(df_plot_class$SampleID, 
                                 levels = sample_order)


# plot
plot_class = ggplot(data = df_plot_class, aes(x = SampleID, 
                                  y = `Relative Abundance`, 
                                  fill = Class)) + 
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_manual(values = final_cols) + 
  facet_wrap(~Group, scales = "free_x", nrow = 1) + 
  # scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = " ", 
       y = "Relative Abundance", 
       title = " ") + 
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 12, face = "italic")
  ) + 
  guides(fill = guide_legend(ncol = 1))
plot_class
# ggsave(plot = plot_class,
#        filename = "plot_class.tiff",
#        device = "tiff",
#        dpi = 600,
#        height = 9,
#        width = 9)
save(top_class_id, 
     file = "top_class_id.rdata")




# Order
load("step1_output.rdata")
rm(list = setdiff(ls(), c("order_rel", "df_sample")))

library(tidyverse)
df_long_order = order_rel %>%
  rownames_to_column(var = "Order") %>% 
  pivot_longer(cols = -Order,          
               names_to = "SampleID",  
               values_to = "Relative Abundance")  %>%
  left_join(df_sample, by = "SampleID")

## Top 计算
top_order = df_long_order %>%
  group_by(Order) %>%
  summarise(mean_rel = mean(`Relative Abundance`)) %>%
  arrange(desc(mean_rel))

top_order_id = top_order %>%
  slice(1:20) %>% 
  pull(Order)

df_plot_order = df_long_order %>%
  mutate(Order = if_else(Order %in% top_order_id, Order, "Others")) %>%
  group_by(SampleID, Order, Group) %>%
  summarise(`Relative Abundance` = sum(`Relative Abundance`), .groups = "drop")

df_plot_order$Order = factor(df_plot_order$Order, 
                             levels = rev(c(top_order_id, "Others")))


# cors
my_cols = ggsci::pal_d3("category20")(length(top_order_id))
names(my_cols) = top_order_id
others_col = c("Others" = "grey60")
final_cols = c(my_cols, others_col)


# 顺序
df_plot_order$Group = factor(df_plot_order$Group, 
                             levels = c("CK", "GW", "CY", "GW_CY"))

sample_order = colnames(order_rel) 
df_plot_order$SampleID = factor(df_plot_order$SampleID, 
                                levels = sample_order)


# plot
plot_order = ggplot(data = df_plot_order, aes(x = SampleID, 
                                              y = `Relative Abundance`, 
                                              fill = Order)) + 
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_manual(values = final_cols) + 
  facet_wrap(~Group, scales = "free_x", nrow = 1) + 
  # scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = " ", 
       y = "Relative Abundance", 
       title = " ") + 
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 12, face = "italic")
  ) + 
  guides(fill = guide_legend(ncol = 1))
plot_order
# ggsave(plot = plot_order,
#        filename = "plot_order.tiff",
#        device = "tiff",
#        dpi = 600,
#        height = 9,
#        width = 9)
save(top_order_id, 
     file = "top_order_id.rdata")




# Family
load("step1_output.rdata")
rm(list = setdiff(ls(), c("family_rel", "df_sample")))

library(tidyverse)
df_long_family = family_rel %>%
  rownames_to_column(var = "Family") %>% 
  pivot_longer(cols = -Family,          
               names_to = "SampleID",  
               values_to = "Relative Abundance")  %>%
  left_join(df_sample, by = "SampleID")

## Top 计算
top_family = df_long_family %>%
  group_by(Family) %>%
  summarise(mean_rel = mean(`Relative Abundance`)) %>%
  arrange(desc(mean_rel))

top_family_id = top_family %>%
  slice(1:20) %>% 
  pull(Family)

df_plot_family = df_long_family %>%
  mutate(Family = if_else(Family %in% top_family_id, Family, "Others")) %>%
  group_by(SampleID, Family, Group) %>%
  summarise(`Relative Abundance` = sum(`Relative Abundance`), .groups = "drop")

df_plot_family$Family = factor(df_plot_family$Family, 
                             levels = rev(c(top_family_id, "Others")))


# cors
my_cols = ggsci::pal_d3("category20")(length(top_family_id))
names(my_cols) = top_family_id
others_col = c("Others" = "grey60")
final_cols = c(my_cols, others_col)


# 顺序
df_plot_family$Group = factor(df_plot_family$Group, 
                             levels = c("CK", "GW", "CY", "GW_CY"))

sample_family = colnames(family_rel) 
df_plot_family$SampleID = factor(df_plot_family$SampleID, 
                                levels = sample_family)


# plot
plot_family = ggplot(data = df_plot_family, aes(x = SampleID, 
                                              y = `Relative Abundance`, 
                                              fill = Family)) + 
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_manual(values = final_cols) + 
  facet_wrap(~Group, scales = "free_x", nrow = 1) + 
  # scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = " ", 
       y = "Relative Abundance", 
       title = " ") + 
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 12, face = "italic")
  ) + 
  guides(fill = guide_legend(ncol = 1))
plot_family
# ggsave(plot = plot_family,
#        filename = "plot_family.tiff",
#        device = "tiff",
#        dpi = 600,
#        height = 9,
#        width = 9)
save(top_family_id, 
     file = "top_family_id.rdata")





# Genus
load("step1_output.rdata")
rm(list = setdiff(ls(), c("genus_rel", "df_sample")))

library(tidyverse)
df_long_genus = genus_rel %>%
  rownames_to_column(var = "Genus") %>% 
  pivot_longer(cols = -Genus,          
               names_to = "SampleID",  
               values_to = "Relative Abundance")  %>%
  left_join(df_sample, by = "SampleID")

## Top 计算
top_genus = df_long_genus %>%
  group_by(Genus) %>%
  summarise(mean_rel = mean(`Relative Abundance`)) %>%
  arrange(desc(mean_rel))

top_genus_id = top_genus %>%
  slice(1:20) %>% 
  pull(Genus)

df_plot_genus = df_long_genus %>%
  mutate(Genus = if_else(Genus %in% top_genus_id, Genus, "Others")) %>%
  group_by(SampleID, Genus, Group) %>%
  summarise(`Relative Abundance` = sum(`Relative Abundance`), .groups = "drop")

df_plot_genus$Genus = factor(df_plot_genus$Genus, 
                               levels = rev(c(top_genus_id, "Others")))


# cors
my_cols = ggsci::pal_d3("category20")(length(top_genus_id))
names(my_cols) = top_genus_id
others_col = c("Others" = "grey60")
final_cols = c(my_cols, others_col)


# 顺序
df_plot_genus$Group = factor(df_plot_genus$Group, 
                              levels = c("CK", "GW", "CY", "GW_CY"))

sample_genus = colnames(genus_rel) 
df_plot_genus$SampleID = factor(df_plot_genus$SampleID, 
                                 levels = sample_genus)


# plot
plot_genus = ggplot(data = df_plot_genus, aes(x = SampleID, 
                                                y = `Relative Abundance`, 
                                                fill = Genus)) + 
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_manual(values = final_cols) + 
  facet_wrap(~Group, scales = "free_x", nrow = 1) + 
  # scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = " ", 
       y = "Relative Abundance", 
       title = " ") + 
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 12, face = "italic")
  ) + 
  guides(fill = guide_legend(ncol = 1))
plot_genus
# ggsave(plot = plot_genus,
#        filename = "plot_genus.tiff",
#        device = "tiff",
#        dpi = 600,
#        height = 9,
#        width = 9)
save(top_genus_id, 
     file = "top_genus_id.rdata")

#### ?: 问， 图片中出现了前面定义的g_Unclassified， 需要删除吗？Un是前面手动赋与的？？

### 删除版
# Genus
load("step1_output.rdata")
rm(list = setdiff(ls(), c("genus_rel", "df_sample")))

library(tidyverse)
df_long_genus = genus_rel %>%
  rownames_to_column(var = "Genus") %>% 
  pivot_longer(cols = -Genus,          
               names_to = "SampleID",  
               values_to = "Relative Abundance")  %>%
  left_join(df_sample, by = "SampleID")

## Top 计算
top_genus = df_long_genus %>%
  group_by(Genus) %>%
  summarise(mean_rel = mean(`Relative Abundance`)) %>%
  arrange(desc(mean_rel))
top_genus = top_genus %>% 
  filter(Genus != "g__Unclassified")

top_genus_id = top_genus %>%
  slice(1:20) %>% 
  pull(Genus)

df_plot_genus = df_long_genus %>%
  mutate(Genus = if_else(Genus %in% top_genus_id, Genus, "Others")) %>%
  group_by(SampleID, Genus, Group) %>%
  summarise(`Relative Abundance` = sum(`Relative Abundance`), .groups = "drop")

df_plot_genus$Genus = factor(df_plot_genus$Genus, 
                             levels = rev(c(top_genus_id, "Others")))


# cors
my_cols = ggsci::pal_d3("category20")(length(top_genus_id))
names(my_cols) = top_genus_id
others_col = c("Others" = "grey60")
final_cols = c(my_cols, others_col)


# 顺序
df_plot_genus$Group = factor(df_plot_genus$Group, 
                             levels = c("CK", "GW", "CY", "GW_CY"))

sample_genus = colnames(genus_rel) 
df_plot_genus$SampleID = factor(df_plot_genus$SampleID, 
                                levels = sample_genus)


# plot
plot_genus = ggplot(data = df_plot_genus, aes(x = SampleID, 
                                              y = `Relative Abundance`, 
                                              fill = Genus)) + 
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_manual(values = final_cols) + 
  facet_wrap(~Group, scales = "free_x", nrow = 1) + 
  # scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = " ", 
       y = "Relative Abundance", 
       title = " ") + 
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.text = element_text(size = 12, face = "italic")
  ) + 
  guides(fill = guide_legend(ncol = 1))
plot_genus
# save(top_genus_id, 
#      file = "top_genus_id.rdata")




# Species
load("step1_output.rdata")
rm(list = setdiff(ls(), c("species_rel", "df_sample")))

library(tidyverse)
df_long_species = species_rel %>%
  rownames_to_column(var = "Species") %>% 
  pivot_longer(cols = -Species,          
               names_to = "SampleID",  
               values_to = "Relative Abundance")  %>%
  left_join(df_sample, by = "SampleID")

## Top 计算
top_species = df_long_species %>%
  group_by(Species) %>%
  summarise(mean_rel = mean(`Relative Abundance`)) %>%
  arrange(desc(mean_rel))

top_species_id = top_species %>%
  slice(1:20) %>% 
  pull(Species)

df_plot_species = df_long_species %>%
  mutate(Species = if_else(Species %in% top_species_id, Species, "Others")) %>%
  group_by(SampleID, Species, Group) %>%
  summarise(`Relative Abundance` = sum(`Relative Abundance`), .groups = "drop")

df_plot_species$Species = factor(df_plot_species$Species, 
                             levels = rev(c(top_species_id, "Others")))


# cors
my_cols = ggsci::pal_d3("category20")(length(top_species_id))
names(my_cols) = top_species_id
others_col = c("Others" = "grey60")
final_cols = c(my_cols, others_col)


# 顺序
df_plot_species$Group = factor(df_plot_species$Group, 
                             levels = c("CK", "GW", "CY", "GW_CY"))

sample_species = colnames(species_rel) 
df_plot_species$SampleID = factor(df_plot_species$SampleID, 
                                levels = sample_species)


# plot
plot_species = ggplot(data = df_plot_species, aes(x = SampleID, 
                                              y = `Relative Abundance`, 
                                              fill = Species)) + 
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_manual(values = final_cols) + 
  facet_wrap(~Group, scales = "free_x", nrow = 1) + 
  # scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = " ", 
       y = "Relative Abundance", 
       title = " ") + 
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.text = element_text(size = 12, face = "italic")
  ) + 
  guides(fill = guide_legend(ncol = 1))
plot_species
# ggsave(plot = plot_species,
#        filename = "plot_species.tiff",
#        device = "tiff",
#        dpi = 600,
#        height = 9,
#        width = 9)
save(top_species_id, 
     file = "top_species_id.rdata")





# Step3 heatmap -------------------------------------------------------------------
rm(list = ls())

# BiocManager::install("ComplexHeatmap")
# install.packages("pheatmap")

library(tidyverse)


# Kingdom
load("step1_output.rdata")
rm(list = setdiff(ls(), c("kingdom_rel", "df_sample")))


library(pheatmap)
annotation_col = as.data.frame(df_sample) %>% 
  select(Group)


mat = as.matrix(kingdom_rel)
# pheatmap::pheatmap(mat, 
#                    cluster_rows = FALSE,
#                    cluster_cols = FALSE,
#                    cellwidth = 12, 
#                    cellheight = 8, 
#                    annotation_col = annotation_col
#                    )

col_fun = circlize::colorRamp2(
  c(0, 0.05, 0.2),
  c("#00468BB2", "white", "#ED0000B2")
)

library(ComplexHeatmap)
ComplexHeatmap::Heatmap(
  matrix = mat, 
  col = col_fun, 
  name = "Relative \nAbundance ",
  cluster_columns = FALSE, 
  cluster_rows = FALSE, 
  row_names_side = "left", 
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  width = unit(12, "cm"), 
  height = unit(6, "cm"),
  )


# 拆分
df_sample$Group = factor(df_sample$Group, levels = c("CK","GW", "CY", "GW_CY"))
group_list = df_sample$Group
ComplexHeatmap::Heatmap(
  matrix = mat, 
  col = col_fun, 
  name = "Relative \nAbundance ",
  # 不聚类
  cluster_columns = FALSE, 
  cluster_rows = FALSE,
  # 调整显示
  row_names_side = "left", 
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  # cell 长宽
  width = unit(12, "cm"), 
  height = unit(6, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.75), 
  # 拆分
  column_split = factor(df_sample$Group, levels = c("CK", "GW", "CY", "GW_CY")), 
  column_title = c("CK", "GW", "CY", "GW_CY"), 
  column_title_gp = gpar(fontsize = 10 , fontface = "bold", fill = "gray90"),
  column_gap = unit(0.5, "cm"), 
  border = TRUE
)


anno_group = HeatmapAnnotation(
  facet = anno_block(
    gp = gpar(fill = "gray90", col = "black"), 
    labels = c("CK", "GW", "CY", "GW_CY"),    # 文字内容
    labels_gp = gpar(col = "black", fontsize = 10, fontface = "bold")
  )
)

ComplexHeatmap::Heatmap(
  matrix = mat, 
  col = col_fun, 
  name = "Relative \nAbundance ",
  # 不聚类
  cluster_columns = FALSE, 
  cluster_rows = FALSE,
  # 调整显示
  row_names_side = "left", 
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  # cell 长宽
  width = unit(15, "cm"), 
  height = unit(5, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.75), 
  # top annotation
  top_annotation = anno_group,
  column_split = group_list,
  column_title = NULL, 
)






## Phylum 
library(tidyverse)
load("step1_output.rdata")
rm(list = setdiff(ls(), c("phylum_rel", "df_sample")))
load("top_phylum_id.rdata")


phylum_rel_top = phylum_rel %>% 
  filter(rownames(phylum_rel) %in% top_phylum_id)


library(ComplexHeatmap)
mat = as.matrix(phylum_rel_top)


col_fun = circlize::colorRamp2(
  c(0, 0.01, 0.05, 0.1, 0.5),
  c("#00468BB2", "#42B540B2", "#EEEEEE", "#FDAF91B2", "#ED0000B2")
)


df_sample$Group = factor(df_sample$Group, levels = c("CK","GW", "CY", "GW_CY"))
group_list = df_sample$Group
anno_group = HeatmapAnnotation(
  facet = anno_block(
    gp = gpar(fill = "gray95", col = "black"), 
    labels = c("CK", "GW", "CY", "GW_CY"),    # 文字内容
    labels_gp = gpar(col = "black", fontsize = 10, fontface = "bold")
  )
)

ht = ComplexHeatmap::Heatmap(
  matrix = mat, 
  col = col_fun, 
  name = "Relative \nAbundance ",
  # 不聚类
  cluster_columns = FALSE, 
  cluster_rows = FALSE,
  # 调整显示
  row_names_side = "left", 
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  # cell 长宽
  width = unit(9, "cm"), 
  height = unit(9, "cm"),
  rect_gp = gpar(col = "white", lwd = 0.75), 
  # top annotation
  top_annotation = anno_group,
  column_split = group_list,
  column_title = NULL, 
)
draw(ht)

pdf("phylum_top_heatmap.pdf", width = 9, height = 9)
draw(ht)
dev.off()






