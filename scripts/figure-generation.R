library(tidyverse)
library(ggh4x)
library(RColorBrewer)
library(patchwork)
library(cowplot)

######################### 
# preprint results figures for summarizing molecule counts and peptides results
######################### 

######################### 
# read in files
######################### 

# genome molecule count summary table
genome_molecule_counts <- read_tsv("results/2025-04-24-mag-results/all_molecule_counts.tsv")
genome_antismash_summary <- read_tsv("results/2025-04-24-mag-results/antismash_summary.tsv")

# genome bioactivity results
genome_bioactivity_results <- read_tsv("results/2025-04-24-mag-results/2025-04-24-mag-bioactivity-all-peptides-predictions.tsv")

# mag metadata
mag_metadata_url <- "https://raw.githubusercontent.com/MicrocosmFoods/fermentedfood_metadata_curation/refs/heads/main/data/2025-05-21-genome-metadata-food-taxonomy.tsv"

mag_metadata <- read_tsv(mag_metadata_url) %>% 
  select(-acid_type, -alcohol_level, -fermentation_temperature, -consistency, -aging_time) %>% 
  mutate(genome_name = mag_id) %>% 
  select(-mag_id) %>% 
  select(genome_name, everything())

# fermented foods peptidomics bioactivity results
proteomics_bioactivity_results <- read_tsv("results/2025-02-20-proteomics-bioactivity-results/ff_peptidomics_peptides_metadata.tsv")

# representative clusters TSV clustered @ 100% identity for genomes peptides results
genomes_peptides_100id_clusters <- read_tsv("raw_data/2025-04-23-combined-batch-results/clustering_results/clusters100_cluster.tsv", col_names = c("cluster_id", "peptide_id"))
genomes_representative_cluster_ids <- genomes_peptides_100id_clusters %>% 
  distinct(cluster_id) %>% 
  pull(cluster_id)

# representative clusters TSV clustered @ 100% identity for proteomics peptides results
proteomics_peptides_100id_clusters <- read_tsv("results/2025-02-20-proteomics-bioactivity-results/clustering/proteomics_100id_cluster.tsv", col_names = c("cluster_id", "peptide_id"))
proteomics_representative_cluster_ids <- proteomics_peptides_100id_clusters %>%
  distinct(cluster_id) %>% 
  pull(cluster_id)

######################### 
# save raw peptide results DFs for just the representative peptides -> what will go on the dashboard for viewing bioactivity
######################### 
# raw genome bioactivity df with rep seqs
rep_genome_bioactivity_results <- genome_bioactivity_results %>% 
  filter(peptide_id %in% genomes_representative_cluster_ids) %>% 
  distinct(peptide_id, .keep_all = TRUE) %>% 
  select(peptide_id, sequence, ends_with("_1")) %>% 
  rename_with(~str_remove(., "_1"), ends_with("_1"))

write_tsv(rep_genome_bioactivity_results, "results/representative_seqs_results/2025-09-08-representative-genome-peptide-seqs-bioactivity-results.tsv")

# metadata df to join with clusters
genome_bioactivity_clusters_metadata <- genome_bioactivity_results %>% 
  mutate(genome_name = str_extract(peptide_id,  "^.*?(?=_id_)")) %>% 
  left_join(mag_metadata) %>% 
  select(peptide_id, taxonomy, food_name, ingredient_group) %>% 
  left_join(genomes_peptides_100id_clusters)

genome_bioactivity_clusters_metadata_food_counts <- genome_bioactivity_clusters_metadata %>% 
  group_by(cluster_id, food_name) %>% 
  count() %>% 
  mutate(n_peptides = n) %>% 
  select(-n)

# raw proteomics bioactivity df with rep seqs
rep_proteomics_bioactivity_results <- proteomics_bioactivity_results %>% 
  filter(peptide_id %in% proteomics_representative_cluster_ids) %>% 
  distinct(peptide_id, .keep_all = TRUE) %>% 
  select(peptide_id, sequence, ends_with("_1")) %>% 
  rename_with(~str_remove(., "_1"), ends_with("_1"))

write_tsv(rep_proteomics_bioactivity_results, "results/representative_seqs_results/2025-09-08-representative-proteomics-peptide-seqs-bioactivity-results.tsv")

# join clusters counts with fermented food sample types
proteomics_bioactivity_clusters_metadata <- proteomics_bioactivity_results %>% 
  select(peptide_id, fermented_food, category) %>% 
  left_join(proteomics_peptides_100id_clusters)

proteomics_bioactivity_clusters_metadata_food_counts <- proteomics_bioactivity_clusters_metadata %>% 
  group_by(cluster_id, fermented_food) %>% 
  count() %>% 
  mutate(n_peptides = n) %>% 
  select(-n)

######################### 
# filter bioactivity results
# remove toxic peptides if have probability > 0.5
# keep a bioactivity only if probability > 0.75
######################### 

# function to filter bioactivity results
prepare_bioactivity_analysis <- function(bioactivity_data, 
                                         id_column = "peptide_id", 
                                         category_column = "food_name") {
  # Select relevant columns and pivot to long format
  bioactivity_df <- bioactivity_data %>% 
    select(!!sym(id_column), !!sym(category_column), ends_with("_1")) %>% 
    distinct() %>% 
    rename_with(~str_remove(., "_1"), ends_with("_1")) %>% 
    pivot_longer(
      -c(!!sym(id_column), !!sym(category_column)),
      names_to = "bioactivity",
      values_to = "probability"
    )
  
  # Identify toxic peptides
  toxic_peptides <- bioactivity_df %>% 
    filter(bioactivity == "TOX") %>% 
    filter(probability >= 0.5) %>% 
    pull(!!sym(id_column))
  
  # Filter out toxic peptides and add binary score
  filtered_bioactivity_df <- bioactivity_df %>% 
    filter(!.data[[id_column]] %in% toxic_peptides) %>% 
    mutate(binary_score = as.numeric(probability >= 0.75))
  
  return(list(
    bioactivity_df = bioactivity_df,
    toxic_peptides = toxic_peptides,
    filtered_bioactivity_df = filtered_bioactivity_df
  ))
}

filtered_peptidomics_bioactivity_results_binary <- prepare_bioactivity_analysis(peptidomics_bioactivity_results, category_column = "fermented_food")

filtered_proteomics_bioactivity_results <- filtered_peptidomics_bioactivity_results_binary[["filtered_bioactivity_df"]] %>% 
  filter(probability > 0.75) %>% 
  select(-binary_score)

metadata_cols <- c("genome_name", "peptide_id", "peptide_type", "taxonomy", "food_name", "ingredient_group")

genome_bioactivity_cleaned <- genome_bioactivity_results %>% 
  select(any_of(metadata_cols), ends_with("_1")) %>% 
  filter(is.na(TOX_1) | TOX_1 < 0.5) %>% 
  { 
    df <- . 
    bio_cols <- names(select(df, ends_with("_1")))
    non_tox_cols <- setdiff(bio_cols, "TOX_1")
    
    df %>%
      mutate(across(all_of(non_tox_cols),
                    ~ ifelse(.x >= 0.75, .x, NA_real_))) %>%
      rename_with(~ sub("_1$", "", .x), all_of(bio_cols))
  } %>% 
  mutate(simplified_food_name = case_when(
    grepl("cheese brine", food_name, ignore.case = TRUE) ~ food_name,
    grepl("cheese", food_name, ignore.case = TRUE) ~ "cheese",  
    grepl("tibicos fig|kefir", food_name, ignore.case = TRUE) ~ "kefir",
    TRUE ~ food_name   
  )) %>% 
  mutate(phylum = gsub(";.*", "", taxonomy)) %>% 
  select(-taxonomy, -genome_name, -food_name, -TOX) %>% 
  pivot_longer(
    !c(peptide_id, peptide_type, simplified_food_name, ingredient_group, phylum),
    names_to = "bioactivity",
    values_to = "probability"
  ) %>% 
  filter(!is.na(probability))


######################### 
# plot molecule counts from genome data
######################### 

# clean food names to simplified names for certain categories
genome_counts_metadata <- left_join(genome_molecule_counts, mag_metadata) %>% 
  select(-bgc_bacteriocin) %>% 
  mutate(simplified_food_name = case_when(
    grepl("cheese brine", food_name, ignore.case = TRUE) ~ food_name,
    grepl("cheese", food_name, ignore.case = TRUE) ~ "cheese",  
    grepl("tibicos fig|kefir", food_name, ignore.case = TRUE) ~ "kefir",
    TRUE ~ food_name   
  )) %>% 
  mutate(phylum = gsub(";.*", "", taxonomy))

genome_antismash_summary_metadata <- genome_antismash_summary %>% 
  mutate(genome_name = mag_id) %>% 
  select(-mag_id) %>% 
  left_join(mag_metadata) %>% 
  mutate(simplified_food_name = case_when(
    grepl("cheese brine", food_name, ignore.case = TRUE) ~ food_name,
    grepl("cheese", food_name, ignore.case = TRUE) ~ "cheese",  
    grepl("tibicos fig|kefir", food_name, ignore.case = TRUE) ~ "kefir",
    TRUE ~ food_name   
  )) %>%
  mutate(phylum = gsub(";.*", "", taxonomy)) %>% 
  select(genome_name, bgc_type, phylum, simplified_food_name)

# create counts of genomes in food categories
food_counts <- genome_counts_metadata %>% 
  select(genome_name, simplified_food_name) %>% 
  filter(!is.na(simplified_food_name)) %>% 
  group_by(simplified_food_name) %>% 
  count()

genome_phylum_counts <- genome_counts_metadata %>% 
  filter(!is.na(simplified_food_name)) %>% 
  filter(phylum %in% c("Actinomycetota", "Bacillota", "Pseudomonadota")) %>%
  select(genome_name, phylum) %>% 
  group_by(phylum) %>% 
  count()
  

food_order <- c("cheese", "kefir", "salami", "chocolate", "sourdough")

counts_vec <- food_counts %>%
  filter(simplified_food_name %in% c("cheese", "chocolate", "kefir", "sourdough", "salami")) %>%
  transmute(simplified_food_name, n = as.integer(n)) %>%
  deframe()

phylum_counts_vec <- genome_phylum_counts %>% 
  filter(phylum %in% c("Actinomycetota", "Bacillota", "Pseudomonadota")) %>%
  transmute(phylum, n = as.integer(n)) %>% 
  deframe()

# get top bgc categories
top_bgc_categories <- genome_antismash_summary_metadata %>% 
  count(bgc_type, sort = TRUE) %>% 
  slice_head(n=8) %>% 
  pull(bgc_type)

rename_bgc <- c(
  "Ripp-Like" = "RiPP-like",
  "Nrps"      = "NRPS",
  "Ras-Ripp"  = "RaS-RiPP"
)

# factored levels to order in the plots
bgc_levels <- sort(unique(c(top_bgc_categories_clean, "Other")))
bgc_levels <- c(setdiff(bgc_levels, "Other"), "Other")

# counts of bgcs in main foods
food_bgc_counts <- genome_antismash_summary_metadata %>%
  filter(simplified_food_name %in% food_order) %>%
  mutate(
    simplified_food_name = factor(simplified_food_name, levels = food_order),
    bgc_type = fct_other(bgc_type, keep = top_bgc_categories, other_level = "Other"),
    bgc_type = str_to_title(bgc_type),
    bgc_type = recode(bgc_type, !!!rename_bgc),
    bgc_type = factor(bgc_type, levels = bgc_levels)
  ) %>%
  count(simplified_food_name, bgc_type) %>%
  group_by(simplified_food_name)

# counts of bgcs within top phyla in the main foods
phylum_bgc_counts <- genome_antismash_summary_metadata %>%
  filter(simplified_food_name %in% food_order) %>%
  mutate(
    simplified_food_name = factor(simplified_food_name, levels = food_order),
    bgc_type = fct_other(bgc_type, keep = top_bgc_categories, other_level = "Other"),
    bgc_type = str_to_title(bgc_type),
    bgc_type = recode(bgc_type, !!!rename_bgc),
    bgc_type = factor(bgc_type, levels = bgc_levels)
  ) %>%
  filter(phylum %in% c("Actinomycetota", "Bacillota", "Pseudomonadota")) %>%
  count(simplified_food_name, phylum, bgc_type) %>%
  arrange(simplified_food_name, phylum, bgc_type) %>%
  mutate(phylum_label = paste0("italic('", phylum, "')"))

# plot main BGC groups in top food groups
bgc_counts_plot <- food_bgc_counts %>% 
  ggplot(aes(x = simplified_food_name, y = n, fill = bgc_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set3", name = "BGC Type") +
  labs(
    x = "Fermented Food",
    y = "Biosynthetic Gene Cluster Count",
    title = "Biosynthetic Gene Clusters by Fermented Food"
  ) +
  scale_x_discrete(labels = function(x) {
    base <- str_to_title(x)
    ifelse(
      x %in% names(counts_vec),
      paste0(base, " (n=", counts_vec[x], ")"),
      base
    )
  }) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic(base_size = 13) +
  theme(
    plot.title   = element_text(size = 15, face = "bold"),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    panel.grid   = element_blank()
  )

# counts of BGCs in the top foods faceted by phylum
phylum_labels <- setNames(
  paste0("italic('", names(phylum_counts_vec), "') ~ '\n(n=", phylum_counts_vec, ")'"),
  paste0("italic('", names(phylum_counts_vec), "')")
)

phylum_bgc_counts_plot <- phylum_bgc_counts %>% 
  ggplot(aes(x = simplified_food_name, y = n, fill = bgc_type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(
    ~ phylum_label,
    scales   = "free_x",
    labeller = labeller(phylum_label = phylum_labels, .default = label_parsed)
  ) +
  scale_fill_brewer(palette = "Set3", name = "BGC Type") +
  labs(
    x = "Fermented Food",
    y = "Biosynthetic Gene Cluster Count",
    title = "Biosynthetic Gene Clusters by Fermented Food"
  ) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.05))) +
  scale_x_discrete(labels = str_to_title) +
  theme_bw(base_size = 13) +
  theme(
    plot.title   = element_text(size = 15, face = "bold"),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 12, angle = 85, hjust = 1),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    panel.grid   = element_blank()
  )

phylum_bgc_inset <- phylum_bgc_counts_plot +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_blank())

# inset the phylum plot within BGC counts
bgc_with_inset <- ggdraw() +
  draw_plot(bgc_counts_plot) +
  draw_plot(phylum_bgc_inset, x = 0.25, y = 0.34, width = 0.58, height = 0.55)

bgc_with_inset

# plot peptide counts in top food groups
pretty_peptide_label <- function(x) {
  case_when(
    x == "smorf" ~ "smORF",
    x == "deeppeptide_Peptide" ~ "Cleavage Peptide",
    TRUE ~ x
  )
}

# peptide counts for top foods
peptide_counts_plot <- genome_counts_metadata %>%
  filter(simplified_food_name %in% food_order) %>%
  mutate(simplified_food_name = factor(simplified_food_name, levels = food_order)) %>%
  select(genome_name, simplified_food_name, smorf, deeppeptide_Peptide) %>%
  pivot_longer(
    !c(genome_name, simplified_food_name),
    names_to = "molecule_type",
    values_to = "count"
  ) %>%
  mutate(molecule_type = pretty_peptide_label(molecule_type)) %>%
  ggplot(aes(x = simplified_food_name, y = count, fill = molecule_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Dark2", name = "Peptide Type") +  # different palette
  labs(
    x = "Fermented Food",
    y = "Counts of Peptide Type",
    title = "Counts of Peptide Types by Fermented Food "
  ) +
  scale_x_discrete(labels = function(x) {
    base <- str_to_title(x)
    ifelse(
      x %in% names(counts_vec),
      paste0(base, " (n=", counts_vec[x], ")"),
      base
    )
  }) +
  scale_y_continuous(expand=c(0,0)) +
  theme_classic(base_size = 13) +
  theme(
    plot.title   = element_text(size = 15, face = "bold"),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    panel.grid   = element_blank()
  )

# peptide counts within top foods and phyla
phylum_peptide_counts <- genome_counts_metadata %>% 
  filter(simplified_food_name %in% food_order) %>%
  mutate(simplified_food_name = factor(simplified_food_name, levels = food_order)) %>%
  select(genome_name, simplified_food_name, phylum, deeppeptide_Peptide, smorf) %>%
  pivot_longer(
    !c(genome_name, simplified_food_name, phylum),
    names_to = "molecule_type",
    values_to = "count"
  ) %>% 
  group_by(simplified_food_name, phylum, molecule_type) %>%
  summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  arrange(simplified_food_name, phylum, molecule_type) %>% 
  filter(phylum %in% c("Actinomycetota", "Bacillota", "Pseudomonadota")) %>% 
  mutate(phylum_label = paste0("italic('", phylum, "')"))

phylum_peptide_counts_plot <- phylum_peptide_counts %>% 
  mutate(molecule_type = pretty_peptide_label(molecule_type)) %>%
  ggplot(aes(x = simplified_food_name, y = total_count, fill = molecule_type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(
    ~ phylum_label,
    scales   = "free_x",
    labeller = labeller(phylum_label = phylum_labels, .default = label_parsed)
  ) +
  scale_fill_brewer(palette = "Set2", name = "Molecule type") +  # prettier colors
  labs(
    x = "Fermented Food",
    y = "Peptide Count",
    title = "Peptide Counts by Fermented Food"
  ) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.05))) +
  scale_x_discrete(labels = str_to_title) +
  theme_bw(base_size = 13) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 12, angle = 85, hjust = 1),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    panel.grid   = element_blank()
  )

phylum_peptide_inset <- phylum_peptide_counts_plot +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_blank())

# inset the phylum plot within BGC counts
peptide_with_inset <- ggdraw() +
  draw_plot(peptide_counts_plot) +
  draw_plot(phylum_peptide_inset, x = 0.25, y = 0.34, width = 0.58, height = 0.55)

peptide_with_inset

# combine inset plots into one figure
combined_plot <- (
  bgc_with_inset / peptide_with_inset
) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave("figures/genomes-molecule-counts-summaries-plot.png", combined_plot, width=11, height=8, units=c("in"))


######################### 
# summary stats of BGC and peptide types within phyla
######################### 

genome_peptide_totals <- genome_counts_metadata %>% 
  select(genome_name, total_length, gc, phylum, species, smorf, deeppeptide_Peptide, simplified_food_name)v

bgc_counts <- genome_antismash_summary_metadata %>%
  group_by(genome_name, phylum) %>%
  summarise(
    n_bgcs = n(),
    n_bgc_types = n_distinct(bgc_type),
    .groups = "drop"
  )
  
genome_bgc_counts <- genome_counts_metadata %>% 
  select(genome_name, phylum) %>% 
  left_join(bgc_counts, by=c("genome_name", "phylum")) %>% 
  mutate(n_bgcs = replace_na(n_bgcs, 0),
         n_bgc_types = replace_na(n_bgc_types, 0)) %>%
  left_join(genome_counts_metadata) %>% 
  select(genome_name, n_bgcs, n_bgc_types)

all_totals <- left_join(genome_peptide_totals, genome_bgc_counts, by="genome_name")

molecule_group_totals_long <- all_totals %>% 
  select(-n_bgc_types) %>% 
  mutate(bgc = n_bgcs) %>% 
  mutate(cleavage_peptide = deeppeptide_Peptide) %>% 
  select(-n_bgcs, -deeppeptide_Propeptide, -deeppeptide_Peptide) %>% 
  filter(phylum %in% c("Actinomycetota", "Bacillota", "Pseudomonadota")) %>% 
  pivot_longer(
    !c(genome_name, total_length, gc, phylum, species),
    names_to = "molecule_type",
    values_to = "count"
  ) %>% 
  filter(count > 0)

plots_by_type <- molecule_group_totals_long %>%
  split(.$molecule_type) %>%
  imap(~ ggplot(.x, aes(x = total_length, y = gc)) +
         geom_point(aes(size = count, color = phylum), alpha = 0.35) + 
         scale_size_area(max_size = 10, name = "Count") +
         guides(color = "none") +  
         labs(title = .y, x = "Total length", y = "GC%") +
         theme_bw() +
         theme(legend.position = "right",
               plot.title = element_text(face = "bold"))
  )

color_legend_plot <-
  ggplot() +
  geom_point(
    data = molecule_group_totals_long,
    aes(x = 1, y = 1, color = phylum),
    size = 0, 
    alpha = 0,  
    show.legend = TRUE
  ) +
  guides(size = "none") +
  labs(color = "Phylum") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

combined <- (wrap_plots(plots_by_type, ncol = 3) | color_legend_plot) +
  plot_layout(widths = c(1, 0.18)) 

# top foods stats 
top_summary <- all_totals %>%
  filter(simplified_food_name %in% top_foods,
         phylum %in% top_phyla) %>%
  group_by(simplified_food_name, phylum) %>%
  summarise(
    n_genomes   = n(),
    total_smorf = sum(smorf, na.rm = TRUE),
    total_pep   = sum(deeppeptide_Peptide, na.rm = TRUE),
    total_bgcs  = sum(n_bgcs, na.rm = TRUE),
    total_bgc_types = sum(n_bgc_types, na.rm = TRUE),
    .groups = "drop"
  )

# phylum stats 
phylum_stats <- all_totals %>%
  filter(!is.na(phylum)) %>% 
  group_by(phylum) %>%
  summarise(
    n_genomes = n(),
    n_smorf = sum(smorf, na.rm =TRUE),
    n_pep = sum(deeppeptide_Peptide, na.rm=TRUE),
    n_bgcs = sum(n_bgcs, na.rm=TRUE),
    median_smorf = median(smorf, na.rm = TRUE),
    median_pep = median(deeppeptide_Peptide, na.rm = TRUE),
    median_bgcs = median(n_bgcs, na.rm = TRUE),
    median_bgc_types = median(n_bgc_types, na.rm = TRUE),
    .groups = "drop"
  )

# food stats
food_stats <- all_totals %>%
  filter(!is.na(simplified_food_name)) %>% 
  group_by(simplified_food_name) %>%
  summarise(
    n_genomes = n(),
    n_smorf = sum(smorf, na.rm =TRUE),
    n_pep = sum(deeppeptide_Peptide, na.rm=TRUE),
    n_bgcs = sum(n_bgcs, na.rm=TRUE),
    median_smorf = median(smorf, na.rm = TRUE),
    median_pep = median(deeppeptide_Peptide, na.rm = TRUE),
    median_bgcs = median(n_bgcs, na.rm = TRUE),
    median_bgc_types = median(n_bgc_types, na.rm = TRUE),
    .groups = "drop"
  )

######################### 
# enrichment tests of # molecule types between food groups
######################### 

# global medians of the dataset
global_medians <- all_totals %>%
  summarise(
    median_smorf = median(smorf, na.rm = TRUE),
    median_pep  = median(deeppeptide_Peptide, na.rm = TRUE),
    median_bgcs  = median(n_bgcs, na.rm = TRUE)
  )

global_means <- all_totals %>%
  summarise(
    mean_smorf = mean(smorf, na.rm = TRUE),
    mean_pep   = mean(deeppeptide_Peptide, na.rm = TRUE),
    mean_bgcs  = mean(n_bgcs, na.rm = TRUE)
  )

# food group stats
food_stats <- all_totals %>%
  group_by(simplified_food_name) %>%
  summarise(
    n_genomes   = n(),
    mean_smorf  = mean(smorf, na.rm = TRUE),
    mean_pep    = mean(deeppeptide_Peptide, na.rm = TRUE),
    mean_bgcs   = mean(n_bgcs, na.rm = TRUE),
    median_smorf = median(smorf, na.rm = TRUE),
    median_pep   = median(deeppeptide_Peptide, na.rm = TRUE),
    median_bgcs  = median(n_bgcs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    smorf_enrichment_mean  = mean_smorf  / global_means %>% pull(mean_smorf),
    pep_enrichment_mean    = mean_pep    / global_means %>% pull(mean_pep),
    bgc_enrichment_mean    = mean_bgcs   / global_means %>% pull(mean_bgcs),
    smorf_enrichment_median = median_smorf / global_medians %>% pull(median_smorf),
    pep_enrichment_median   = median_pep   / global_medians %>% pull(median_pep),
    bgc_enrichment_median   = median_bgcs  / global_medians %>% pull(median_bgcs)
  )

# food tests
food_tests <- all_totals %>%
  filter(!is.na(simplified_food_name)) %>%
  nest(data = -simplified_food_name) %>%
  mutate(
    smorf_test = map(data, ~ safe_wilcox(.x$smorf, global_smorf)),
    pep_test   = map(data, ~ safe_wilcox(.x$deeppeptide_Peptide, global_pep)),
    bgc_test   = map(data, ~ safe_wilcox(.x$n_bgcs, global_bgcs))
  ) %>%
  transmute(
    simplified_food_name,
    smorf_p = map_dbl(smorf_test, ~ .x$result$p.value %||% NA_real_),
    pep_p   = map_dbl(pep_test,   ~ .x$result$p.value %||% NA_real_),
    bgc_p   = map_dbl(bgc_test,   ~ .x$result$p.value %||% NA_real_)
  ) %>%
  mutate(
    smorf_p_adj = p.adjust(smorf_p, method = "fdr"),
    pep_p_adj   = p.adjust(pep_p,   method = "fdr"),
    bgc_p_adj   = p.adjust(bgc_p,   method = "fdr")
  )

food_enrichment <- food_stats %>%
  left_join(food_tests, by = "simplified_food_name")

plot_median_df <- food_enrichment %>%
  select(simplified_food_name,
         smorf_enrichment_median, pep_enrichment_median, bgc_enrichment_median,
         smorf_p_adj, pep_p_adj, bgc_p_adj) %>%
  pivot_longer(
    cols = -simplified_food_name,
    names_to = c("molecule", ".value"),
    names_pattern = "(.*)_(enrichment_median|p_adj)"
  ) %>% 
  filter(!is.na(simplified_food_name))


plot_mean_df <- food_enrichment %>%
  select(simplified_food_name,
         smorf_enrichment_mean, pep_enrichment_mean, bgc_enrichment_mean,
         smorf_p_adj, pep_p_adj, bgc_p_adj) %>%
  pivot_longer(
    cols = -simplified_food_name,
    names_to = c("molecule", ".value"),
    names_pattern = "(.*)_(enrichment_mean|p_adj)"
  ) %>% 
  filter(!is.na(simplified_food_name))


ggplot(plot_mean_df, aes(x = simplified_food_name, y = molecule, fill = enrichment_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = case_when(
    p_adj < 0.05 ~ "*",
    TRUE ~ ""
  )), color = "black", size = 4) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 1,
    name = "Enrichment (median ratio)"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",                    # legend at bottom
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), # vertical labels
    panel.grid = element_blank()
  ) +
  labs(
    x = "Food group",
    y = "Molecule type",
    title = "Enrichment of molecule types by food group",
    subtitle = "Color = enrichment vs global median, stars = Wilcoxon FDR significance"
  )


######################### 
# plot summaries of filtered bioactivity results
######################### 

select_bioactivities <- c("AB", "ACE", "ANIF", "AOX", "AV", "DPPIV", "IMM")
top_foods <- c("cheese", "kefir", "salami", "chocolate", "sourdough")
top_phyla <- c("Actinomycetota", "Bacillota", "Pseudomonadota")

custom_palette <- c(
  "#80B1D3", # soft blue-green
  "#66C2A5", # soft teal
  "#E6C229", # warm gold
  "#8DA0CB", # muted violet
  "#FC8D62", # light coral
  "#E78AC3", # dusty rose
  "#A6D854"  # pale olive
)

highlight_colors <- setNames(
  rep_len(custom_palette, length(select_bioactivities)),
  select_bioactivities
)

custom_colors <- c(highlight_colors, "Other" = "#D3D3D3")

bioact_labels <- c(
  "AB"    = "Antibiotic",
  "ACE"   = "ACE inhibitor",
  "ANIF"  = "Anti-inflammatory",
  "AOX"   = "Antioxidant",
  "AV"    = "Antiviral",
  "DPPIV" = "DPP IV inhibitors",
  "IMM"   = "Immunomodulatory",
  "Other" = "Other"
)

legend_breaks <- c(select_bioactivities, "Other")

bioactivity_summary_plot <- genome_bioactivity_cleaned %>%
  filter(simplified_food_name %in% top_foods) %>% 
  mutate(color_group = ifelse(bioactivity %in% select_bioactivities, bioactivity, "Other")) %>%
  group_by(simplified_food_name, color_group) %>% 
  summarise(count = n(), .groups="drop") %>%
  ggplot(aes(x = simplified_food_name, y = count, fill = color_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = custom_colors,               # ensure names match the values in color_group
    breaks = legend_breaks,
    labels = bioact_labels[legend_breaks],
    name   = "Bioactivity"
  ) +
  scale_y_continuous(
    expand = c(0,0),
    labels = scales::comma,
    breaks = seq(0, 200000, by = 25000)
  ) +
  labs(
    x = "Fermented Food",
    y = "Peptide Count",
    title = "Highlighted Bioactivities of Genome-Encoded\nPeptides in Top Fermented Foods"
  ) +
  scale_x_discrete(labels = str_to_title) +
  theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 13), 
    panel.grid   = element_blank(),
    legend.position      = c(0.6, 0.75),
    legend.justification = c("left","top"),
    legend.background    = element_rect(fill = scales::alpha("white", 0.85), color = NA)
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))

bioactivity_summary_plot
  
bioactivity_phylum_plot <- genome_bioactivity_cleaned %>%
  filter(simplified_food_name %in% top_foods) %>% 
  filter(phylum %in% top_phyla) %>% 
  mutate(color_group = ifelse(bioactivity %in% select_bioactivities, bioactivity, "Other")) %>%
  select(simplified_food_name, phylum, color_group) %>% 
  group_by(simplified_food_name, phylum, color_group) %>% 
  summarise(count = n(), .groups="drop") %>%
  mutate(phylum_label = paste0("italic('", phylum, "')")) %>% 
  ggplot(aes(x=simplified_food_name, y=count, fill = color_group)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ phylum_label, scales = "free_x", labeller = label_parsed) +
  scale_fill_manual(values = custom_colors, name = "Bioactivity") +
  scale_y_continuous(expand=c(0,0), labels = scales::comma, breaks = seq(0, 125000, by = 25000)) +
  labs(
    x = "Fermented Food",
    y = "Peptide Count",
    title = "Highlighted Bioactivities in Top Fermented Foods"
  ) +
  scale_x_discrete(labels = str_to_title) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face="bold"),
    axis.title = element_text(size = 14, face="bold"),
    panel.grid = element_blank(),
    axis.text.x  = element_text(size = 12, angle = 85, hjust = 1),
    legend.position = "none"
  )

phylum_bioactivity_inset <- bioactivity_phylum_plot +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_blank())

# inset the phylum plot within BGC counts
bioactivity_with_inset <- ggdraw() +
  draw_plot(bioactivity_summary_plot) +
  draw_plot(phylum_bioactivity_inset, x = 0.33, y = 0.32, width = 0.65, height = 0.55)

ggsave("figures/genome-bioactivity-summary-plot.png", bioactivity_with_inset, width=10, height=11, units=c("in"))

# proteomics bioactivity results

library(stringr)

proteomics_bioactivity_summary_plot <- filtered_proteomics_bioactivity_results %>% 
  filter(fermented_food != "donkey_milk") %>% 
  mutate(color_group = ifelse(bioactivity %in% select_bioactivities, bioactivity, "Other")) %>%
  select(fermented_food, color_group) %>% 
  group_by(fermented_food, color_group) %>% 
  summarise(count = n(), .groups="drop") %>%
  ggplot(aes(x = fermented_food, y = count, fill = color_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    values = custom_colors,
    breaks = legend_breaks,
    labels = bioact_labels[legend_breaks],
    name   = "Bioactivity"
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    labels = scales::comma,
    breaks = seq(0, 80000, by = 25000)
  ) +
  labs(
    x = "Fermented Food",
    y = "Peptide Count",
    title = "Highlighted Bioactivities in Fermented Foods\nfrom Proteomics Experiments"
  ) +
  # clean underscores and wrap words onto new lines
  scale_x_discrete(labels = function(x) {
    x %>%
      str_replace_all("_", " ") %>%   # replace underscores with spaces
      str_to_title() %>%              # capitalize
      str_replace_all(" ", "\n")      # force newlines at spaces
  }) +
  theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 13, lineheight = 0.9), # improve readability
    panel.grid   = element_blank(),
    legend.position      = c(0.6, 0.75),
    legend.justification = c("left", "top"),
    legend.background    = element_rect(fill = scales::alpha("white", 0.85), color = NA)
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))


proteomics_bioactivity_summary_plot

ggsave("figures/genomes-bioactivity-summary-plot.png", bioactivity_summary_plot, width=8, height=6, units=c("in"))
ggsave("figures/proteomics-bioactivity-summary-plot.png", proteomics_bioactivity_summary_plot, width=8, height=6, units=c("in"))
