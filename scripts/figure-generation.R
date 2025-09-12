library(tidyverse)
library(ggh4x)
library(RColorBrewer)
library(patchwork)
library(cowplot)
library(UpSetR)

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
genomes_peptides_100id_clusters <- read_tsv("raw_data/2025-04-23-combined-batch-results/clusters100/clusters.tsv", col_names = c("cluster_id", "peptide_id"))
genomes_representative_cluster_ids <- genomes_peptides_100id_clusters %>% 
  distinct(cluster_id) %>% 
  pull(cluster_id)

# representative clusters TSV clustered @ 100% identity for proteomics peptides results
proteomics_peptides_100id_clusters <- read_tsv("results/2025-02-20-proteomics-bioactivity-results/clusters100/clusters.tsv", col_names = c("cluster_id", "peptide_id"))
proteomics_representative_cluster_ids <- proteomics_peptides_100id_clusters %>%
  distinct(cluster_id) %>% 
  pull(cluster_id)

# peptipedia DB curated metadata
peptipedia_metadata <- read_tsv("metadata/2025-08-11-peptipedia-validated-metadata.tsv") %>% 
  mutate(peptipedia_id = peptide_id) %>% 
  select(-peptide_id) %>% 
  mutate(peptipedia_sequence = sequence) %>% 
  select(-sequence)

######################### 
# save raw peptide results DFs for just the representative peptides
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
# peptide bioactivity results
# Analyze the representative peptides (clustered at 100% sequence identity)
######################### 

# remove toxic peptides and pivot to long for counts
prep_nr_long <- function(df = rep_genome_bioactivity_results,
                         id_cols = c("peptide_id","sequence"),
                         tox_col = "TOX", tox_cut = 0.5) {
  bio_cols <- df |> select(where(is.numeric)) |> names()
  bio_cols <- setdiff(bio_cols, id_cols)
  
  if (tox_col %in% names(df)) {
    df <- df |> filter(.data[[tox_col]] < tox_cut | is.na(.data[[tox_col]]))
    bio_cols <- setdiff(bio_cols, tox_col)
  }
  
  long <- df |>
    select(any_of(id_cols), any_of(bio_cols)) |>
    pivot_longer(cols = any_of(bio_cols),
                 names_to = "bioactivity",
                 values_to = "probability")
  
  long
}

# plot bioactivity results
# create bar plot of number of peptides with a bioactivity given probability threshold
# plot UpSet plot of overlapping categories
# create supplementary bar plot figure at different thresholds 0.5 > 1 at 0.10 intervals
plot_bioactivity_suite <- function(long_df,
                                   id_col = "peptide_id",
                                   main_thr = 0.75,
                                   save_dir = "figures",
                                   save = TRUE) {
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Main bar
  bar_main_df <- long_df |>
    summarise(count = sum(probability >= main_thr, na.rm = TRUE), .by = bioactivity) |>
    arrange(desc(count)) |>
    mutate(bioactivity = forcats::fct_reorder(bioactivity, count))
  
  p_bar <- ggplot(bar_main_df, aes(bioactivity, count)) +
    geom_col() +
    labs(x = "Bioactivity", y = "Number of Peptides with Predicted Bioactivity",
         title = sprintf("Non-redundant peptides per bioactivity (TOX≥0.5 removed, p≥%.2f)", main_thr)) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    scale_y_continuous(expand=c(0,0))
  
  if (save) ggsave(file.path(save_dir, "nrpeptides_main_bar_p075.png"), p_bar,
                   width = 8, height = 5, units = "in", dpi = 300)
  
  # Prep UpSet data only (no plotting here)
  upset_data <- NULL
  wide_bin <- long_df |>
    mutate(hit = as.integer(probability >= main_thr)) |>
    select(all_of(id_col), bioactivity, hit) |>
    distinct() |>
    tidyr::pivot_wider(names_from = bioactivity, values_from = hit, values_fill = 0)
  
  cols <- setdiff(names(wide_bin), id_col)
  keep <- cols[colSums(as.data.frame(wide_bin[, cols, drop = FALSE])) > 0]
  if (length(keep)) {
    up_df <- as.data.frame(wide_bin[, keep, drop = FALSE])
    for (nm in names(up_df)) up_df[[nm]] <- as.integer(up_df[[nm]] != 0)
    upset_data <- up_df
  }
  
  # Supplemental bars across thresholds
  thresholds <- seq(0.5, 1.0, by = 0.1)
  supp_df <- lapply(thresholds, function(t) {
    long_df |>
      summarise(count = sum(probability >= t, na.rm = TRUE), .by = bioactivity) |>
      mutate(threshold = sprintf("%.2f", t))
  }) |>
    dplyr::bind_rows() |>
    mutate(threshold = factor(threshold, levels = sprintf("%.2f", thresholds)))
  
  p_supp <- ggplot(supp_df, aes(bioactivity, count)) +
    geom_col() +
    facet_wrap(~ threshold, ncol = 3, scales = "free_y") +
    labs(x = "Bioactivity", y = "Peptide count",
         title = "Non-redundant peptide counts by bioactivity across cutoffs (TOX≥0.5 removed)") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          strip.text = element_text(face = "bold")) +
    scale_y_continuous(expand=c(0,0))
  
  list(
    bar_plot   = p_bar,
    bar_data   = bar_main_df,
    upset_data = upset_data,  # ready for UpSetR outside the function
    wide_hits  = wide_bin,
    p_supp_data = supp_df,
    p_supp_plot = p_supp
  )
}


# run results for the representative peptides dataframe - using main default probability threshold of 0.75, and for supplemental figures at thresholds ranging from 0.5 to 1.0 and .10 increments

# representative peptides from genomes dataset
long_nr_genomes <- prep_nr_long(rep_genome_bioactivity_results) 
genomes_peptides_results <- plot_bioactivity_suite(long_nr_genomes, main_thr = 0.90)
genomes_peptides_results[["bar_plot"]]
genomes_peptides_results[["p_supp_plot"]]

# upset plot
keep_sets <- colnames(genomes_peptides_results$upset_data)[
  colSums(genomes_peptides_results$upset_data) > 1500
]

# run UpSetR only on those
UpSetR::upset(
  genomes_peptides_results$upset_data,
  sets = keep_sets,
  order.by = "freq"
)

# representative peptides from proteomics dataset
long_nr_proteomics <- prep_nr_long(rep_proteomics_bioactivity_results)
proteomics_peptides_results <- plot_bioactivity_suite(long_nr_proteomics, main_thr=0.90)

proteomics_peptides_results[["bar_plot"]]
# upset plot
keep_sets <- colnames(proteomics_peptides_results$upset_data)[
  colSums(proteomics_peptides_results$upset_data) > 85
]

UpSetR::upset(
  proteomics_peptides_results$upset_data,
  sets = keep_sets,
  order.by = "freq"
)

######################### 
# peptide BLAST results
# analyze peptides with a BLAST hit to FermFooDB/Peptipedia, how many times it's detected in the entire dataset
######################### 

representative_genome_seqs_blast_results <- genome_bioactivity_results %>% 
  filter(peptide_id %in% genomes_representative_cluster_ids) %>% 
  filter(!is.na(sseqid)) %>% 
  mutate(peptipedia_id = sseqid) %>% 
  select(peptide_id, sequence, peptipedia_id, full_sseq, pident, length, qlen, slen, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) %>% 
  left_join(peptipedia_metadata, by = "peptipedia_id")

representative_proteomics_seqs_blast_results <- proteomics_bioactivity_results %>% 
  filter(peptide_id %in% proteomics_representative_cluster_ids) %>% 
  filter(!is.na(sseqid)) %>% 
  mutate(peptipedia_id = sseqid) %>% 
  select(peptide_id, sequence, peptipedia_id, full_sseq, pident, length, qlen, slen, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) %>% 
  left_join(peptipedia_metadata, by="peptipedia_id")
  

