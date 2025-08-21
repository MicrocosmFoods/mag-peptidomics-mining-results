library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(cowplot)

######################### 
# vignette figures for counts of molecules and bioactivity results
######################### 

######################### 
# read in files
######################### 

# genome molecule count summary table
genome_molecule_counts <- read_tsv("results/2025-04-24-mag-results/all_molecule_counts.tsv")

# genome bioactivity results
genome_bioactivity_results <- read_tsv("results/cleaned_results/2025-08-05-mag-bioactivity-info.tsv.gz")

# mag metadata
mag_metadata_url <- "https://raw.githubusercontent.com/MicrocosmFoods/fermentedfood_metadata_curation/refs/heads/main/data/2025-05-21-genome-metadata-food-taxonomy.tsv"

mag_metadata <- read_tsv(mag_metadata_url) %>% 
  select(-acid_type, -alcohol_level, -fermentation_temperature, -consistency, -aging_time) %>% 
  mutate(genome_name = mag_id) %>% 
  select(-mag_id) %>% 
  select(genome_name, everything())

# fermented foods peptidomics bioactivity results
peptidomics_bioactivity_results <- read_tsv("results/2025-02-20-proteomics-bioactivity-results/ff_peptidomics_peptides_metadata.tsv")

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

filtered_genome_bioactivity_results <- prepare_bioactivity_analysis(genome_bioactivity_results)

filtered_peptidomics_bioactivity_results <- prepare_bioactivity_analysis(peptidomics_bioactivity_results, category_column = "fermented_food")

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

# create counts of genomes in food categories
food_counts <- genome_counts_metadata %>% 
  select(genome_name, simplified_food_name) %>% 
  filter(!is.na(simplified_food_name)) %>% 
  group_by(simplified_food_name) %>% 
  count()

food_order <- c("cheese", "kefir", "salami", "chocolate", "sourdough")

counts_vec <- food_counts %>%
  filter(simplified_food_name %in% c("cheese", "chocolate", "kefir", "sourdough", "salami")) %>%
  transmute(simplified_food_name, n = as.integer(n)) %>%
  deframe()

pretty_bgc_label <- function(x) {
  x <- gsub("^bgc_", "", x)  
  keep <- c("NRPS", "T1PKS", "T3PKS", "RiPP-like")
  out <- ifelse(x %in% keep, x, str_to_sentence(x))
  out
}

# plot main BGC groups in top food groups
bgc_counts_plot <- genome_counts_metadata %>%
  filter(simplified_food_name %in% food_order) %>%
  mutate(simplified_food_name = factor(simplified_food_name, levels = food_order)) %>%
  select(genome_name, simplified_food_name, bgc_NRPS, `bgc_RiPP-like`, bgc_T1PKS,
         bgc_T3PKS, bgc_betalactone, bgc_other, bgc_terpene) %>%
  pivot_longer(
    !c(genome_name, simplified_food_name),
    names_to = "molecule_type",
    values_to = "count"
  ) %>%
  mutate(molecule_type = pretty_bgc_label(molecule_type)) %>%
  ggplot(aes(x = simplified_food_name, y = count, fill = molecule_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set2", name = "BGC Type") +  # prettier colors
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
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    panel.grid   = element_blank()
  )

# counts of BGCs in the top foods faceted by phylum
phylum_bgc_counts <- genome_counts_metadata %>% 
  filter(simplified_food_name %in% food_order) %>%
  mutate(simplified_food_name = factor(simplified_food_name, levels = food_order)) %>%
  select(genome_name, simplified_food_name, phylum, bgc_NRPS, `bgc_RiPP-like`, bgc_T1PKS,
         bgc_T3PKS, bgc_betalactone, bgc_other, bgc_terpene) %>%
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
  
phylum_bgc_counts_plot <- phylum_bgc_counts %>% 
  mutate(molecule_type = pretty_bgc_label(molecule_type)) %>%
  ggplot(aes(x = simplified_food_name, y = total_count, fill = molecule_type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ phylum_label, scales = "free_x", labeller = label_parsed) +
  scale_fill_brewer(palette = "Set2", name = "BGC Type") +  # prettier colors
  labs(
    x = "Fermented Food",
    y = "Biosynthetic Gene Cluster Count",
    title = "Biosynthetic Gene Clusters by Fermented Food"
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

phylum_bgc_inset <- phylum_bgc_counts_plot +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_blank())

# inset the phylum plot within BGC counts
bgc_with_inset <- ggdraw() +
  draw_plot(bgc_counts_plot) +
  draw_plot(phylum_bgc_inset, x = 0.25, y = 0.35, width = 0.60, height = 0.55)

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
    y = "Peptide Count",
    title = "Peptide Counts by Fermented Food"
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
    plot.title   = element_text(size = 18, face = "bold"),
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
  facet_wrap(~ phylum_label, scales = "free_x", labeller = label_parsed) +
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
  draw_plot(phylum_peptide_inset, x = 0.25, y = 0.35, width = 0.60, height = 0.55)

peptide_with_inset

# combine inset plots into one figure
combined_plot <- bgc_with_inset / peptide_with_inset +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

combined_plot
