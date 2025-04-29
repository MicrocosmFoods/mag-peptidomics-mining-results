library(tidyverse)
library(RColorBrewer)
library(kableExtra)
library(UpSetR)
library(DT)

# Helper functions to be used in analysis and summary notebooks of fermented food results
# Summarizing molecule counts from MAG mining workflow runs
# Summarizing bioactive peptide results from genome-encoded peptides and peptidomics fermented food experiments

#' Read and prepare metadata and results
#' 
#' @param results_path Path to the results TSV file
#' @param metadata_path Path to the metadata TSV file
#' @param genome_id_col Name of the column in metadata containing genome identifiers
#' @param category_col Name of the column containing substrate/food categories
#' @return A dataframe with combined results and metadata
prepare_data <- function(results_path, metadata_path, genome_id_col = "genome_name", category_col = "substrate_category") {
  # Read data
  results <- read_tsv(results_path)
  metadata <- read_tsv(metadata_path)
  
  # Ensure genome_name column exists in metadata
  if (genome_id_col != "genome_name") {
    metadata <- metadata %>%
      mutate(genome_name = !!sym(genome_id_col)) %>%
      select(genome_name, everything())
  }
  
  # Join results with metadata
  combined_data <- left_join(results, metadata) %>%
    filter(!is.na(!!sym(category_col)))
  
  return(combined_data)
}

#' Generate summary statistics by category
#' 
#' @param data Combined data from prepare_data function
#' @param category_col Name of the column to group by (e.g., substrate_category)
#' @return A dataframe with summary statistics by category
generate_summary_stats <- function(data, category_col = "substrate_category") {
  summary_stats <- data %>% 
    group_by(!!sym(category_col)) %>% 
    summarize(
      `Total Genomes` = n(),
      `Avg Completeness` = mean(completeness),
      `Avg Contamination` = mean(contamination),
      `Median Contig Size` = median(contigs),
      `Min Contig Size` = min(contigs),
      `Max Contig Size` = max(contigs),
      `Total smORFs` = sum(smorf),
      `Total Cleavage Peptides` = sum(deeppeptide_Peptide),
      `Total Cleavage Propeptides` = sum(deeppeptide_Propeptide),
      `Total BGCs` = sum(across(starts_with("bgc_"))),
    ) %>% 
    arrange(desc(`Total Genomes`)) %>% 
    mutate_if(is.numeric, round, 2) %>%
    mutate(
      `Median Contig Size` = as.integer(`Median Contig Size`)) %>% 
    mutate(`Category` = !!sym(category_col)) %>% 
    select(-!!sym(category_col)) %>% 
    select(Category, everything())
  
  return(summary_stats)
}

#' Filter data for high-quality genomes
#' 
#' @param data Combined data from prepare_data function
#' @param min_completeness Minimum completeness percentage (default: 90)
#' @param max_contamination Maximum contamination percentage (default: 10)
#' @param max_contigs Maximum number of contigs (default: 200)
#' @return Filtered dataframe with only high-quality genomes
filter_high_quality <- function(data, min_completeness = 90, max_contamination = 10, max_contigs = 200) {
  filtered_data <- data %>% 
    filter(completeness > min_completeness & 
             contamination < max_contamination & 
             contigs < max_contigs)
  
  return(filtered_data)
}

#' Create genome quality scatter plot
#' 
#' @param data Data to plot
#' @param category_col Name of the column for coloring points (default: "substrate_category")
#' @param title Plot title
#' @return A ggplot object
plot_genome_quality <- function(data, category_col = "substrate_category", title = "Genome Quality Statistics") {
  # Combine colors from multiple palettes
  custom_colors <- c(
    brewer.pal(8, "Set2"),
    brewer.pal(6, "Dark2")
  )
  
  p <- data %>% 
    ggplot(aes(x = completeness, y = contamination, size = contigs)) +
    geom_point(aes(color = !!sym(category_col))) +
    theme_classic() +
    theme(legend.position = "right") + 
    scale_color_manual(values = custom_colors) +
    ggtitle(title)
  
  return(p)
}

#' Create BGC distribution plot
#' 
#' @param data Data to plot (should be filtered for high quality)
#' @param category_col Name of the category column (default: "substrate_category")
#' @param include_peptides Whether to include smORF and cleavage peptides (default: FALSE)
#' @return A ggplot object
plot_molecule_distribution <- function(data, category_col = "substrate_category", include_peptides = FALSE) {
  # Get molecule columns
  molecule_columns <- c(
    names(data)[grep("^bgc_", names(data))]
  )
  
  if (include_peptides) {
    molecule_columns <- c(molecule_columns, "smorf", "deeppeptide_Peptide")
  }
  
  # Get genome counts
  genome_counts <- data %>%
    count(!!sym(category_col)) %>%
    rename(genome_count = n)
  
  # Prepare counts of each molecule
  molecule_counts <- data %>%
    group_by(!!sym(category_col)) %>%
    summarise(across(all_of(molecule_columns), sum)) %>%
    pivot_longer(
      cols = all_of(molecule_columns),
      names_to = "molecule_type",
      values_to = "count"
    ) %>%
    mutate(
      molecule_type = case_when(
        startsWith(molecule_type, "bgc_") ~ str_replace(molecule_type, "bgc_", ""),
        molecule_type == "smorf" ~ "smORF",
        molecule_type == "deeppeptide_Peptide" ~ "Cleavage_Peptides",
        TRUE ~ molecule_type
      )
    )
  
  molecule_counts <- molecule_counts %>%
    left_join(genome_counts, by = category_col) %>%
    mutate(category_label = paste0(!!sym(category_col), " (n=", genome_count, ")")) %>%
    mutate(
      molecule_type = case_when(
        molecule_type == "other" ~ "Other BGC",
        TRUE ~ molecule_type
      )
    )
  
  # Filter out peptides if not included
  if (!include_peptides) {
    molecule_counts <- molecule_counts %>%
      filter(molecule_type != 'smORF' & molecule_type != 'Cleavage_Peptides')
  }
  
  # Create plot
  p <- molecule_counts %>%
    ggplot(aes(x = category_label, y = count, fill = molecule_type)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    labs(
      x = "Category",
      y = "Count",
      fill = "Molecule Type",
      title = "Distribution of BGC Types by Category"
    ) +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(expand = c(0, 0))
  
  return(p)
}

#' Format summary table with kableExtra
#' 
#' @param summary_data Summary data from generate_summary_stats
#' @return A formatted kable table
format_summary_table <- function(summary_data) {
  formatted_table <- summary_data %>%
    kbl() %>%
    kable_styling(bootstrap_options = c("striped", "hover"),
                  full_width = FALSE) %>%
    row_spec(0, bold = TRUE)
  
  return(formatted_table)
}

#' Process peptide bioactivity data and extract peptide types
#' 
#' @param bioactivity_results Path to bioactivity results file
#' @param deeppeptide_results Path to DeepPeptide results file
#' @param smorfinder_results Path to SMORFinder results file
#' @param ripp_results Path to RiPP results file
#' @return A dataframe with peptide bioactivity data and peptide types
process_peptide_bioactivity_info <- function(bioactivity_results, 
                                        deeppeptide_results = NULL, 
                                        smorfinder_results = NULL, 
                                        ripp_results = NULL) {
  # Read bioactivity results
  bioactivity_data <- bioactivity_results
  
  # Process peptide type information if provided
  if (!is.null(deeppeptide_results) && !is.null(smorfinder_results) && !is.null(ripp_results)) {
    # Read in info for each peptide type
    deeppeptide_info <- read_tsv(deeppeptide_results) %>% 
      mutate(peptide_type = paste0(peptide_type, "_", tolower(peptide_class))) %>% 
      mutate(peptide_id = paste0(genome_name, "_id_", peptide_id)) %>% 
      select(peptide_id, peptide_type)
    
    smorfinder_info <- read_tsv(smorfinder_results) %>% 
      mutate(peptide_id = paste0(genome_name, "_id_", seqid)) %>% 
      mutate(peptide_type = "smorf") %>% 
      select(peptide_id, peptide_type)
    
    ripp_info <- read_tsv(ripp_results) %>% 
      mutate(bgc_name = gsub(".gbk", "", bgc_file)) %>% 
      mutate(peptide_id = paste0(mag_id, "_id_", scaffold_name, "_", bgc_name)) %>% 
      mutate(peptide_type = "ripp") %>% 
      select(peptide_id, peptide_type)
    
    # Combine all peptide types
    all_peptide_types_info <- rbind(deeppeptide_info, smorfinder_info, ripp_info)
    
    # Clean peptide IDs and join with peptide type info
    bioactivity_data <- bioactivity_data %>% 
      mutate(peptide_id = gsub("\\.gbk.*", "", peptide_id)) %>% 
      left_join(all_peptide_types_info) %>% 
      filter(peptide_type != "cleavage_propeptide")
  }
  
  return(bioactivity_data)
}

#' Prepare bioactivity data for analysis, remove toxic peptide predictions from dataset
#' 
#' @param bioactivity_data Bioactivity data from process_peptide_bioactivity
#' @param id_column Column name containing identifiers
#' @param category_column Column name containing categories
#' @return A dataframe in long format with bioactivity probabilities
prepare_bioactivity_analysis <- function(bioactivity_data, 
                                         id_column = "peptide_id", 
                                         category_column = "substrate_category") {
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
    mutate(binary_score = as.numeric(probability >= 0.5))
  
  return(list(
    bioactivity_df = bioactivity_df,
    toxic_peptides = toxic_peptides,
    filtered_bioactivity_df = filtered_bioactivity_df
  ))
}

#' Generate summary statistics for bioactivity data
#' 
#' @param bioactivity_analysis Output from prepare_bioactivity_analysis
#' @param id_column Column name containing identifiers
#' @param category_column Column name containing categories
#' @param genome_id_column Column name containing genome identifiers (optional)
#' @return A dataframe with summary statistics
generate_bioactivity_summary <- function(bioactivity_analysis, 
                                         id_column = "peptide_id", 
                                         category_column = "substrate_category",
                                         genome_id_column = NULL) {
  
  bioactivity_df <- bioactivity_analysis$bioactivity_df
  toxic_peptides <- bioactivity_analysis$toxic_peptides
  filtered_bioactivity_df <- bioactivity_analysis$filtered_bioactivity_df
  
  # Base summary
  if (is.null(genome_id_column)) {
    summary_stats <- bioactivity_df %>%
      select(!!sym(id_column), !!sym(category_column)) %>%
      distinct() %>%
      group_by(!!sym(category_column)) %>%
      summarise(
        total_peptides = n(),
        non_toxic_peptides = sum(!.data[[id_column]] %in% toxic_peptides)
      )
  } else {
    summary_stats <- bioactivity_df %>%
      select(!!sym(id_column), !!sym(category_column), !!sym(genome_id_column)) %>%
      distinct() %>%
      group_by(!!sym(category_column)) %>%
      summarise(
        total_genomes = n_distinct(!!sym(genome_id_column)),
        total_peptides = n(),
        non_toxic_peptides = sum(!.data[[id_column]] %in% toxic_peptides)
      )
  }
  
  # Calculate confident bioactivity count
  confident_bioactivity <- filtered_bioactivity_df %>%
    filter(probability == 1) %>%
    select(!!sym(id_column), !!sym(category_column)) %>%
    distinct() %>%
    group_by(!!sym(category_column)) %>%
    summarise(confident_bioactivity_count = n())
  
  # Calculate most common bioactivity
  most_common_bioactivity <- filtered_bioactivity_df %>%
    filter(binary_score == 1) %>%
    group_by(!!sym(category_column), bioactivity) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(!!sym(category_column)) %>%
    slice_max(count, n = 1) %>%
    slice_head(n = 1) %>%
    select(!!sym(category_column), 
           most_common_bioactivity = bioactivity, 
           most_common_bioactivity_count = count)
  
  # Combine all summaries
  summary_stats <- summary_stats %>%
    left_join(confident_bioactivity, by = category_column) %>%
    left_join(most_common_bioactivity, by = category_column) %>%
    arrange(desc(total_peptides))
  
  return(summary_stats)
}

#' Plot bioactivity distribution by category
#' 
#' @param filtered_bioactivity_df Filtered bioactivity dataframe
#' @param category_column Column name containing categories
#' @param probability_threshold Minimum probability threshold (default: 0.75)
#' @param title Plot title
#' @return A ggplot object
plot_bioactivity_distribution <- function(filtered_bioactivity_df, 
                                          category_column = "substrate_category",
                                          probability_threshold = 0.75,
                                          title = "Absolute Counts of Bioactivity Labels") {
  
  # Create plot
  p <- filtered_bioactivity_df %>% 
    filter(probability >= probability_threshold) %>% 
    ggplot(aes(x = bioactivity)) +
    geom_bar(aes(fill = !!sym(category_column))) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      x = "Bioactivity",
      y = "Count",
      fill = "Category",
      title = title
    )
  
  # Use appropriate color palette based on number of categories
  n_categories <- length(unique(filtered_bioactivity_df[[category_column]]))
  
  if (n_categories <= 12) {
    p <- p + scale_fill_brewer(palette = "Paired")
  } else {
    # For more than 12 categories, use a combination of palettes
    brewer_palette <- brewer.pal(n = 12, name = "Paired")
    extra_colors <- colorRampPalette(brewer_palette)(n_categories)
    p <- p + scale_fill_manual(values = extra_colors)
  }
  
  return(p)
}

#' Plot bioactivity distribution by peptide type
#' 
#' @param filtered_bioactivity_df Filtered bioactivity dataframe with peptide_type column
#' @param probability_threshold Minimum probability threshold (default: 0.75)
#' @param title Plot title
#' @return A ggplot object
plot_peptide_type_bioactivity <- function(filtered_bioactivity_df,
                                          probability_threshold = 0.75,
                                          title = "Bioactivity Distribution by Peptide Type") {
  
  # Check if peptide_type column exists
  if (!"peptide_type" %in% colnames(filtered_bioactivity_df)) {
    stop("peptide_type column not found in the dataframe")
  }
  
  # Create plot
  p <- filtered_bioactivity_df %>% 
    filter(probability >= probability_threshold) %>% 
    ggplot(aes(x = bioactivity)) +
    geom_bar(aes(fill = peptide_type)) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      x = "Bioactivity",
      y = "Count",
      fill = "Peptide Type",
      title = title
    )
  
  return(p)
}

#' Process sequence clustering results
#' 
#' @param cluster_files Vector of paths to cluster files
#' @param pattern Pattern to extract sample name from file name
#' @param remove_pattern Pattern to remove from sample name
#' @return A dataframe with clustering information
process_sequence_clusters <- function(cluster_files, 
                                      pattern = "*_cluster.tsv",
                                      remove_pattern = "_clustering_cluster.tsv$") {
  
  all_clusters <- map_df(cluster_files, ~{
    sample_name <- basename(.x) %>%
      str_remove(remove_pattern)
    
    read_tsv(.x, col_names = c("cluster_name", "peptide_id")) %>% 
      mutate(sample_name = gsub("_peptides", "", sample_name))
  })
  
  return(all_clusters)
}

#' Plot cluster size distribution
#' 
#' @param cluster_data Output from process_sequence_clusters
#' @param min_total Minimum total sequences to include a sample
#' @param title Plot title
#' @return A ggplot object
plot_cluster_sizes <- function(cluster_data, 
                               min_total = 100,
                               title = "Proportion of Cluster Sizes within Samples") {
  
  # Calculate sample totals
  sample_totals <- cluster_data %>%
    group_by(sample_name) %>%
    summarize(total = n()) %>%
    filter(total > min_total) %>% 
    mutate(sample = paste0(sample_name, "\n(n=", total, ")"))
  
  # Calculate cluster counts
  cluster_counts <- cluster_data %>% 
    group_by(cluster_name) %>% 
    count() %>% 
    mutate(protein_count = n) %>% 
    select(cluster_name, protein_count)
  
  # Create plot
  p <- cluster_data %>%
    # Join with the totals to get the new labels
    left_join(sample_totals, by = "sample_name") %>%
    left_join(cluster_counts, by = "cluster_name") %>% 
    drop_na() %>% 
    mutate(size_category = case_when(
      protein_count == 1 ~ "Singleton",
      protein_count <= 5 ~ "Small (2-5)",
      protein_count <= 10 ~ "Medium (6-10)",
      TRUE ~ "Large (>10)"
    )) %>%
    mutate(size_category = factor(size_category, 
                                  levels = c("Singleton", "Small (2-5)", "Medium (6-10)", "Large (>10)")
    )) %>%
    ggplot(aes(y = sample, fill = size_category)) +
    geom_bar(position = "fill") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.9),
          panel.grid = element_blank(),
          legend.position = "bottom") +
    labs(
      title = title,
      x = "Proportion",
      y = "Sample",
      fill = "Cluster Size"
    ) +
    scale_fill_brewer(palette = 7)
  
  return(p)
}

#' Process Peptipedia BLAST hits
#' 
#' @param bioactivity_data Bioactivity data with BLAST hit information
#' @param peptipedia_metadata Peptipedia metadata
#' @param category_column Column name containing categories
#' @return A dataframe with processed BLAST hits
process_peptipedia_hits <- function(bioactivity_data, 
                                    peptipedia_metadata,
                                    category_column = "substrate_category") {
  
  # Extract BLAST hits
  blastp_hits <- bioactivity_data %>% 
    filter(!is.na(sseqid)) %>% 
    mutate(peptipedia_peptide_id = sseqid) %>% 
    select(peptide_id, sequence, peptipedia_peptide_id, pident, length, 
           qlen, slen, evalue, !!sym(category_column))
  
  # Join with Peptipedia metadata
  blastp_hits_metadata <- blastp_hits %>% 
    left_join(peptipedia_metadata)
  
  # Calculate summary statistics
  stats_summary <- blastp_hits_metadata %>%
    group_by(peptipedia_peptide_id) %>%
    summarise(
      n = n(),
      median_pident = median(pident, na.rm = TRUE),
      median_length = median(length, na.rm = TRUE)
    )
  
  # Select bioactivity columns
  bioactivity_columns <- blastp_hits_metadata %>%
    select(peptipedia_peptide_id, starts_with("aceinhibitors"):ends_with("fermfoodb")) %>%
    distinct()
  
  # Combine summary stats with bioactivity columns
  full_summary <- stats_summary %>%
    left_join(bioactivity_columns, by = "peptipedia_peptide_id")
  
  return(list(
    blastp_hits_metadata = blastp_hits_metadata,
    stats_summary = full_summary
  ))
}

#' Plot Peptipedia hit distribution
#' 
#' @param blastp_hits_metadata BLAST hits metadata from process_peptipedia_hits
#' @param stats_summary Summary statistics from process_peptipedia_hits
#' @param category_column Column name containing categories
#' @param min_hits Minimum number of hits to include a peptide ID
#' @param fermfoodb_only Whether to only include FermFooDB peptides
#' @param title Plot title
#' @return A ggplot object
plot_peptipedia_hits <- function(blastp_hits_metadata, 
                                 stats_summary,
                                 category_column = "substrate_category",
                                 min_hits = 10,
                                 fermfoodb_only = FALSE,
                                 title = "Distribution of Peptipedia Hits") {
  
  # Filter peptide IDs based on criteria
  if (fermfoodb_only) {
    peptide_ids <- blastp_hits_metadata %>% 
      filter(fermfoodb == TRUE) %>% 
      pull(peptipedia_peptide_id)
  } else {
    peptide_ids <- stats_summary %>% 
      filter(n > min_hits) %>% 
      pull(peptipedia_peptide_id)
  }
  
  # Create plot
  p <- blastp_hits_metadata %>% 
    filter(peptipedia_peptide_id %in% peptide_ids) %>%
    ggplot(aes(x = as.character(peptipedia_peptide_id))) +
    geom_bar(aes(fill = !!sym(category_column))) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      x = "Peptipedia Peptide ID",
      y = "Count",
      fill = "Category",
      title = title
    )
  
  return(p)
}

#' Format Peptipedia hits as interactive table
#' 
#' @param stats_summary Summary statistics from process_peptipedia_hits
#' @param fermfoodb_ids Vector of FermFooDB peptide IDs (optional)
#' @param title Table caption
#' @return A DT datatable object
format_peptipedia_table <- function(stats_summary, 
                                    fermfoodb_ids = NULL,
                                    title = "Summary of Peptipedia Diamond Blastp Hits") {
  
  # Filter for FermFooDB IDs if provided
  if (!is.null(fermfoodb_ids)) {
    table_data <- stats_summary %>%
      filter(peptipedia_peptide_id %in% fermfoodb_ids)
  } else {
    table_data <- stats_summary
  }
  
  # Create datatable
  dt <- table_data %>%
    datatable(
      colnames = c(
        "Peptipedia Peptide ID",
        "Hit Counts",
        "Median % identity of Hit",
        "Median Aln Length of Hit",
        "Ace-Inhibitor",
        "Anticancer",
        "Anti-hypertensive",
        "Antiinflammatory",
        "Antimicrobial",
        "Antioxidative",
        "Cytotoxic",
        "Hormonal", 
        "Immunomodulatory",
        "Metabolic",
        "FermFooDB"
      ),
      caption = title,
      options = list(
        pageLength = 20,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        order = list(list(1, 'desc')),
        autoWidth = TRUE
      ),
      rownames = FALSE,
      filter = 'top',
      extensions = c('Buttons', 'FixedHeader'),
      class = 'cell-border stripe'
    )
  
  return(dt)
}
