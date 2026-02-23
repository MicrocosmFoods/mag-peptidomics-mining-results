library(tidyverse)

################################################## 
# Fix peptide bioactivity results files from different runs
# Update with new BLAST results to corrected database
# Remove IMM model
# Update with new ANIF model
##################################################

################################################## 
# genome peptide results
################################################## 

# original genome bioactivity results
original_genome_bioactivity_results <- read_tsv("results/2025-04-24-mag-results/2025-04-24-mag-bioactivity-all-peptides-predictions.tsv")

# drop old models and out of date BLAST results
# distinct peptide ID since won't have multiple rows per peptide for the BLAST results
modified_original_bioactivity_columns <- original_genome_bioactivity_results %>% 
  select(-IMM_1, -ANIF_1, -ANIF_2_BMdata, -sseqid, -full_sseq, -pident, -length, -qlen, -slen, -mismatch, -gapopen, -qstart, -qend, -sstart, -send, -evalue, -bitscore) %>% 
  distinct(peptide_id, .keep_all = TRUE)

# updated results with new BLAST results and models
updated_genome_bioactivity_results <- read_tsv("results/2025-04-24-mag-results/2026-02-23-updated-peptide-predictions-results.tsv")

# modify to just those columns
modified_updated_genome_bioactivity_columns <- updated_genome_bioactivity_results %>% 
  select(peptide_id, sseqid, full_sseq, pident, length, qlen, slen, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, ANIF2_20260220) %>% 
  mutate(ANIF_1 = ANIF2_20260220) %>% 
  select(-ANIF2_20260220)

# merged 
merged_genome_bioactivity_results <- left_join(modified_updated_genome_bioactivity_columns, modified_original_bioactivity_columns)

# write out merged TSV
write_tsv(merged_genome_bioactivity_results, "results/2025-04-24-mag-results/2026-02-23-updated-merged-genome-peptide-results.tsv")

################################################## 
# proteomics peptide results
################################################## 

# original proteomics bioactivity results
original_proteomics_bioactivity_results <- read_tsv("results/2025-02-20-proteomics-bioactivity-results/ff_peptidomics_peptides_metadata.tsv")

# drop old models and out of date BLAST results
# distinct peptide ID since won't have multiple rows per peptide for the BLAST results
modified_original_proteomics_columns <- original_proteomics_bioactivity_results %>% 
  select(-IMM_1, -ANIF_1, -sseqid, -full_sseq, -pident, -length, -qlen, -slen, -mismatch, -gapopen, -qstart, -qend, -sstart, -send, -evalue, -bitscore) %>% 
  distinct(peptide_id, .keep_all = TRUE)

# updated results with new BLAST results and models
updated_proteomic_bioactivity_results <- read_tsv("results/2025-02-20-proteomics-bioactivity-results/2026-02-23-ff-peptidomics-upated-run.tsv")

# modify to just the BLAST and updated model columns
modified_updated_proteomics_bioactivity_columns <- updated_proteomic_bioactivity_results %>% 
  select(peptide_id, sseqid, full_sseq, pident, length, qlen, slen, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, ANIF2_20260220) %>% 
  mutate(ANIF_1 = ANIF2_20260220) %>% 
  select(-ANIF2_20260220)

# merged
merged_proteomics_bioactivity_results <- left_join(modified_updated_proteomics_bioactivity_columns, modified_original_proteomics_columns)

# write out merged TSV
write_tsv(merged_proteomics_bioactivity_results, "results/2025-02-20-proteomics-bioactivity-results/2026-02-23-merged-proteomics-peptide-results.tsv")
