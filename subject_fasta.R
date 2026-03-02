library(rentrez) 
# Get all subject files from Example_Files 
subject_files <- list.files("Example_Files", pattern = "^subject_.*\\.txt$", full.names = TRUE)

# Function to clean accession IDs (remove trailing | etc.)
clean_accession <- function(x) {
  # Remove any database prefix and trailing |
  gsub("^[a-z]+\\|([A-Z0-9_.]+)\\|$", "\\1", x)
}

# Function to append [Organism] to SwissProt headers
add_organism_to_swissprot <- function(fasta_text, accs) {
  # Split sequences by lines starting with ">"
  seqs <- unlist(strsplit(fasta_text, "\n(?=>)", perl = TRUE))
  seqs <- seqs[seqs != ""]  # remove empty
  
  for (i in seq_along(seqs)) {
    lines <- unlist(strsplit(seqs[i], "\n"))
    
    # Make sure header has leading ">"
    header_idx <- which(grepl("^>", lines))[1]
    if (is.na(header_idx)) next
    header <- lines[header_idx]
    
    # Only SwissProt sequences
    if (!grepl("^>sp\\|", header)) next
    
    # Extract cleaned accession
    acc <- gsub("^>sp\\|([^|]+)\\|.*", "\\1", header)
    
    # Fetch organism
    summary <- tryCatch(entrez_summary(db="protein", id=acc), error=function(e) NULL)
    if (is.null(summary) || is.null(summary$organism)) next
    org <- summary$organism
    
    # Append [Organism] to header
    lines[header_idx] <- paste0(header, " [", org, "]")
    
    # Recombine sequence
    seqs[i] <- paste(lines, collapse="\n")
    Sys.sleep(0.3)
  }
  
  # Recombine all sequences
  paste(seqs, collapse="\n")
}



# ============================
# Main loop
# ============================
for (file in subject_files) {
  cat("Processing", file, "...\n")
  accs <- readLines(file)
  accs <- accs[accs != ""]
  accs <- clean_accession(accs)
  
  if (length(accs) == 0) next
  
  fasta_out <- c()
  batch_size <- 100
  for (i in seq(1, length(accs), by=batch_size)) {
    batch <- accs[i:min(i + batch_size - 1, length(accs))]
    seqs <- entrez_fetch(db="protein", id=batch, rettype="fasta", retmode="text")
    
    # Append organism names for SwissProt entries in this batch
    seqs <- add_organism_to_swissprot(seqs, batch)
    
    fasta_out <- c(fasta_out, seqs)
    Sys.sleep(0.4)
  }
  
  out_file <- file.path(sub("^subject_(.*)\\.txt$", "fasta_\\1.fasta", basename(file)))
  writeLines(fasta_out, out_file)
  cat("  -> Saved", length(accs), "sequences to", out_file, "\n")
}
