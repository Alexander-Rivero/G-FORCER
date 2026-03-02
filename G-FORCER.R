# Copyright (c) 2026 Rivero, Alexander Gabriel; Manifesto, María Marcela & Soria, Marcelo Abel
# Licensed under the MIT License. See LICENSE file for details.

#MIT License
#Copyright (c) 2026 Rivero, Alexander Gabriel; Manifesto, María Marcela & Soria, Marcelo Abel
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

# =============================
# G-FORCER (per-MSA outputs, automated multi-run)
# Gap-FORcing sequenCe itErative Remover
# =============================
# === Capture full command used to run G-FORCER ===
full_cmd <- paste(commandArgs(), collapse = " ")

cat("\n========================================\n")
cat("        G-FORCER v1.2\n")
cat("  Gap-FORcing sequenCe itErative Remover\n")
cat("        by RAG, MMM & SMA\n")
cat("========================================\n\n")

# ================================
# Load G-FORCER configuration
# ================================

config_file <- "gforcer.conf"

if (!file.exists(config_file)) {
  stop(
    "ERROR: G-FORCER configuration file not found.\n\n",
    "Please create a file named 'gforcer.conf' in the same directory\n",
    "with the following content:\n\n",
    "  MAFFT_PATH=/full/path/to/mafft\n\n",
    "Example:\n",
    "  echo \"MAFFT_PATH=$(which mafft)\" > gforcer.conf\n"
  )
}

config_lines <- readLines(config_file)

mafft_path <- sub("^MAFFT_PATH=", "", config_lines[grepl("^MAFFT_PATH=", config_lines)])

if (length(mafft_path) == 0 || mafft_path == "" || !file.exists(mafft_path)) {
  stop(
    "ERROR: MAFFT_PATH is not defined or invalid in gforcer.conf\n\n",
    "Please edit gforcer.conf and set:\n",
    "  MAFFT_PATH=/full/path/to/mafft\n"
  )
}

mafft_exec <- shQuote(mafft_path)

# =============================
# Command-line argument parsing
# =============================

args <- commandArgs(trailingOnly = TRUE)

# Default values for G-FORCER parameters
# -------------------------------------
# These correspond to the original hard-coded values used before CLI options.
# Users may override them from the command line, but if no CLI argument is
# provided, the pipeline will fall back to these defaults.

cli <- list(
  msa_file = NULL,   # Path to a single .msa file OR to a folder full of .msa files
  
  min_consecutive = 3,           
  # Minimum number of consecutive “gap-forcing residues” required to flag 
  # a sequence as a Gap Forcing Sequence. 
  
  gap_forcing_threshold = 0.95,              
  # Threshold for gap-forcing residues identification. If >=95% of sequences have a 
  # gap in a column, that column is considered to have a gap-forcing residue.
  
  core_deletion_threshold = 0.05,
  # Threshold for core_deletion_cols identification.
  # Columns where <=5% of sequences contain a gap.
  # Used to detect suspicious continuous regions with missing residues.
  
  consec_frac = 0.2,
  # Fraction of the total ungapped sequence length required to define
  # long consecutive blocks of core_deletion_cols gaps. Seqs with
  # block >= consec_frac are seqs with deletions.
  # Example: if a sequence is 400 aa (ungapped), 0.2 → 80 aa.
  
  PROTECTED_POI_idx = NULL,
  # Indices of the Protein(s) of Interest (POI), i.e., "PROTECTED_POI".
  # Allows explicit manual protection instead of relying on the assumption
  # that POIs appear first in the MSA.
  # User provides a range like: -PROTECTED_POI_idx "1-4"  (protect sequences 1..4)
  # If NULL → no sequence protected
  
  PROTECTED_REF_idx = NULL,
  # Indices of reference sequences to protect (e.g., SwissProt).
  # User provides a range: -PROTECTED_REF_idx "78-95"
  # If NULL → no sequence protected
  
  spp = FALSE,
  # Enable species-level protection based on FASTA headers [Species]
  # Default: FALSE (no species protection).
  # Use -spp to activate.
  
  max_len_factor = 1.3,
  # Maximum allowed length relative to reference/median length.
  # Sequence length > mean_ref * max_len_factor → TOO_LONG
  
  min_len_factor = 0.7,
  # Minimum allowed length relative to reference/median length.
  # Sequence length < mean_ref * min_len_factor → TOO_SHORT
  
  gap_frac_threshold = 0.3   
  # NEW: overall gap fraction threshold for removal
)

# =============================
# Help (-h / -help) and Version (-v / -version) screens
# =============================
if (length(args) > 0) {
  if (any(args %in% c("-h", "-help", "--help"))) {
    cat("\n========================================\n")
    cat(" G-FORCER v1.2 — HELP\n")
    cat("========================================\n\n")
  
  cat("Usage:\n")
  cat("  Rscript G-FORCER.R -msa_file <path> [options]\n\n")
  
  cat("Flags:\n")
  cat(" -h, -help, --help     Show this help message\n")
  cat(" -v, -version, --version  Show version information\n\n")
  
  cat("Required argument:\n")
  cat("  -msa_file <path>\n")
  cat("        Path to ONE .msa file OR a directory containing multiple .msa files\n\n")
  
  cat("Optional parameters (default values shown):\n")
  cat(sprintf(
    "  -max_len_factor <float>\n        Maximum allowed length relative to reference or median length.\n        Sequences longer than mean_ref × factor are flagged TOO_LONG (default: %.2f)\n\n",
    cli$max_len_factor
  ))
  cat(sprintf(
    "  -min_len_factor <float>\n        Minimum allowed length relative to reference or median length.\n        Sequences shorter than mean_ref × factor are flagged TOO_SHORT (default: %.2f)\n\n",
    cli$min_len_factor
  ))
  cat(sprintf(
    "  -min_con <int>\n        Minimum number of consecutive gap-forcing residues\n        required to flag a sequence (default: %d)\n\n",
    cli$min_consecutive
  ))
  
  cat(sprintf(
    "  -th <float>\n        Gap-forcing residue threshold.\n        A column is considered gap-forcing if >= this fraction of sequences\n        contain a gap (default: %.2f)\n\n",
    cli$gap_forcing_threshold
  ))
  
  cat(sprintf(
    "  -th2 <float>\n        Core-deletion threshold.\n        Columns where <= this fraction of sequences contain a gap\n        are used to detect deletion blocks (default: %.2f)\n\n",
    cli$core_deletion_threshold
  ))
  
  cat(sprintf(
    "  -cf <float>\n        Fraction of ungapped sequence length defining a long deletion block.\n        Sequences with deletion blocks >= this fraction are removed\n        (default: %.2f)\n\n",
    cli$consec_frac
  ))
  
  cat(sprintf(
    " -gap_frac <float>\n Overall gap fraction threshold.\n Sequences with > this fraction of gaps across the entire alignment\n are removed in the final step (default: %.2f)\n\n",
    cli$gap_frac_threshold
  ))
  
  cat("Sequence protection options (ALL disabled by default):\n")
  cat("  -PROTECTED_POI_idx x-y,a,b\n")
  cat("        Protect Protein(s) of Interest by MSA index (e.g. 1-3)\n\n")
  
  cat("  -PROTECTED_REF_idx x-y,a\n")
  cat("        Protect reference sequences by MSA index (e.g. 4-10)\n\n")
  
  cat("  -spp\n")
  cat("        Enable species-level protection using [Species] tags in FASTA headers\n\n")
  
  cat("Examples:\n")
  cat("  Rscript G-FORCER.R -msa_file 3_MSA\n")
  cat("  Rscript G-FORCER.R -msa_file data.msa -max_len_factor 1.5 -min_len_factor 0.6\n")
  cat("  Rscript G-FORCER.R -msa_file input.msa -min_con 4\n")
  cat("  Rscript G-FORCER.R -msa_file x.msa -PROTECTED_POI_idx 1-3\n")
  cat("  Rscript G-FORCER.R -msa_file x.msa -spp\n\n")
  
  cat("Notes:\n")
  cat("  • All sequence protection mechanisms are OFF by default.\n")
  cat("  • If both PROTECTED_POI_idx and PROTECTED_REF_idx are provided, they must NOT overlap.\n")
  cat("  • FASTA files must share the same basename as the .msa file.\n")
  cat("  • G-FORCER removes SEQUENCES, not columns. Column logic is internal.\n\n")
  
  cat("========================================\n\n")
  quit(save = "no")
}
  if (any(args %in% c("-v", "-version", "--version"))) {
    cat("\n========================================\n")
    cat(" G-FORCER v1.2\n")
    cat(" Gap-FORcing sequenCe itErative Remover\n")
    cat(" Authors: RAG, MMM & SMA\n")
    cat(" Release date: December 2025\n")  # or whatever you prefer
    cat("========================================\n\n")
    quit(save = "no")
  }
}

# Parse args: -key value
if (length(args) > 0) {
  # First pass: detect flags
  if ("-spp" %in% args) cli$spp <- TRUE
  # Second pass: parse key-value args
  for (i in seq(1, length(args), by = 2)) {
    key <- args[i]
    # skip flags
    if (key %in% c("-spp")) next
    if (i + 1 <= length(args)) {
      val <- args[i + 1]
    } else {
      stop(paste("Missing value for argument", key))
    }
    if (key == "-msa_file") cli$msa_file <- val
    if (key == "-min_con")  cli$min_consecutive <- as.numeric(val)
    if (key == "-th")       cli$gap_forcing_threshold <- as.numeric(val)
    if (key == "-th2")      cli$core_deletion_threshold <- as.numeric(val)
    if (key == "-cf")       cli$consec_frac <- as.numeric(val)
    if (key == "-gap_frac") cli$gap_frac_threshold <- as.numeric(val)
    if (key == "-PROTECTED_POI_idx")   cli$PROTECTED_POI_idx <- val
    if (key == "-PROTECTED_REF_idx")  cli$PROTECTED_REF_idx <- val
  }
}
cli$PROTECTED_POI_idx_raw  <- cli$PROTECTED_POI_idx
cli$PROTECTED_REF_idx_raw <- cli$PROTECTED_REF_idx


# =============================
# Validate numeric parameters
# =============================

if (is.na(cli$min_consecutive) || cli$min_consecutive < 1) {
  stop("ERROR: -min_con must be an integer >= 1")
}

if (is.na(cli$gap_forcing_threshold) || cli$gap_forcing_threshold < 0 || cli$gap_forcing_threshold > 1) {
  stop("ERROR: -th must be a float between 0 and 1")
}

if (is.na(cli$core_deletion_threshold) || cli$core_deletion_threshold < 0 || cli$core_deletion_threshold > 1) {
  stop("ERROR: -th2 must be a float between 0 and 1")
}

if (is.na(cli$consec_frac) || cli$consec_frac < 0 || cli$consec_frac > 1) {
  stop("ERROR: -cf must be a float between 0 and 1")
}
if (is.na(cli$max_len_factor) || cli$max_len_factor <= 1) {
  stop("ERROR: -max_len_factor must be > 1")
}

if (is.na(cli$min_len_factor) || cli$min_len_factor <= 0 || cli$min_len_factor >= 1) {
  stop("ERROR: -min_len_factor must be between 0 and 1")
}

if (is.na(cli$gap_frac_threshold) || cli$gap_frac_threshold < 0 || cli$gap_frac_threshold > 1) {
  stop("ERROR: -gap_frac must be a float between 0 and 1")
}

# Helper: parse index ranges "1-5,7,10-12"
parse_idx <- function(txt) {
  parts <- unlist(strsplit(txt, ","))
  out <- c()
  for (p in parts) {
    if (grepl("-", p)) {
      rng <- as.numeric(unlist(strsplit(p, "-")))
      out <- c(out, seq(rng[1], rng[2]))
    } else {
      out <- c(out, as.numeric(p))
    }
  }
  unique(out)
}
if (!is.null(cli$PROTECTED_POI_idx)) {
  cli$PROTECTED_POI_idx <- parse_idx(cli$PROTECTED_POI_idx)
}
if (!is.null(cli$PROTECTED_REF_idx)) {
  cli$PROTECTED_REF_idx <- parse_idx(cli$PROTECTED_REF_idx)
}
if (!is.null(cli$PROTECTED_POI_idx) && !is.null(cli$PROTECTED_REF_idx)) {
  if (length(intersect(cli$PROTECTED_POI_idx, cli$PROTECTED_REF_idx)) > 0) {
    stop("ERROR: PROTECTED_POI_idx and PROTECTED_REF_idx overlap — FIX ARGUMENTS.")
  }
}

# assign parameter variables used later in the script
min_consecutive <- cli$min_consecutive
gap_forcing_threshold <- cli$gap_forcing_threshold
core_deletion_threshold <- cli$core_deletion_threshold
consec_frac <- cli$consec_frac
max_len_factor <- cli$max_len_factor
min_len_factor <- cli$min_len_factor
gap_frac_threshold <- cli$gap_frac_threshold

# Load library for reading FASTA/MSA
library(Biostrings)

# ---------- Build 'msa_files' list ----------
# If the user provides -msa_file it may be a directory or a single file.
# If not provided, we STOP because we have no default 3_MSA anymore.

if (is.null(cli$msa_file)) {
  stop("ERROR: You must provide -msa_file pointing to a folder or one .msa file\n")
}

# expand ~ and make absolute
user_path <- tryCatch(normalizePath(cli$msa_file, mustWork = FALSE),
                      error = function(e) cli$msa_file)

if (dir.exists(user_path)) {
  # list ONLY .msa (per your new rule)
  msa_files <- list.files(user_path,
                          pattern = "\\.msa$",
                          full.names = TRUE,
                          ignore.case = TRUE)
  if (length(msa_files) == 0) {
    stop("No .msa files found in directory: ", cli$msa_file)
  }
} else if (file.exists(user_path)) {
  # single file path provided
  if (!grepl("\\.msa$", user_path, ignore.case = TRUE)) {
    stop("Provided file is not a .msa file: ", cli$msa_file)
  }
  msa_files <- user_path
} else {
  stop("Provided -msa_file does not exist: ", cli$msa_file)
}
# Define msa_folder (where MSA+FASTA files live, and where all outputs must go)
if (dir.exists(user_path)) {
  msa_folder <- user_path
} else {
  msa_folder <- dirname(user_path)
}
## ---- OPEN LOGGING HERE ----
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_filename <- sprintf("GFORCER_run_%s.log", timestamp)
log_con <- file(file.path(msa_folder, log_filename), open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
cat("\n========================================\n")
cat("        G-FORCER v1.2\n")
cat("  Gap-FORcing sequenCe itErative Remover\n")
cat("        by RAG, MMM & SMA\n")
cat("========================================\n\n")
cat("=== Command used ===\n")
cat(full_cmd, "\n")
cat("====================\n\n")
# Print parameter summary
cat("=== G-FORCER PARAMETERS ===\n")

params <- list(
  msa_file               = cli$msa_file,
  min_con                = cli$min_consecutive,
  th                     = cli$gap_forcing_threshold,
  th2                    = cli$core_deletion_threshold,
  cf                     = cli$consec_frac,
  gap_frac              = cli$gap_frac_threshold,
  PROTECTED_POI_idx      = cli$PROTECTED_POI_idx_raw,
  PROTECTED_REF_idx      = cli$PROTECTED_REF_idx_raw,
  spp                    = cli$spp,
  max_len_factor         = cli$max_len_factor,
  min_len_factor         = cli$min_len_factor
)

for (nm in names(params)) {
  val <- params[[nm]]
  if (is.null(val)) val <- "NULL"
  cat(sprintf(
    "%-25s = %s\n",
    nm,
    paste(val, collapse = ", ")
  ))
}

cat("===========================\n")
# ---- After building msa_files, verify associated FASTA ----
cat("=== Checking associated FASTA files ===\n")

for (m in msa_files) {
  base <- tools::file_path_sans_ext(m)
  fasta_file <- paste0(base, ".fasta")
  
  cat("MSA:", basename(m), " → ")
  
  if (file.exists(fasta_file)) {
    cat("associated FASTA found:", basename(fasta_file), "\n")
  } else {
    stop("ERROR: No associated FASTA found for ", basename(m),
         "\nExpected: ", fasta_file, "\n")
  }
}

cat("All associated FASTA files verified.\n")
cat("=====================================\n\n")

# ---------- done building msa_files ----------

# --- Updated function with species protection ---
identify_problematic_seqs <- function(seq_mat,
                                      gap_forcing_threshold,
                                      min_consecutive,
                                      core_deletion_threshold,
                                      consec_frac,
                                      gap_frac_threshold,
                                      PROTECTED_POI_idx,
                                      PROTECTED_REF_idx,
                                      max_len_factor,
                                      min_len_factor,
                                      enable_species_protection = FALSE)
{
  
  n_seq <- nrow(seq_mat)
  n_col <- ncol(seq_mat)
  
  seq_names <- rownames(seq_mat)
  if (is.null(seq_names)) seq_names <- as.character(1:n_seq)
  
  is_CG  <- seq_len(n_seq) %in% PROTECTED_POI_idx
  is_ref <- seq_len(n_seq) %in% PROTECTED_REF_idx
  is_protected <- is_CG | is_ref
  
  seq_len_no_gap <- rowSums(seq_mat != "-")
  
  if (any(is_ref)) {
    PROTECTED_REF_len_no_gap <- seq_len_no_gap[is_ref]
    mean_ref_length <- ceiling(mean(PROTECTED_REF_len_no_gap))
    max_len <- ceiling(mean_ref_length * max_len_factor)
    min_len <- floor(mean_ref_length * min_len_factor)
      } else {
    median_len <- ceiling(median(seq_len_no_gap))
    max_len <- ceiling(median_len * max_len_factor)
    min_len <- floor(median_len * min_len_factor)
    PROTECTED_REF_len_no_gap <- numeric(0)
    mean_ref_length <- NA
  }
  
  threshold_count <- floor(gap_forcing_threshold * n_seq)
  if (threshold_count <= 1) threshold_count <- 1
  gap_counts <- colSums(seq_mat == "-")
  gap_forcing_cols <- which(gap_counts >= threshold_count)
  
  if (any(is_protected) && length(gap_forcing_cols) > 0) {
    prot_has_res <- colSums(seq_mat[is_protected, gap_forcing_cols, drop = FALSE] != "-") > 0
    gap_forcing_cols <- gap_forcing_cols[!prot_has_res]
  }
  
  gap_forcing_count <- integer(n_seq)
  final_action       <- rep("KEEP", n_seq)
  primary_reason     <- rep(NA_character_, n_seq)
  rescued_by_species <- rep(FALSE, n_seq)
  
    # --- Species extraction (updated: include protected sequences) ---
  seq_species <- rep(NA_character_, n_seq)
  for (i in seq_len(n_seq)) {
    sp_match <- regmatches(seq_names[i], regexpr("\\[.*\\]$", seq_names[i]))
    if (length(sp_match) > 0) seq_species[i] <- gsub("\\[|\\]", "", sp_match)
  }
  
  for (i in 1:n_seq) {
    gap_forcing_count[i] <- sum(seq_mat[i, gap_forcing_cols] != "-")
    if (is_protected[i]) next
    
    # Step 2: consecutive rare columns
    seq_flags <- as.integer(seq_mat[i, gap_forcing_cols] != "-")
    if (length(seq_flags) == 0 || all(seq_flags == 0)) {
      gap_forcing_run_flag <- FALSE
    } else {
      rle_flags <- rle(seq_flags)
      gap_forcing_run_flag <- any(rle_flags$values == 1 & rle_flags$lengths >= min_consecutive)
    }
    
    if (gap_forcing_run_flag) {
      final_action[i]   <- "REMOVE"
      primary_reason[i] <- "gap_forcing"
    } else if (seq_len_no_gap[i] > max_len) {
      final_action[i]   <- "REMOVE"
      primary_reason[i] <- "length_check"
    } else if (seq_len_no_gap[i] < min_len) {
      final_action[i]   <- "REMOVE"
      primary_reason[i] <- "length_check"
    } else {
      final_action[i] <- "KEEP"
      primary_reason[i] <- NA_character_
    }
  }
  total_problem_gaps <- gap_forcing_count
  # Step 3: problem_col_2
  non_prot_indices <- which(!is_protected)
  deletion_block_run_executed <- FALSE
  if (length(non_prot_indices) > 0) deletion_block_run_executed <- all(final_action[non_prot_indices] == "KEEP")

  deletion_gap_count <- integer(n_seq)
  
  if (deletion_block_run_executed) {
    col_gap_count <- colSums(seq_mat == "-")
    min_count <- floor(core_deletion_threshold * n_seq)
    if (min_count < 1) min_count <- 1
    core_deletion_cols <- which(col_gap_count <= min_count & col_gap_count > 0)
    
    if (any(is_protected) && length(core_deletion_cols) > 0) {
      prot_has_res2 <- colSums(seq_mat[is_protected, core_deletion_cols, drop = FALSE] != "-") > 0
      core_deletion_cols <- core_deletion_cols[!prot_has_res2]
    }
    
    for (i in 1:n_seq) {
      if (is_protected[i]) {
        if (length(core_deletion_cols) > 0) deletion_gap_count[i] <- sum(as.integer(seq_mat[i, core_deletion_cols] == "-"))
        next
      }
      
      flags2_full <- integer(n_col)
      if (length(core_deletion_cols) > 0) flags2_full[core_deletion_cols] <- as.integer(seq_mat[i, core_deletion_cols] == "-")
      
      rle2 <- rle(flags2_full)
      X <- ceiling(consec_frac * seq_len_no_gap[i])
      if (X < 1) X <- 1
      
      if (any(rle2$values == 1 & rle2$lengths >= X)) {
        final_action[i] <- "REMOVE"
        primary_reason[i] <- "deletion_block"
      }
      
      if (length(core_deletion_cols) > 0) deletion_gap_count[i] <- sum(flags2_full[core_deletion_cols])
    }
    
    total_problem_gaps <- gap_forcing_count + deletion_gap_count
  }
  
  # Step 4: high gap fraction
  step4_executed <- FALSE
  if (length(non_prot_indices) > 0)  step4_executed <- all(final_action[non_prot_indices] == "KEEP")

  if (step4_executed) {
    for (i in non_prot_indices) {
      if ((sum(seq_mat[i, ] == "-") / n_col) > gap_frac_threshold) {
        final_action[i] <- "REMOVE"
        primary_reason[i] <- "gap_fraction"
      }
    }
  }
  
  # --- Species protection logic (updated: consider protected sequences and skip if any unflagged exist) ---
  if (enable_species_protection &&
      any(!is_protected & !is.na(seq_species))) {    
    # Identify species already protected (CG or sp)
    protected_species <- seq_species[is_protected]
    
    species_list <- unique(seq_species[!is.na(seq_species)])
    for (sp in species_list) {
      # Skip species already protected
      if (sp %in% protected_species) next
      
      # indices of non-protected sequences for this species
      idx <- which(seq_species == sp & !is_protected)
      if (length(idx) == 0) next
      
      # If any non-protected sequence is already KEEP, do NOT force-keep any
      if (any(final_action[idx] == "KEEP")) next
      
      # Now all non-protected seqs for this species are flagged -> we must keep one
      if (length(idx) == 1) {
        final_action[idx]       <- "KEEP"
        primary_reason[idx]     <- NA_character_
        rescued_by_species[idx] <- TRUE
      } else {
        # robust mapping of step_source to numeric order
        step_order <- c( "gap_forcing"   = 2,"deletion_block" = 3, "gap_fraction" = 4, "NA"= 0)
        ss <- as.character(primary_reason[idx])
        ss[is.na(ss)] <- "NA"
        steps <- step_order[ss]
        max_step <- max(steps, na.rm = TRUE)
        candidates <- idx[which(steps == max_step)]
        
        if (length(candidates) == 1) {
          final_action[candidates] <- "KEEP"
          primary_reason[candidates] <- NA_character_
          rescued_by_species[candidates] <- TRUE
          } else {
          if (max_step == 2) {
            lengths <- seq_len_no_gap[candidates]
            med_len <- median(seq_len_no_gap[non_prot_indices])
            #picks the sequence whose length is closest to the median length.
            keep_idx <- candidates[which.min(abs(lengths - med_len))]
          } else if (max_step == 3) {
            rle_lengths <- sapply(candidates, function(i) {
              flags2_full <- integer(n_col)
              if (length(core_deletion_cols) > 0) flags2_full[core_deletion_cols] <- as.integer(seq_mat[i, core_deletion_cols] == "-")
              rle_vals <- rle(flags2_full)
              if (any(rle_vals$values == 1 & rle_vals$lengths >= ceiling(consec_frac * seq_len_no_gap[i]))) {
                min(rle_vals$lengths[rle_vals$values == 1 & rle_vals$lengths >= ceiling(consec_frac * seq_len_no_gap[i])])
              } else Inf
            })
            #keep the seq with the shortest deletion
            keep_idx <- candidates[which.min(rle_lengths)]
          } else if (max_step == 4) {
            gap_frac <- sapply(candidates, function(i) sum(seq_mat[i, ] == "-") / n_col)
            keep_idx <- candidates[which.min(gap_frac)]
          } else {
            keep_idx <- candidates[1]
          }
          final_action[keep_idx] <- "KEEP"
          primary_reason[keep_idx] <- NA_character_
          rescued_by_species[keep_idx] <- TRUE
        }
      }
    }
  }
  # ------------------------------------------------------
  
  # --- Build final result table ---
  result <- data.frame(
    Seq_Name           = seq_names,
    final_action       = final_action,
    primary_reason     = primary_reason,
    rescued_by_species = rescued_by_species,
    total_problem_gaps = total_problem_gaps,
    GapForcingResidues = gap_forcing_count,
    stringsAsFactors = FALSE
  )
  
  # Include deletion block gaps if executed
  if (deletion_block_run_executed) {
    result$DeletionBlockGaps <- deletion_gap_count
  }
  
  # --- Protected sequence info (like original code) ---
  result$is_protected <- is_protected
  # Numeric codes for ordering (internal)
  protected_type <- rep(NA_integer_, n_seq)
  protected_type[is_CG]  <- 1
  protected_type[is_ref] <- 2
  # Human-readable labels for the table
  result$protected_label <- rep(NA_character_, n_seq)
  result$protected_label[is_CG]  <- "POI"
  result$protected_label[is_ref] <- "REF"
  
  # --- Ordering for presentation ---
  ord_remove    <- result$final_action == "REMOVE"
  ord_protected <- result$is_protected
  ord_prot_type <- protected_type
  ord_species   <- result$rescued_by_species
  ord_gaps      <- result$total_problem_gaps
  
  result <- result[
    order(
      ord_remove,
      ord_protected,
      ord_prot_type,
      ord_species,
      -ord_gaps
    ),
  ]
  
  return(result)
}

# --- Main loop remains unchanged ---
for (msa_file in msa_files) {
  base_name <- tools::file_path_sans_ext(basename(msa_file))
  fasta_base <- file.path(msa_folder, paste0(base_name, ".fasta"))
  msa_input <- msa_file
  
  run <- 1
  repeat {
    cat("========== RUN", run, "==========\n")
    cat("Processing:", msa_input, "\n")
    
    msa <- readAAStringSet(msa_input)
    # Count unique species at start
    seq_species_all <- sapply(names(msa), function(h) {
      sp_match <- regmatches(h, regexpr("\\[.*\\]$", h))
      if(length(sp_match) > 0) gsub("\\[|\\]", "", sp_match) else NA
    })
    unique_species_start <- length(unique(seq_species_all[!is.na(seq_species_all)]))
    cat("Unique species at start:", unique_species_start, "\n")
    seq_mat <- do.call(rbind, strsplit(as.character(msa), ""))
    rownames(seq_mat) <- names(msa)
    
    seq_headers <- names(msa)
    total_seqs <- length(seq_headers)
    
    if (is.null(cli$PROTECTED_POI_idx)) {
      PROTECTED_POI_idx <- numeric(0)
    }
    if (is.null(cli$PROTECTED_REF_idx)) {
      PROTECTED_REF_idx <- numeric(0)
    }
    
    if (!is.null(cli$PROTECTED_POI_idx) || !is.null(cli$PROTECTED_REF_idx)) {
      # ===== USER-DEFINED INDEX MODE=====
      if (!is.null(cli$PROTECTED_POI_idx)) {
        PROTECTED_POI_idx <- cli$PROTECTED_POI_idx
      } else {
        PROTECTED_POI_idx <- numeric(0)
      }
      if (!is.null(cli$PROTECTED_REF_idx)) {
        PROTECTED_REF_idx <- cli$PROTECTED_REF_idx
      } else {
        PROTECTED_REF_idx <- numeric(0)
      }
      
      # Validar rangos en cada MSA
      if (any(PROTECTED_POI_idx > total_seqs)) {
        stop("ERROR: some PROTECTED_POI_idx exceed number of sequences in MSA: ", total_seqs)
      }
      if (any(PROTECTED_REF_idx > total_seqs)) {
        stop("ERROR: some PROTECTED_REF_idx exceed number of sequences in MSA: ", total_seqs)
      }
      
    } 
    
    result_table <- identify_problematic_seqs(
      seq_mat,
      gap_forcing_threshold,
      min_consecutive,
      core_deletion_threshold,
      consec_frac,
      gap_frac_threshold,
      PROTECTED_POI_idx,
      PROTECTED_REF_idx,
      max_len_factor,
      min_len_factor,
      enable_species_protection = cli$spp
    )
    
    # Count unique species among sequences that will be kept
    keep_idx <- result_table$final_action == "KEEP"
    seq_species_kept <- sapply(result_table$Seq_Name[keep_idx], function(h) {
    sp_match <- regmatches(h, regexpr("\\[.*\\]$", h))
      if (length(sp_match) > 0) gsub("\\[|\\]", "", sp_match) else NA
    })
    unique_species_end <- length(unique(seq_species_kept[!is.na(seq_species_kept)]))
        if (cli$spp) {
      cat("Unique species after protection:", unique_species_end, "\n")
    } else {
      cat("Unique species after filtering:", unique_species_end, "\n")
    }
    
      out_file <- file.path(
      msa_folder, paste0("seq_to_remove_", base_name,
             ifelse(run==1, "", paste0("_run", run-1)),".tsv"))
    write.table(result_table, out_file, sep="\t", quote=FALSE, row.names=FALSE)
    cat("Written:", out_file, "\n")
    
    final_counts <- table(result_table$final_action)
    cat("Final actions:\n")
    cat("  REMOVE :", ifelse(!is.na(final_counts["REMOVE"]), final_counts["REMOVE"], 0), "\n")
    cat("  KEEP   :", ifelse(!is.na(final_counts["KEEP"]),   final_counts["KEEP"],   0), "\n\n")
    reason_counts <- table(result_table$primary_reason)
        cat("Primary reasons detected:\n")
    for (r in names(reason_counts)) {
      cat(" ", r, ":", reason_counts[r], "\n")
    }
    cat("\n")
    prot_counts <- table(result_table$protected_label)
    cat("Protection summary:\n")
    cat("  PROTECTED_POI :", ifelse("POI" %in% names(prot_counts),
                                    prot_counts["POI"], 0), "\n")
    cat("  PROTECTED_REF :", ifelse("REF" %in% names(prot_counts),
                                    prot_counts["REF"], 0), "\n\n")
    
      has_removals <- any(result_table$final_action == "REMOVE")
    if (!has_removals) {
      cat("No sequence to REMOVE found. Stopping at run", run-1, "\n")
      break
    }
    
    fasta_in  <- if (run == 1) fasta_base else file.path(msa_folder, paste0(base_name, ".run", run - 1, ".fasta"))
    fasta_out <- paste0(base_name, ".run", run, ".fasta")
    msa_out   <- paste0(base_name, ".run", run, ".msa")
    
    # --- Read report and get sequences to remove ---
    report_df <- read.table(out_file,
                            header = TRUE,
                            sep = "\t",
                            quote = "",
                            comment.char = "",
                            stringsAsFactors = FALSE)
    
    remove_ids <- report_df$Seq_Name[
      report_df$final_action == "REMOVE"]
    
    # --- Filter FASTA using Biostrings ---
    fasta_seqs <- readAAStringSet(fasta_in)
    
    keep <- !names(fasta_seqs) %in% remove_ids
    fasta_filtered <- fasta_seqs[keep]
    
    writeXStringSet(fasta_filtered,
                    file.path(msa_folder, fasta_out))
    mafft_cmd <- paste0(mafft_exec," --auto ",shQuote(file.path(msa_folder, fasta_out)),
      " > ",shQuote(file.path(msa_folder, msa_out)))
    
    cat("Running mafft.bat...\n")
    system(mafft_cmd)
    
    msa_input <- file.path(msa_folder, msa_out)
    run <- run + 1
  }

}

## ---- CLOSE LOGGING HERE ----
sink(type = "message")
sink(type = "output")
close(log_con)


