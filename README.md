# G-FORCER
Automated, iterative curation of protein multiple sequence alignments. Removes alignment-distorting sequences (gap-forcing insertions, long deletions, high gap fraction) using biologically-aware logic, protected sequences and species-level rescue. R + MAFFT.
# G-FORCER v1.2

**Gap-FORcing sequenCe itErative Remover**

**Authors:** RAG, MMM & SMA

\- Rivero, Alexander Gabriel (main developer)

https://orcid.org/0009-0001-5813-7865

\- Manifesto, María Marcela (supervisor / co-author)

https://orcid.org/0000-0001-7333-0826

\- Soria, Marcelo Abel (supervisor / co-author)

https://orcid.org/0000-0001-8556-147X

## Overview

**G-FORCER** is an R-based command-line tool designed to iteratively
remove problematic protein sequences from multiple sequence alignments
(MSAs).

It targets sequences that:

- introduce large, alignment-distorting insertions (*gap-forcing
  sequences*),
- contain long deletion blocks in otherwise conserved regions,
- exhibit excessive overall gap content,
- or show abnormal length relative to the alignment core.

G-FORCER applies a **biologically aware, iterative strategy** while
supporting **explicit protection rules** for:

- proteins of interest (POI),
- reference sequences,
- and optional species-level protection.

⚠️ **Important:**  
G-FORCER removes **sequences**, not alignment columns.  
Column-based logic is used **only** to decide which sequences should be
removed.

## Key Concepts

G-FORCER operates in **iterative runs** over one or more MSAs:

1.  Identify *gap-forcing columns* (columns dominated by gaps).
2.  Detect sequences with **consecutive residues** in those columns
    (likely divergent insertions).

3.  Identify **core columns** (highly conserved positions).
4.  Detect **long deletion blocks** across core columns.
5.  Flag sequences with a **high overall gap fraction**.
6.  Remove flagged sequences.
7.  Re-align the remaining sequences with MAFFT.
8.  Repeat automatically until no further sequences are removed.

Protection mechanisms ensure that biologically relevant sequences are
never lost.

### The Steps Are Progressive (Hierarchical Filtering)

G-FORCER does **not** apply all removal criteria at once in a single pass. Instead, it uses a **strict hierarchical order** inside each iteration — this is one of the most important design choices.

The logic is:

1. **First check for gap-forcing sequences** (divergent insertions that create columns dominated by gaps).  
   → These are considered the **primary culprits** — they force gaps in many other sequences, artificially inflate alignment length, and cause secondary problems (long branches in trees, excessive overall gaps, apparent deletions elsewhere).

2. **Only if no sequences were flagged by gap-forcing** → check for **length outliers** (too long or too short).

3. **Only if still no removals** → check for **long deletion blocks** in conserved core columns.

4. **Only if still no removals** → check for **excessive overall gap fraction**.

**Why this order is crucial:**

- Gap-forcing sequences are the **root cause** of many alignment artifacts.  
  Removing them first often fixes secondary issues (length extremes, high gap %, deletion blocks) in the next iteration after re-alignment with MAFFT.

- If we removed everything at once, **almost all sequences could be flagged** because the presence of one or two divergent/gap-forcing sequences makes others look bad (too many gaps, abnormal relative length, etc.).

- By being progressive and re-aligning after each removal round, G-FORCER gives the alignment a chance to “heal” — new conserved regions may appear, gap-forcing columns may disappear, and fewer sequences end up being removed overall.

This mimics (but automates) the ideal manual workflow: remove the worst offenders first → re-align → re-evaluate → repeat until stable.

**Result**: Cleaner, more reliable MSAs with much less risk of over-filtering biologically meaningful sequences.

## Input Requirements

For each dataset, G-FORCER expects **paired files with the same
basename**:

example.msa

example.fasta

- .msa → MAFFT-aligned protein MSA
- .fasta → unaligned protein sequences (used for re-alignment)

### FASTA Headers and Species Tags (optional)

FASTA headers may optionally include species tags:

\>Seq123 some description \[Arabidopsis thaliana\]

Species tags are only used if **species-level protection** is enabled
(-spp).

More info in “Adding species annotations to protein sequences from BLAST
output” section.

## Installation

### **Required Software**

- **R ≥ 4.0**
- **MAFFT**
- **R package:** Biostrings

### Install Biostrings

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("Biostrings")

## One-Time Configuration (MAFFT)

G-FORCER requires the path to the MAFFT executable.  
This is defined **once** using a configuration file.

### Step 1: Create gforcer.conf

From the directory containing G-FORCER.R:

**Linux**

echo "MAFFT_PATH=$(which mafft.bat)" \> gforcer.conf

**Windows (using MAFFT for Windows)**

MAFFT_PATH=C:\path_to\_\mafft.bat

mafft.bat must be accessible from the working directory or system PATH.

### Step 2: Done

- G-FORCER automatically reads gforcer.conf.

- If the file is missing or invalid, execution stops with a clear error
  message.

## Usage

Rscript G-FORCER.R -msa_file \<path_or_file\> \[options\]

-msa_file can be:

- a **single** .msa **file**, or

- a **directory** containing multiple .msa files

## Windows Usage (Important)

On Windows, **run G-FORCER from** cmd.exe, not PowerShell.

Example invocation (tested):

"C:\Path\To\R\x64\Rscript" G-FORCER.R -msa_file 1.1.1.198.msa -PROTECTED_POI_idx 1-3 -PROTECTED_REF_idx 4-10 -spp

## Parameters

### Core Parameters

| Option | Description | Default |
|:--------------------|:------------------------------------------|:-------|
| -min_con | Minimum number of consecutive residues in gap-forcing columns required to flag a sequence | 3 |
| -th | Fraction of sequences with gaps required for a column to be considered gap-forcing | 0.95 |
| -th2 | Columns where ≤ this fraction of sequences contain gaps are considered core columns | 0.05 |
| -cf | Fraction of ungapped sequence length defining a long deletion block | 0.2 |
| -max_len_factor | Maximum allowed sequence length relative to reference mean (if defined) or median length | 1.3 |
| -min_len_factor | Minimum allowed sequence length relative to reference mean (if defined) or median length | 0.7 |
| -gap_frac | Overall gap fraction threshold. Sequences with more than this fraction of gaps across the entire alignment are removed | 0.3 |

### Sequence Protection Options (all optional)

| Option | Description |
|:-----------------------------|:----------------------------------------|
| -PROTECTED_POI_idx x-y,a | Protect Protein(s) of Interest by MSA index |
| -PROTECTED_REF_idx x-y,a | Protect reference sequences by MSA index |
| -spp | Enable species-level protection using \[Species\] tags |

**Default behavior:**  
No sequences are protected unless explicitly specified.

## Protection Logic Summary

- **PROTECTED_POI**  
  User-defined indices, never removed.

- **PROTECTED_REF**  
  User-defined reference sequences, never removed.

- **Species protection (**-spp**)**

  - Disabled by default

  - Activated only if *all* sequences of a species would be removed

  - Guarantees at least **one representative per species**

  - Selects the least problematic sequence based on internal metrics

Protected sequences are **included** when counting species.

## Output Files

All outputs are written next to the input MSA(s):

- seq_to_remove\_\<msa_name\>.tsv  
  Per-run table describing sequence status and removal logic

- \<msa_name\>.runX.fasta  
  FASTA file after sequence removal (per iteration)

- \<msa_name\>.runX.msa  
  Re-aligned MSA (MAFFT)

- GFORCER_run.log  
  Full execution log, including parameters and MAFFT commands

## Sequence Status in Output Tables (.tsv)

Each row in seq_to_remove\_\<msa\>.tsv represents **one sequence** and
is described by **four independent concepts**:

1.  **Final action** – what happened to the sequence

2.  **Primary removal reason** – why it was removed (if applicable)

3.  **Species rescue flag** – whether it was force-kept to preserve a
    species

4.  **Protection status** – whether it was explicitly protected by the
    user

These must be interpreted **together**.

### 1. Final Action (final_action column)

This column reports the **final outcome** for each sequence.

| Value  | Meaning                                |
|:-------|:---------------------------------------|
| REMOVE | Sequence was removed in this iteration |
| KEEP   | Sequence was retained                  |

Notes:

- KEEP includes **clean sequences**, **protected sequences**, and
  **species-rescued sequences**

- Iteration stops automatically when no sequences have final_action =
  REMOVE

### 2. Primary Removal Reason (primary_reason column)

This column is **only meaningful for sequences marked as** REMOVE.

For sequences that are kept, this value is always NA.

| Value | Meaning |
|:--------------|:--------------------------------------------------------|
| gap_forcing | Sequence contains consecutive residues in gap-forcing columns (divergent insertions) |
| length_check | Sequence length is outside the allowed range relative to the reference or alignment core (too long or too short) |
| core_deletion | Sequence contains long deletion blocks in conserved (core) columns |
| gap_fraction | Sequence exceeds the allowed overall gap fraction |
| NA | Sequence was **not removed** |

Important:

- primary_reason = NA **does NOT mean “no problem”**  
  It simply means *the sequence was not removed*

- Protected and species-rescued sequences always have primary_reason =
  NA

### Clarification for length_check

- The allowed length range is defined by:

  - the **mean ungapped length of protected reference sequences**, if
    PROTECTED_REF_idx is provided, or

  - the **median ungapped length of all sequences**, otherwise

- A sequence is removed with length_check if its ungapped length is:

  - **greater than** max_len_factor × reference mean length, or

  - **less than** min_len_factor × reference mean length

### 3. Species Rescue Flag (rescued_by_species column)

This column indicates whether a sequence was **force-kept to preserve
species representation**.

| Value | Meaning                                                    |
|:------|:-----------------------------------------------------------|
| TRUE  | Sequence was rescued to ensure ≥1 sequence for its species |
| FALSE | Normal decision logic applied                              |

Notes:

- This column is only relevant when -spp is enabled

### 4. Protection Status (is_protected and protected_label columns)

These columns report **explicit user protection**, independent of
filtering logic.

#### is_protected

| Value | Meaning                           |
|:------|:----------------------------------|
| TRUE  | Sequence was explicitly protected |
| FALSE | Not protected                     |

#### protected_label

| Value | Meaning                                  |
|:------|:-----------------------------------------|
| POI   | Protein of Interest (-PROTECTED_POI_idx) |
| REF   | Reference sequence (-PROTECTED_REF_idx)  |
| NA    | Not protected                            |

Notes:

- Protected sequences are **never removed**

- Protected sequences always have:

  - final_action = KEEP

  - primary_reason = NA

- Protection does **not** mean the sequence is gap-free — only that it
  cannot be removed

## Examples

Process all MSAs in a folder:

Rscript G-FORCER.R -msa_file 3_MSA/

Protect first three sequences:

Rscript G-FORCER.R -msa_file example.msa -PROTECTED_POI_idx 1-3

Enable species-level protection:

Rscript G-FORCER.R -msa_file example.msa -spp

## Notes

- Species are **always counted at the start**, regardless of -spp.

- Species protection affects **rescue logic**, not initial filtering.

- Iteration stops automatically when no more sequences are removed.

- Designed for **large protein MSAs** and iterative phylogenetic
  curation.

## **Adding species annotations to protein sequences from BLAST output**

To retrieve protein FASTA sequences with species/organism names appended
(in the format \[Organism\] at the end of each header) from a BLAST
tabular output:

1.  Extract the subject accession column from your BLAST result (outfmt
    6 or similar) into a single-column text file, one accession per
    line. Supported formats include GenBank (gb\|), RefSeq (ref\|), and
    UniProt/Swiss-Prot (sp\|) accessions. See the files in
    Example_files/ for reference examples.

2.  Run the custom R script subject_fasta.R located in the repository
    root.

    - **Recommended**: Execute it inside **RStudio** to easily monitor
      progress, view console messages, and debug if needed.

    - The script uses the rentrez package to:

      - Fetch sequences in batches from NCBI's protein database.

      - Automatically append the organism name to Swiss-Prot entries
        (which do not include it by default in FASTA headers).

      - Save the annotated FASTA files.

**Output** For each input file named subject\_\*.txt, the script
generates a corresponding fasta\_\*.fasta file with organism annotations
added where applicable.

**Notes**

- Ensure you have an active internet connection (NCBI Entrez API access
  required).

- Rate limiting is handled with short delays between requests to respect
  NCBI usage guidelines.

- If you encounter API errors or timeouts with very large lists, reduce
  the batch size (currently set to 100) in the script.

**Citation**

If you use G-FORCER in published work, please acknowledge the authors
(RAG, MMM & SMA) and reference the repository:

https://github.com/Alexander-Rivero/G-FORCER

## License

This project is licensed under the MIT License- see the
\[LICENSE\](LICENSE) file for details.

MIT License

Copyright (c) 2026 Rivero, Alexander Gabriel; Manifesto, María Marcela &
Soria, Marcelo Abel

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
