# extract-orthocaller-results

## Overview

`extract-orthocaller-results` is a collection of Python and SLURM workflows for extracting, processing, and preparing orthogroup-derived protein and CDS alignments from Orthocaller/GeneRax outputs for downstream comparative genomics and molecular evolution analyses.

The pipeline was developed for large-scale teleost comparative genomics datasets and is designed to support codon-aware alignments and downstream selection analyses using tools such as HyPhy (e.g. BUSTED-E and RELAX).

The workflow performs:

* Extraction of orthogroup protein sequences from Orthocaller outputs
* Recovery of CDS sequences corresponding to extracted proteins
* Translation and frame validation
* Automated correction and filtering of frame-shifted sequences
* Generation of codon-aware alignments
* Cleaning and formatting of alignments for downstream selection analyses

Although originally developed for large-scale cavefish comparative genomics datasets, the workflow is broadly applicable to orthology-based comparative genomics projects.

---

# Pipeline Overview

## Orthocaller results ? HyPhy-ready codon alignments

### 1. Generate master summary files

Generate summary tables from Orthocaller output directories.

Scripts:

* `GenerateMasterSummaryFile2.py`
* `GenerateMasterSummaryFile2.sh`

---

### 2. Generate GeneRax key files

Parse GeneRax logs and generate lookup/key files used throughout downstream processing.

Scripts:

* `generax_log_to_key.sh`

---

### 3. Extract orthogroup protein sequences

Extract protein FASTA files corresponding to Orthocaller orthogroups.

Scripts:

* `extract_orthocaller_proteins.py`
* `extract_orthocaller_proteins.sh`

Outputs:

* Orthogroup-level protein FASTA files

---

### 4. Recover CDS sequences

Retrieve CDS sequences corresponding to extracted proteins.

Scripts:

* `GetCDSFromProteinsV5_1.py`
* `GetCDSFromProteinsV5_array.sh`

Additional rerun workflow for unmatched sequences:

* `GetCDSFromProteinsV5_1_symbol_safe.py`
* `GetCDSFromProteinsV5_1_symbol_safe_ReRun.sh`

Outputs:

* CDS FASTA files for each orthogroup

---

### 5. Translate CDS sequences

Translate CDS FASTA files into protein sequences for frame validation.

Scripts:

* `translate_CDS_folder.py`
* `translate_CDS_folder.sh`

---

### 6. Recover original protein sequences

Retrieve original source protein sequences for frame comparison and validation.

Scripts:

* `GetOriginalProts_V1.py`
* `GetOriginalProts_V1.sh`

---

### 7. Detect frame errors

Identify frame shifts and translation inconsistencies between translated CDS sequences and original proteins.

Scripts:

* `CheckFrameErrors3.py`
* `CheckFrameErrors3.sh`

Outputs:

* CSV summaries of frame failures
* Lists of problematic sequences

---

### 8. Correct translation and CDS frame errors

Attempt automated correction of frame-shifted sequences.

Protein correction:

* `FIXTRANSLATIONFRAMES.py`
* `FIXTRANSLATIONFRAMES.sh`

CDS correction:

* `correctCDSFrame.py`
* `correctCDSFrame.sh`

---

### 9. Filter remaining problematic sequences

Remove sequences that continue to fail translation/frame validation after correction.

Scripts:

* `sift_failed_sequences.py`
* `sift_failed_sequences.sh`

---

### 10. Generate codon-aware alignments

Protein alignments are generated using MAFFT and converted into codon alignments using PAL2NAL.

Scripts:

* `MAFFT_arrayScript.sh`

Outputs:

* Codon-aware nucleotide alignments

---

### 11. Clean alignments for selection analyses

Remove problematic headers, stop codons, and remaining formatting issues prior to HyPhy analyses.

Scripts:

* `cln.py`
* `cln.sh`

Outputs:

* HyPhy-ready codon alignments

---

### 12. Alignment statistics and species filtering

Generate summary statistics and species counts for downstream filtering.

Scripts:

* `count_species_pal2nal.py`
* `count_species_pal2nal.sh`
* `Alignment_Length_Stats.py`
* `Alignment_Length_Stats.sh`

Typical downstream filtering thresholds:

* =3 foreground cavefish species
* =30 background species

---

# Repository Structure

```text
extract-orthocaller-results/
+-- config/
+-- docs/
+-- helpers/
+-- scripts/
¦   +-- python/
¦   +-- slurm/
+-- README.md
```

---

# Requirements

## Software

Typical dependencies include:

* Python 3
* MAFFT
* PAL2NAL
* HyPhy
* GeneRax
* SLURM
* Conda

## Python packages

Common packages used include:

* argparse
* pathlib
* csv
* pandas
* Biopython

---

# Configuration

This repository uses configurable path variables rather than hardcoded HPC-specific paths whenever possible.

Example configuration files are located in:

```text
config/config.example.yaml
```

Users should create their own local configuration file:

```bash
cp config/config.example.yaml config/config.yaml
```

---

# Notes

This repository reflects an actively developed comparative genomics workflow originally optimized for large-scale teleost datasets generated from Orthocaller and GeneRax pipelines.

Several scripts were developed iteratively during large HPC production runs and continue to be refactored for portability and modularity.

The workflow is currently optimized for SLURM-based HPC environments.
