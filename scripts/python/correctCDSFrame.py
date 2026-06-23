```python
#!/usr/bin/env python3
"""
Adjust CDS sequences so their translated proteins best match a reference
protein FASTA.

For each CDS sequence:
    1) Match the corresponding protein sequence.
       - Exact FASTA ID match preferred.
       - Fallback to species prefix before first underscore if unique.

    2) Test alternative CDS interpretations:
       - 5' trim lengths from 0..max_trim5.
       - Optional reverse-complement orientation.

    3) Translate each candidate CDS.

    4) Select the best translation:
       - Exact protein match preferred.
       - Otherwise highest pairwise protein identity.

    5) Write corrected CDS sequences and a CSV report describing
       the chosen trim/orientation.

Outputs:
    - Corrected CDS FASTA.
    - CSV report summarizing chosen corrections.

Typical use:
    Repair CDS files after frame-shift artifacts, extra leading bases,
    strand issues, or translation inconsistencies discovered during
    orthology and comparative genomics workflows.
"""
```

Then add these docstrings:

```python
def clean_protein(s: str) -> str:
    """
    Normalize a translated protein sequence.

    Removes whitespace, converts to uppercase, and strips all trailing
    stop codons ('*') from the protein sequence.

    Returns:
        Cleaned protein sequence.
    """
```

```python
def clean_cds(s: str) -> str:
    """
    Normalize a CDS sequence.

    Removes whitespace, converts to uppercase, and replaces any
    non-ACGTN character with 'N'.

    Returns:
        Cleaned nucleotide sequence.
    """
```

```python
def species_key(seq_id: str) -> str:
    """
    Extract a species identifier from a FASTA record ID.

    Uses the first whitespace-delimited token and then keeps everything
    before the first underscore.

    Returns:
        Species key string.
    """
```

```python
def build_aligner() -> PairwiseAligner:
    """
    Create the global protein aligner used for identity comparisons.

    Returns:
        Configured PairwiseAligner object.
    """
```

```python
def align_to_gapped_strings(aligner: PairwiseAligner, s1: str, s2: str) -> Tuple[str, str]:
    """
    Perform a global alignment and reconstruct aligned strings.

    Args:
        aligner: Configured PairwiseAligner.
        s1: Query sequence.
        s2: Target sequence.

    Returns:
        Tuple of aligned sequence strings containing gap characters.
    """
```

```python
def identity_ignore_gaps(aln_a: str, aln_b: str) -> float:
    """
    Calculate protein identity ignoring gap-containing columns.

    Returns:
        Fractional identity from 0.0 to 1.0.
    """
```

```python
def translate_with_trim(cds: Seq, trim5: int, table: int, strand: str) -> Tuple[str, str]:
    """
    Translate a CDS after applying strand selection and 5' trimming.

    Processing steps:
        1) Clean CDS.
        2) Reverse-complement if requested.
        3) Remove trim5 nucleotides from the 5' end.
        4) Trim sequence to a complete codon boundary.
        5) Translate using the requested genetic code.

    Args:
        cds: Input CDS sequence.
        trim5: Number of nucleotides to remove from the 5' end.
        table: NCBI genetic code table number.
        strand: '+' or '-'.

    Returns:
        (translated_protein, processed_cds_sequence)
    """
```

```python
def load_fasta_as_dict(path: str) -> Dict[str, Seq]:
    """
    Load a FASTA file into a dictionary keyed by record ID.

    Returns:
        Dict of record_id -> sequence.
    """
```

```python
def build_prefix_index(ids: List[str]) -> Dict[str, List[str]]:
    """
    Build an index of FASTA IDs by species prefix.

    Returns:
        Dict mapping species_key -> list of FASTA IDs.
    """
```

```python
def resolve_protein_id_for_cds(
    cds_id: str,
    prot_dict: Dict[str, Seq],
    prot_prefix_index: Dict[str, List[str]],
) -> Tuple[Optional[str], str]:
    """
    Determine which protein sequence should be paired with a CDS.

    Matching strategy:
        1) Exact FASTA ID match.
        2) Unique species-prefix match.

    Returns:
        (matched_protein_id, match_mode)

    match_mode:
        EXACT_ID
        PREFIX_ID
        NO_MATCH
        AMBIGUOUS_PREFIX
    """
```

```python
def main():
    """
    Run the CDS frame-correction workflow.

    For each CDS:
        - Find corresponding protein.
        - Test trim and strand combinations.
        - Select best translation match.
        - Write corrected CDS FASTA.
        - Record decisions in CSV report.

    Summary statistics are printed at completion.
    """
```