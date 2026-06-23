```python
#!/usr/bin/env python3

"""
Build a master Orthocaller summary from *_generax summary.txt files.

This script scans all *_generax directories under a base input directory,
extracts orthogroups meeting user-defined species thresholds, and writes a
single consolidated summary file.

For each summary.txt file:
    1) Parse orthogroup entries containing species counts.
    2) Extract total species, cavefish species, and background species.
    3) Apply minimum cavefish and background thresholds.
    4) Rewrite the orthogroup ID using the parent *_generax directory name.
    5) Write passing orthogroups to a master summary file.

Typical use:
    Generate a filtered master orthogroup list for downstream Orthocaller,
    PAL2NAL, selection-test, or sequence-extraction workflows.
"""

import argparse
import os
import re


def parse_args():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace containing input paths and filtering thresholds.
    """
    parser = argparse.ArgumentParser(
        description="Build a master summary from Orthocaller *_generax/summary.txt files."
    )
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Base input directory containing *_generax folders"
    )
    parser.add_argument(
        "-o", "--output_file",
        required=True,
        help="Path to write master summary file"
    )
    parser.add_argument(
        "--min-cavefish",
        type=int,
        default=8,
        help="Minimum number of cavefish species"
    )
    parser.add_argument(
        "--min-background",
        type=int,
        default=31,
        help="Minimum number of background species"
    )
    return parser.parse_args()


def main():
    """
    Build a filtered master orthogroup summary.

    Workflow:
        1) Find all *_generax directories.
        2) Open each summary.txt file.
        3) Parse orthogroup species-count entries.
        4) Filter by cavefish and background thresholds.
        5) Write passing orthogroups to the output file.

    Summary statistics are printed at completion.
    """
    args = parse_args()

    input_dir = os.path.abspath(args.input_dir)
    output_file = os.path.abspath(args.output_file)

    matched_lines = 0
    opened_files = 0
    missing_summary = 0

    summary_pattern = re.compile(
        r"^\s*(\S*Gene-(\d+)):\s+(\d+)\s+species\s+\((\d+)\s+cavefish,\s+(\d+)\s+background\)"
    )

    generax_dirs = sorted(
        d for d in os.listdir(input_dir)
        if d.endswith("_generax")
        and os.path.isdir(os.path.join(input_dir, d))
    )

    with open(output_file, "w", encoding="utf-8") as out:
        for folder_name in generax_dirs:
            summary_path = os.path.join(input_dir, folder_name, "summary.txt")

            if not os.path.isfile(summary_path):
                missing_summary += 1
                continue

            opened_files += 1

            with open(summary_path, "r", encoding="utf-8", errors="replace") as f:
                for raw_line in f:
                    line = raw_line.strip()

                    if "Gene-" not in line or "cavefish" not in line:
                        continue

                    match = summary_pattern.search(line)
                    if not match:
                        continue

                    gene_number = match.group(2)
                    total_species = match.group(3)
                    cavefish = int(match.group(4))
                    background = int(match.group(5))

                    if cavefish < args.min_cavefish:
                        continue

                    if background < args.min_background:
                        continue

                    fixed_line = (
                        f"{folder_name}-Gene-{gene_number}: "
                        f"{total_species} species "
                        f"({cavefish} cavefish, {background} background)"
                    )

                    out.write(fixed_line + "\n")
                    matched_lines += 1

    print("Input directory:", input_dir)
    print("Processed summary files:", opened_files)
    print("Folders missing summary.txt:", missing_summary)
    print("Matching lines written:", matched_lines)
    print("Output written to:", output_file)


if __name__ == "__main__":
    main()
```