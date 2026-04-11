#!/usr/bin/env python3

"""
multi_species_rnaseq_compare.py
-----
Compare RNAseq expression of a gene of interest using raw count files
and a metadata table.  Two comparison modes are available:

  within   – compare two tissues within each species (one panel per species)
             e.g. liver vs gonads for Crotalus adamanteus, Crotalus tigris,Thamnophis sirtalis …

  across   – compare one tissue across all species (species on the x-axis)
             e.g. gonad expression in Crotalus adamanteus vs Crotalus tigris vs Thamnophis sirtalis

Usage
-----
# Within-species (tissue1 vs tissue2 for every species):
python3 multi_species_rnaseq_compare.py \\
    --counts *.tsv --metadata metadata.csv \\
    --gene "PRDM9" --mode within \\
    --tissue1 liver --tissue2 brain \\
    --output results.csv --plot results.png

# Across-species (one tissue, all species):
python3 multi_species_rnaseq_compare.py \\
    --counts *.tsv --metadata metadata.csv \\
    --gene "PRDM9" --mode across \\
    --tissue1 liver \\
    --output results.csv --plot results.png

# Across-species, restrict to a subset of species:
python3 multi_species_rnaseq_compare.py \\
    --counts *.tsv --metadata metadata.csv \\
    --gene "PRDM9" --mode across \\
    --tissue1 gonads --species "Crotalus adamanteus" "Crotalus tigris" "Thamnophis sirtalis" \\
    --output results.csv
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


# Loads and validates the metadata file, ensuring required columns are present and 
# normalised to lower-case for consistent access. 
# The sample_id column is also stripped of whitespace and converted

def load_metadata(metadata_path: str) -> pd.DataFrame:
    """
    Load the metadata file (CSV or TSV).

    Required columns (case-insensitive):
        species   – species name
        sample_id – SRR / sample accession matching count-file column headers
        tissue    – tissue label (e.g. gonad, liver,skin, etc.)

    Returns a DataFrame with columns normalised to lower-case.
    """
    path = Path(metadata_path)
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep)
    df.columns = df.columns.str.strip().str.lower()

    required = {"species", "sample_id", "tissue"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Metadata file is missing required column(s): {missing}\n"
            f"Found columns: {list(df.columns)}"
        )

    df["sample_id"] = df["sample_id"].astype(str).str.strip()
    return df


# Loads a count file, ensuring the first column is 'gene_id' and setting it as the index.
def load_count_file(filepath: str) -> pd.DataFrame:
    """
    Load a single count file. Auto-detects comma vs tab separator
    based on file extension (.csv → comma, .tsv/.txt → tab).
    """
    path = Path(filepath)
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    df.index = df.index.astype(str)
    df.index.name = "gene_id"
    return df


# Search the gene_id index of a count DataFrame for rows matching a regex pattern.

def find_gene_rows(count_df: pd.DataFrame, pattern: str) -> pd.DataFrame:
    """
    Search the gene_id index of *count_df* with *pattern* (regex,
    case-insensitive). Returns all matching rows.
    """
    regex = re.compile(pattern, re.IGNORECASE)
    matched = [idx for idx in count_df.index if regex.search(idx)]
    return count_df.loc[matched] if matched else pd.DataFrame()

# If multiple rows match the gene pattern, prompt the user to select one.
def resolve_gene(count_df: pd.DataFrame, pattern: str) -> tuple[pd.Series, str]:
    """
    Return (row Series, gene_id string) for the gene matching *pattern*.

    Prompts the user to pick if multiple rows match.
    Exits with an error if nothing matches.
    """
    hits = find_gene_rows(count_df, pattern)

    if hits.empty:
        print(f"[ERROR] No gene matching '{pattern}' found.", file=sys.stderr)
        sys.exit(1)

    if len(hits) == 1:
        return hits.iloc[0], hits.index[0]

    print(f"\n[INFO] Multiple genes matched '{pattern}':")
    for i, gene_id in enumerate(hits.index):
        print(f"  [{i}] {gene_id}")
    while True:
        try:
            choice = int(input("Enter the number of the gene to use: "))
            if 0 <= choice < len(hits):
                return hits.iloc[choice], hits.index[choice]
        except ValueError:
            pass
        print("  Invalid choice, try again.")


# Extract expression values for the gene of interest from all count files,
# filtering by tissue and optionally by species. Returns a tidy DataFrame.

def extract_expression(
    count_files: list[str],
    metadata: pd.DataFrame,
    gene_pattern: str,
    tissues: list[str],
    species_filter: list[str] | None = None,
) -> pd.DataFrame:
    """
    Collect expression values for *gene_pattern* from all count files.

    Parameters
    ----------
    count_files    : paths to count TSV files
    metadata       : loaded metadata DataFrame
    gene_pattern   : regex pattern for the gene_id column
    tissues        : list of tissue labels to retain
    species_filter : optional list of species to retain (None = all)

    Returns a tidy DataFrame: gene_id | species | tissue | sample_id | count
    """
    tissues_lower = [t.lower() for t in tissues]
    meta_sub = metadata[metadata["tissue"].str.lower().isin(tissues_lower)].copy()

    if species_filter:
        species_lower = [s.lower() for s in species_filter]
        meta_sub = meta_sub[meta_sub["species"].str.lower().isin(species_lower)]

    if meta_sub.empty:
        print(
            "[ERROR] Metadata filter returned no samples. "
            "Check tissue and species names.",
            file=sys.stderr,
        )
        sys.exit(1)

    records = []
    gene_label = None

    for filepath in count_files:
        print(f"[INFO] Processing: {filepath}")
        try:
            counts = load_count_file(filepath)
        except Exception as exc:
            print(f"[WARNING] Skipping '{filepath}': {exc}", file=sys.stderr)
            continue

        hits = find_gene_rows(counts, gene_pattern)
        if hits.empty:
            print(f"[WARNING] '{gene_pattern}' not found in {filepath} — skipping.")
            continue
        row, matched_gene_id = resolve_gene(counts, gene_pattern)

        if gene_label is None:
            gene_label = matched_gene_id

        available = set(counts.columns)
        for _, meta_row in meta_sub[meta_sub["sample_id"].isin(available)].iterrows():
            sid = meta_row["sample_id"]
            records.append(
                {
                    "gene_id": matched_gene_id,
                    "species": meta_row["species"],
                    "tissue": meta_row["tissue"].lower(),
                    "sample_id": sid,
                    "count": row[sid],
                }
            )

    if not records:
        print(
            "[ERROR] No expression data collected. "
            "Verify sample IDs match count-file column headers.",
            file=sys.stderr,
        )
        sys.exit(1)

    return pd.DataFrame(records)


# Summarise expression by species × tissue, calculating mean and SD of counts.

def summarise_expression(expr_df: pd.DataFrame) -> pd.DataFrame:
    """Mean ± SD of counts grouped by species × tissue."""
    summary = (
        expr_df.groupby(["species", "tissue"])["count"]
        .agg(n="count", mean="mean", std="std")
        .reset_index()
    )
    summary["std"] = summary["std"].fillna(0)
    return summary


# Plotting – within mode: tissue1 vs tissue2 for every species (one panel per species).

def plot_within(
    expr_df: pd.DataFrame,
    gene_label: str,
    tissue1: str,
    tissue2: str,
    output_plot: str | None = None,
) -> None:
    """
    Within-species comparison: one grouped pair of bars per species,
    coloured by tissue (tissue1 vs tissue2).

    Each species panel shows tissue1 and tissue2 side-by-side so you can
    immediately read the within-species difference, while the x-axis lets
    you compare across species at the same time.
    """
    summary = summarise_expression(expr_df)
    species_list = sorted(expr_df["species"].unique())
    t1, t2 = tissue1.lower(), tissue2.lower()
    colours = {"t1": "steelblue", "t2": "darkorange"}

    x = np.arange(len(species_list))
    bar_width = 0.35
    fig, ax = plt.subplots(figsize=(max(6, len(species_list) * 1.8), 5))

    for i, (tissue, colour, label) in enumerate(
        [(t1, colours["t1"], tissue1), (t2, colours["t2"], tissue2)]
    ):
        sub = (
            summary[summary["tissue"] == tissue]
            .set_index("species")
        )
        means = [sub.loc[sp, "mean"] if sp in sub.index else 0 for sp in species_list]
        stds  = [sub.loc[sp, "std"]  if sp in sub.index else 0 for sp in species_list]
        offset = x + (i - 0.5) * bar_width
        ax.bar(offset, means, bar_width, yerr=stds,
               label=label.capitalize(), color=colour, capsize=4, alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(species_list, rotation=30, ha="right")
    ax.set_ylabel("Read count (mean ± SD)")
    ax.set_title(f"{gene_label} — within-species: {tissue1} vs {tissue2}")
    ax.legend(title="Tissue")
    plt.tight_layout()
    _save_or_show(output_plot)


# Plotting – across mode: species on the x-axis, one (or two) tissue(s) compared across all species.

def plot_across(
    expr_df: pd.DataFrame,
    gene_label: str,
    tissues: list[str],
    output_plot: str | None = None,
) -> None:
    """
    Across-species comparison: species on the x-axis.

    If a single tissue was requested, bars are coloured by species.
    If two tissues were requested, bars are grouped by tissue (same as
    within-mode layout but the intent is inter-species, not intra-species).
    """
    summary = summarise_expression(expr_df)
    species_list = sorted(expr_df["species"].unique())
    tissues_lower = [t.lower() for t in tissues]

    bar_width = 0.35
    x = np.arange(len(species_list))
    fig, ax = plt.subplots(figsize=(max(6, len(species_list) * 1.8), 5))

    if len(tissues_lower) == 1:
        # Single tissue – one bar per species, coloured distinctly
        tissue = tissues_lower[0]
        sub = summary[summary["tissue"] == tissue].set_index("species")
        colours = cm.tab10(np.linspace(0, 0.9, len(species_list)))
        means = [sub.loc[sp, "mean"] if sp in sub.index else 0 for sp in species_list]
        stds  = [sub.loc[sp, "std"]  if sp in sub.index else 0 for sp in species_list]
        bars = ax.bar(x, means, bar_width * 1.4, yerr=stds,
                      color=colours, capsize=4, alpha=0.85)
        ax.set_title(
            f"{gene_label} — across species: {tissues[0].capitalize()}"
        )
        ax.legend(bars, species_list, title="Species",
                  bbox_to_anchor=(1.01, 1), loc="upper left")

    else:
        # Two tissues requested in across mode – grouped bars
        palette = ["steelblue", "darkorange"]
        for i, (tissue, colour) in enumerate(zip(tissues_lower, palette)):
            sub = summary[summary["tissue"] == tissue].set_index("species")
            means = [sub.loc[sp, "mean"] if sp in sub.index else 0 for sp in species_list]
            stds  = [sub.loc[sp, "std"]  if sp in sub.index else 0 for sp in species_list]
            offset = x + (i - 0.5) * bar_width
            ax.bar(offset, means, bar_width, yerr=stds,
                   label=tissues[i].capitalize(), color=colour, capsize=4, alpha=0.85)
        ax.legend(title="Tissue")
        ax.set_title(
            f"{gene_label} — across species: "
            f"{tissues[0].capitalize()} vs {tissues[1].capitalize()}"
        )

    ax.set_xticks(x)
    ax.set_xticklabels(species_list, rotation=30, ha="right")
    ax.set_ylabel("Read count (mean ± SD)")
    plt.tight_layout()
    _save_or_show(output_plot)


def _save_or_show(output_plot: str | None) -> None:
    if output_plot:
        plt.savefig(output_plot, dpi=150)
        print(f"[INFO] Plot saved to: {output_plot}")
    else:
        plt.show()


# Save the tidy expression DataFrame to a CSV file.

def save_results(expr_df: pd.DataFrame, output_path: str) -> None:
    """Write the tidy expression table to CSV."""
    expr_df.to_csv(output_path, index=False)
    print(f"[INFO] Results saved to: {output_path}")


# Build the command-line argument parser with detailed help messages and validation rules.

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq_compare",
        description=(
            "Compare RNAseq expression of a gene across tissues and/or species.\n\n"
            "Modes\n"
            "-----\n"
            "  within  – tissue1 vs tissue2 for every species\n"
            "  across  – one (or two) tissue(s) compared across all species"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "-c", "--counts",
        nargs="+", required=True, metavar="COUNT_FILE",
        help="One or more tab-separated count files.",
    )
    parser.add_argument(
        "-m", "--metadata",
        required=True, metavar="METADATA_FILE",
        help="CSV or TSV metadata file with columns: species, sample_id, tissue.",
    )
    parser.add_argument(
        "-g", "--gene",
        required=True, metavar="PATTERN",
        help="Regex to search the gene_id column (e.g. 'PRDM9').",
    )

    # Comparison mode
    parser.add_argument(
        "--mode",
        choices=["within", "across"],
        required=True,
        help=(
            "'within': compare tissue1 vs tissue2 inside each species.  "
            "'across': compare one tissue (or two) across species."
        ),
    )

    # Tissue selection
    parser.add_argument(
        "--tissue1",
        required=True,
        help="Primary tissue label (required for both modes).",
    )
    parser.add_argument(
        "--tissue2",
        default=None,
        help=(
            "Second tissue label. "
            "Required for 'within' mode; optional for 'across'."
        ),
    )

    # Species filter (across mode)
    parser.add_argument(
        "--species",
        nargs="+", default=None, metavar="SPECIES",
        help=(
            "Restrict 'across' mode to these species "
            "(default: all species in metadata). "
            "Use quotes for names with spaces, e.g. --species 'Crotalus tigris' 'Crotalus adamanteus'."
        ),
    )

    # Output options
    parser.add_argument(
        "-o", "--output",
        default="expression_results.csv",
        help="Path for the output CSV file (default: expression_results.csv).",
    )
    parser.add_argument(
        "--plot",
        metavar="PLOT_FILE", default=None,
        help="Save the chart to this file instead of displaying it interactively.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip plotting entirely.",
    )

    return parser


def validate_args(args: argparse.Namespace) -> None:
    """Cross-argument validation."""
    if args.mode == "within" and args.tissue2 is None:
        print(
            "[ERROR] --tissue2 is required when --mode within is used.",
            file=sys.stderr,
        )
        sys.exit(1)
    if args.mode == "across" and args.species and args.tissue2 is None:
        # Fine – user wants one tissue across a subset of species
        pass


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args)

    # Load metadata
    print(f"[INFO] Loading metadata from: {args.metadata}")
    metadata = load_metadata(args.metadata)
    print(f"[INFO] {len(metadata)} samples found in metadata.")

    # Determine tissues to collect
    tissues = [args.tissue1]
    if args.tissue2:
        tissues.append(args.tissue2)

    # Extract expression
    expr_df = extract_expression(
        count_files=args.counts,
        metadata=metadata,
        gene_pattern=args.gene,
        tissues=tissues,
        species_filter=args.species if args.mode == "across" else None,
    )

    gene_label = expr_df["gene_id"].iloc[0]
    print(f"\n[INFO] Gene matched : {gene_label}")
    print(f"[INFO] Mode         : {args.mode}")
    print(f"[INFO] Observations : {len(expr_df)} sample-level records\n")

    # Summary table (printed to console)
    summary = summarise_expression(expr_df)
    print(">>>>>Expression Summary (mean ± SD)<<<<<<<")
    print(summary.to_string(index=False))
    print()

    # Save tidy results
    save_results(expr_df, args.output)

    # ---- Plot ----
    if not args.no_plot:
        if args.mode == "within":
            plot_within(
                expr_df=expr_df,
                gene_label=gene_label,
                tissue1=args.tissue1,
                tissue2=args.tissue2,
                output_plot=args.plot,
            )
        else:
            plot_across(
                expr_df=expr_df,
                gene_label=gene_label,
                tissues=tissues,
                output_plot=args.plot,
            )


if __name__ == "__main__":
    main()
