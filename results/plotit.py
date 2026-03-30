import os
import re
import shutil
import tarfile
import glob
import subprocess
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def detect_layout(results_dir):
    """Infer 'pe' (paired) or 'se' (single) from subfolder names inside results_dir."""
    for entry in os.listdir(results_dir):
        if os.path.isdir(os.path.join(results_dir, entry)):
            name = entry.lower()
            if "paired" in name:
                return "pe"
            if "single" in name:
                return "se"
    return "unknown"


def detect_mode(archive_name):
    """Infer 'dd', 'cd', or 'cc' from the archive filename."""
    match = re.search(r'(dd|cd|cc)', archive_name.lower())
    return match.group(1) if match else "unknown"


def estimate_natural_duplication(num_reads: int) -> float:
    """
    Estimates the natural duplication count and rate for reads generated
    by wgsim from the GRCh38 reference genome.

    Args:
        num_reads (int): The total number of unique base reads generated.

    Returns:
        tuple: (estimated_duplicate_count, duplication_rate_percentage)
    """
    if num_reads <= 0:
        return 0.0

    # Empirical constant 'k' derived from: C = k * N^2
    # Based on the observation of ~60,000 duplicates at 270,000,000 reads.
    k = 60_000 / (270_000_000**2)

    # Expected duplicates scale quadratically with dataset size
    expected_duplicates = int(k * (num_reads ** 2))

    # Calculate the percentage rate
    duplication_rate_pct = (expected_duplicates / num_reads) * 100

    return duplication_rate_pct


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def build_palette(tool_lists):
    """
    Given a list of iterables of tool names, return a single dict {tool: color}
    with a perfectly stable assignment across all plots and different archives.
    """
    current_tools = {t for tools in tool_lists for t in tools}

    # Master list of tools — exact names as they appear in the CSVs.
    # Order is fixed so colors never shift between plots or archives.
    master_tools = [
        'fdedup_h64',
        'fdedup_h128',
        'fastp',
        'fastuniq',
        'clumpify',
        'seqkit',
        'prinseqpp',
        'cdhitdup',
    ]

    # Deduplicate the master list while preserving the order
    all_tools = list(dict.fromkeys(master_tools))

    # Append any unforeseen tools found in the CSV dynamically
    unexpected_tools = sorted(list(current_tools - set(all_tools)))
    all_tools.extend(unexpected_tools)

    # Use 'tab20' to ensure we have enough distinct colors for all variations
    colors = sns.color_palette('tab20', n_colors=max(20, len(all_tools)))

    # Map each tool to its permanently assigned color
    return {tool: colors[i] for i, tool in enumerate(all_tools)}


def plot_and_save(csv_path, prefix, dup_rates_dir=None):
    """
    Run all three plots and save SVGs:
      <prefix>_runtime_and_ram.svg
      <prefix>_total_io.svg
      <prefix>_duplication_rate.svg   (if dup_rates_dir is provided)
    """
    df = pd.read_csv(csv_path)

    is_compressed = prefix.endswith('_cc') or prefix.endswith('_cd')
    if is_compressed:
        df = df[~df['tool'].astype(str).str.lower().isin(['cdhitdup'])]
        df.to_csv(csv_path, index=False)

    for col in ['s', 'max_rss', 'io_in', 'io_out']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df['total_io'] = df['io_in'] + df['io_out']
    df['reads_abs'] = df['reads'] * 10**6

    # --- Load duplication data early so we can include its tools in the shared palette ---
    agg = None
    if dup_rates_dir is not None:
        tsv_files = glob.glob(os.path.join(dup_rates_dir, '*', 'duplication_rates.tsv'))
        if not tsv_files:
            print('  WARNING: no duplication_rates.tsv files found, skipping plot 3.')
        else:
            frames = []
            for tsv in tsv_files:
                try:
                    frames.append(pd.read_csv(tsv, sep='\t'))
                except Exception as e:
                    print(f'  WARNING: could not read {tsv}: {e}')
            if frames:
                dup_df = pd.concat(frames, ignore_index=True)

                if is_compressed:
                    dup_df = dup_df[~dup_df['tool'].astype(str).str.lower().isin(['cdhitdup'])]

                out_csv = f'{prefix}_duplication.csv'
                dup_df.to_csv(out_csv, index=False)
                print(f'  Saved: {out_csv}')

                dup_df['tool_rate'] = pd.to_numeric(dup_df['tool_rate'], errors='coerce')
                dup_df['reads'] = pd.to_numeric(dup_df['reads'], errors='coerce')
                dup_df = dup_df.dropna(subset=['tool_rate', 'reads'])
                agg = (
                    dup_df.groupby(['reads', 'tool'])['tool_rate']
                    .agg(['mean', 'std'])
                    .reset_index()
                )
                agg.columns = ['reads', 'tool', 'mean_rate', 'std_rate']
                agg['std_rate'] = agg['std_rate'].fillna(0)

    # Build one shared palette from the union of all tools across all three plots
    perf_tools = df['tool'].dropna().unique().tolist()
    dup_tools = agg['tool'].unique().tolist() if agg is not None else []
    palette = build_palette([perf_tools, dup_tools])

    # --- PLOT 1: Performance (Runtime & RAM) ---
    fig1, axes1 = plt.subplots(1, 2, figsize=(14, 6))

    sns.lineplot(data=df, x='reads_abs', y='s', hue='tool', marker='o',
                 palette=palette, ax=axes1[0])
    axes1[0].set_title('Execution Time vs Number of Reads')
    axes1[0].set_ylabel('Time (s)')
    axes1[0].set_xlabel('Reads')
    axes1[0].set_xscale('log')
    axes1[0].set_yscale('log')

    sns.lineplot(data=df, x='reads_abs', y='max_rss', hue='tool', marker='o',
                 palette=palette, ax=axes1[1])
    axes1[1].set_title('Max Memory (RSS) vs Number of Reads')
    axes1[1].set_xlabel('Reads')
    axes1[1].set_ylabel('Max RSS')
    axes1[1].set_xscale('log')
    axes1[1].set_yscale('log')

    fig1.tight_layout()
    out1 = f'{prefix}_runtime_and_ram.svg'
    fig1.savefig(out1)
    print(f'  Saved: {out1}')
    plt.close(fig1)

    # --- PLOT 2: Throughput (Total I/O) ---
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    sns.lineplot(data=df, x='reads_abs', y='total_io', hue='tool', marker='o',
                 palette=palette, ax=ax2)
    ax2.set_title('Total I/O vs Number of Reads')
    ax2.set_ylabel('Total I/O (in + out)')
    ax2.set_xlabel('Reads')
    ax2.set_xscale('log')
    ax2.set_yscale('symlog', linthresh=1.0)

    fig2.tight_layout()
    out2 = f'{prefix}_total_io.svg'
    fig2.savefig(out2)
    print(f'  Saved: {out2}')
    plt.close(fig2)

    # --- PLOT 3: Duplication Rate ---
    if agg is None:
        return

    fig3, ax3 = plt.subplots(figsize=(10, 6))

    # Plot each tool's lines and standard deviation fill
    for tool in sorted(agg['tool'].unique()):
        color = palette[tool]
        sub = agg[agg['tool'] == tool].sort_values('reads')
        ax3.plot(sub['reads']* 10**6, sub['mean_rate'], marker='o', label=tool, color=color)
        ax3.fill_between(
            sub['reads'] * 10**6,
            sub['mean_rate'] - sub['std_rate'],
            sub['mean_rate'] + sub['std_rate'],
            alpha=0.15,
            color=color,
        )

    # ---------------------------------------------------------
    # Target Line & Tolerance Band (Outside the tool loop)
    # ---------------------------------------------------------
    unique_reads = sorted(agg['reads'].unique())
    unique_reads = [u * 10 ** 6 for u in unique_reads]
    is_single_end = prefix.startswith('se')

    lower_bounds = []
    upper_bounds = []
    target_rates = []

    for r in unique_reads:
        # Only apply the natural duplication math if this is a single-end run
        nat_dup = estimate_natural_duplication(r) if is_single_end else 0.0

        lower_bounds.append(9.99 + nat_dup)
        upper_bounds.append(10.01 + nat_dup)
        target_rates.append(10.0 + nat_dup)

    # 1. Draw the top, bottom, and target lines
    ax3.plot(unique_reads, upper_bounds, color='grey', linewidth=1, linestyle='--')
    ax3.plot(unique_reads, lower_bounds, color='grey', linewidth=1, linestyle='--')
    ax3.plot(unique_reads, target_rates, color='black', linewidth=1.2, linestyle='--', label='Target (10% + natural)')

    # 2. Fill the tolerance band
    ax3.fill_between(
        unique_reads,
        lower_bounds,
        upper_bounds,
        color='grey',
        alpha=0.12
    )
    # ---------------------------------------------------------

    ax3.set_xscale('log')

    # Adjust Y-limits dynamically so single-end graphs don't clip the top of the band
    ax3.set_ylim(9.98, 10.05)
    ax3.set_xlim(min(unique_reads), max(unique_reads))

    ax3.set_xlabel('Reads')
    ax3.set_ylabel('Detected Duplication Rate (%)')
    ax3.set_title(f'Duplication Rate Accuracy vs Number of Reads')
    ax3.legend(loc='upper left', fontsize=8)
    ax3.grid(True, which='both', linestyle=':', linewidth=0.5, alpha=0.7)

    fig3.tight_layout()
    out3 = f'{prefix}_duplication_rate.svg'
    fig3.savefig(out3)
    print(f'  Saved: {out3}')
    plt.close(fig3)

# ---------------------------------------------------------------------------
# Archive processing
# ---------------------------------------------------------------------------

def run_extraction(results_dir):
    """Run extract_rates_standalone.py on the given results dir if it exists."""
    script = 'extract_rates_standalone.py'
    if not os.path.exists(script):
        print(f'  WARNING: {script} not found — skipping duplication rate extraction.')
        return False
    print(f'  Running {script} on {results_dir} …')
    ret = subprocess.run(
        [sys.executable, script, results_dir],
        capture_output=True, text=True
    )
    if ret.stdout:
        for line in ret.stdout.strip().splitlines():
            print(f'    {line}')
    if ret.returncode != 0:
        print(f'  WARNING: extraction script exited with code {ret.returncode}')
        if ret.stderr:
            print(ret.stderr[:500])
        return False
    return True


def process_archive(archive_path):
    archive_name = os.path.basename(archive_path)
    print(f'\n=== Processing: {archive_name} ===')

    mode = detect_mode(archive_name)

    tmp_dir = f'_tmp_extract_{os.getpid()}'
    os.makedirs(tmp_dir, exist_ok=True)

    try:
        print(f'  Extracting {archive_name} → {tmp_dir}/')
        with tarfile.open(archive_path, 'r:gz') as tar:
            tar.extractall(tmp_dir)

        extracted_dirs = [
            os.path.join(tmp_dir, d)
            for d in os.listdir(tmp_dir)
            if os.path.isdir(os.path.join(tmp_dir, d))
        ]
        if not extracted_dirs:
            print('  ERROR: no directory found after extraction, skipping.')
            return

        results_dir = extracted_dirs[0]
        print(f'  Detected results dir: {os.path.basename(results_dir)}')

        layout = detect_layout(results_dir)
        prefix = f'{layout}_{mode}'
        print(f'  Detected prefix: {prefix}')

        # Find and copy the performance CSV
        csv_candidates = glob.glob(os.path.join(results_dir, 'all_report*.csv'))
        if not csv_candidates:
            print('  ERROR: all_reports.csv not found, skipping.')
            return

        src_csv = csv_candidates[0]
        dst_csv = f'{prefix}.csv'
        shutil.copy2(src_csv, dst_csv)
        print(f'  CSV saved as: {dst_csv}')

        # Re-extract duplication rates from the unpacked logs
        run_extraction(results_dir)

        plot_and_save(dst_csv, prefix, dup_rates_dir=results_dir)

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        print(f'  Cleaned up temp dir: {tmp_dir}')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    archives = sorted(glob.glob('*.tar.gz'))
    if not archives:
        print('No .tar.gz archives found in the current directory.')
        return

    print(f'Found {len(archives)} archive(s): {archives}')
    for archive in archives:
        process_archive(archive)

    print('\nAll done.')


if __name__ == '__main__':
    main()