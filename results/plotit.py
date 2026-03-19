import os
import re
import shutil
import tarfile
import glob

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


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


def plot_and_save(csv_path, prefix):
    """Run the two plots and save SVGs named <prefix>_runtime_and_ram.svg / <prefix>_total_io.svg."""
    df = pd.read_csv(csv_path)

    for col in ['s', 'max_rss', 'io_in', 'io_out']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df['total_io'] = df['io_in'] + df['io_out']

    # --- PLOT 1: Performance (Runtime & RAM) ---
    fig1, axes1 = plt.subplots(1, 2, figsize=(14, 6))

    sns.lineplot(data=df, x='reads', y='s', hue='tool', marker='o', ax=axes1[0])
    axes1[0].set_title(f'Execution Time vs Number of Reads')
    axes1[0].set_ylabel('Time (s)')
    axes1[0].set_xscale('log')
    axes1[0].set_yscale('log')

    sns.lineplot(data=df, x='reads', y='max_rss', hue='tool', marker='o', ax=axes1[1])
    axes1[1].set_title(f'Max Memory (RSS) vs Number of Reads')
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
    sns.lineplot(data=df, x='reads', y='total_io', hue='tool', marker='o', ax=ax2)
    ax2.set_title(f'Total I/O vs Number of Reads')
    ax2.set_ylabel('Total I/O (in + out)')
    ax2.set_xlabel('Reads')
    ax2.set_xscale('log')
    ax2.set_yscale('symlog', linthresh=1.0)

    fig2.tight_layout()
    out2 = f'{prefix}_total_io.svg'
    fig2.savefig(out2)
    print(f'  Saved: {out2}')
    plt.close(fig2)


def process_archive(archive_path):
    archive_name = os.path.basename(archive_path)
    print(f'\n=== Processing: {archive_name} ===')

    mode = detect_mode(archive_name)

    # Extract into a temp dir so we never collide with an existing "results" folder
    tmp_dir = f'_tmp_extract_{os.getpid()}'
    os.makedirs(tmp_dir, exist_ok=True)

    try:
        print(f'  Extracting {archive_name} → {tmp_dir}/')
        with tarfile.open(archive_path, 'r:gz') as tar:
            tar.extractall(tmp_dir)

        # Find the results directory (could be named 'results' or similar)
        extracted_dirs = [
            os.path.join(tmp_dir, d)
            for d in os.listdir(tmp_dir)
            if os.path.isdir(os.path.join(tmp_dir, d))
        ]
        if not extracted_dirs:
            print('  ERROR: no directory found after extraction, skipping.')
            return

        results_dir = extracted_dirs[0]   # e.g. _tmp_.../results
        print(f'  Detected results dir: {os.path.basename(results_dir)}')

        layout = detect_layout(results_dir)
        prefix = f'{layout}_{mode}'
        print(f'  Detected prefix: {prefix}')

        # Find all_reports.csv (handles both all_reports.csv and all_report.csv)
        csv_candidates = glob.glob(os.path.join(results_dir, 'all_report*.csv'))
        if not csv_candidates:
            print('  ERROR: all_reports.csv not found, skipping.')
            return

        src_csv = csv_candidates[0]
        dst_csv = f'{prefix}.csv'
        shutil.copy2(src_csv, dst_csv)
        print(f'  CSV saved as: {dst_csv}')

        plot_and_save(dst_csv, prefix)

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        print(f'  Cleaned up temp dir: {tmp_dir}')


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