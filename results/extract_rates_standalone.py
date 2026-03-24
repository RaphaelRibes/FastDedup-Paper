import re
import csv
import os
import sys

RESULTS_DIR = sys.argv[1] if len(sys.argv) > 1 else "results"

def strip_ansi(text):
    return re.sub(r'\x1b\[[0-9;]*[mGKHF]', '', text)

def extract_rate(tool, tool_log, total_fragments, mode):
    tool_rate = "NaN"
    try:
        with open(tool_log, 'r') as f:
            raw_content = f.read()
    except FileNotFoundError:
        return tool_rate

    content = strip_ansi(raw_content)
    lines = content.split('\n')

    if tool == "fastp" or tool.startswith("fdedup"):
        match = re.search(r"(?:Duplication rate|Duplicates removed)[^:]*:\s*([\d.]+)%", content)
        if match:
            tool_rate = match.group(1)

    elif tool == "clumpify":
        matches = re.findall(r"Duplicates Found:\s+\d+\s+([\d.]+)%", content)
        if matches:
            tool_rate = matches[-1]

    elif tool == "seqkit":
        match = re.search(r"\[INFO\]\s+(\d+)\s+duplicated records removed", content)
        if match:
            removed = int(match.group(1))
            tool_rate = f"{(removed / total_fragments) * 100:.3f}"

    elif tool == "prinseqpp":
        # Extract numeric lines from log
        num_lines = [line.strip() for line in lines if line.strip().isdigit()]
        if len(num_lines) >= 9:
            removed = int(num_lines[8])
            # If paired, prinseq counts R1 and R2 separately, so we compare against total reads
            denominator = total_fragments * 2 if mode == "paired" else total_fragments
            if removed > 0 and denominator > 0:
                tool_rate = f"{(removed / denominator) * 100:.3f}"

    elif tool == "cdhitdup":
        reads_match = re.search(r"Number of reads:\s*(\d+)", content)
        clusters_match = re.search(r"Number of clusters found:\s*(\d+)", content)
        if reads_match and clusters_match:
            reads = int(reads_match.group(1))
            clusters = int(clusters_match.group(1))
            if reads > 0:
                tool_rate = f"{((reads - clusters) / reads) * 100:.3f}"

    elif tool == "fastuniq":
        matches = re.findall(r"Input Reads:\s*(\d+).*?Unique Reads:\s*(\d+)", content, re.DOTALL)
        if matches:
            total, unique = map(int, matches[-1])
            if total > 0:
                tool_rate = f"{((total - unique) / total) * 100:.3f}"

    return tool_rate

def discover_tools_and_runs(condition_dir):
    tools = []
    max_run = 0
    for entry in sorted(os.listdir(condition_dir)):
        entry_path = os.path.join(condition_dir, entry)
        if not os.path.isdir(entry_path):
            continue
        run_dirs = [d for d in os.listdir(entry_path) if re.match(r"run_\d+$", d)]
        if run_dirs:
            tools.append(entry)
            current_max = max(int(d.split("_")[1]) for d in run_dirs)
            max_run = max(max_run, current_max)
    return tools, max_run

def process_condition(condition_dir, condition_name):
    mreads_match = re.match(r"^([\d.]+)M_", condition_name)
    if not mreads_match:
        return

    mreads = mreads_match.group(1)
    total_fragments = float(mreads) * 1_000_000

    parts = condition_name.split("_")
    paired_mode = parts[1] if len(parts) >= 2 else "single" # Default to single if unknown
    in_out = parts[2] if len(parts) >= 3 else "unknown"

    tools, runs = discover_tools_and_runs(condition_dir)
    if not tools:
        return

    print(f"  Processing {condition_name} ({paired_mode})...")

    results = []
    for tool in tools:
        for run in range(1, runs + 1):
            log_path = os.path.join(condition_dir, tool, f"run_{run}", "log.txt")
            # Pass paired_mode to the extractor
            rate = extract_rate(tool, log_path, total_fragments, paired_mode)
            results.append({
                "reads": mreads,
                "mode": paired_mode,
                "in_out": in_out,
                "tool": tool,
                "run": run,
                "tool_rate": rate,
            })

    out_tsv = os.path.join(condition_dir, "duplication_rates.tsv")
    with open(out_tsv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["reads", "mode", "in_out", "tool", "run", "tool_rate"], delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

def main():
    if not os.path.isdir(RESULTS_DIR):
        print(f"Error: {RESULTS_DIR} not found.")
        sys.exit(1)

    conditions = sorted([d for d in os.listdir(RESULTS_DIR) if re.match(r"[\d.]+M_", d)])
    for condition in conditions:
        process_condition(os.path.join(RESULTS_DIR, condition), condition)
    print("\nExtraction complete.")

if __name__ == "__main__":
    main()