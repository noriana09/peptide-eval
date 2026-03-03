import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


def run_pepnet_amp(project_root: Path, fasta_path: Path, out_dir: Path) -> Path:
    """
    Runs PepNet fast-mode AMP prediction and returns the path to PepNet's output CSV.
    """
    pepnet_script_dir = project_root / "tools" / "PepNet" / "script"
    pepnet_predict = pepnet_script_dir / "predict_fast.py"

    if not pepnet_predict.exists():
        raise FileNotFoundError(f"PepNet predict_fast.py not found at: {pepnet_predict}")

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    out_dir.mkdir(parents=True, exist_ok=True)

    # We run from PepNet/script so the relative paths to ../datasets/... work.
    cmd = [
        sys.executable,  # uses the current venv python
        str(pepnet_predict.name),  # run predict_fast.py from cwd
        "--type", "AMP",
        "--test_fasta", str(fasta_path),
        "--output_path", str(out_dir),
    ]

    result = subprocess.run(
        cmd,
        cwd=str(pepnet_script_dir),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    print(result.stdout)

    if result.returncode != 0:
        raise RuntimeError("PepNet failed. See logs above.")

    pepnet_csv = out_dir / "AMP_prediction_result.csv"
    if not pepnet_csv.exists():
        raise FileNotFoundError(f"PepNet finished but output CSV not found: {pepnet_csv}")

    return pepnet_csv


def build_clean_amp_csv(pepnet_csv: Path, fasta_path: Path, final_csv: Path):
    """
    Creates a clean CSV with: id, sequence, AMP_prob, AMP
    id comes from FASTA headers.
    """
    # Parse FASTA ids in order
    ids = []
    seqs = []
    cur_id = None
    cur_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    ids.append(cur_id)
                    seqs.append("".join(cur_seq))
                cur_id = line[1:].split()[0]  # take first token after '>'
                cur_seq = []
            else:
                cur_seq.append(line)

    if cur_id is not None:
        ids.append(cur_id)
        seqs.append("".join(cur_seq))

    # Read PepNet output
    df = pd.read_csv(pepnet_csv)

    # PepNet CSV usually includes an unnamed index col; keep only what we need
    if "sequence" not in df.columns or "probability" not in df.columns or "Binary" not in df.columns:
        raise ValueError(f"Unexpected PepNet CSV columns: {list(df.columns)}")

    # Sanity: match sequences length
    if len(df) != len(seqs):
        raise ValueError(f"FASTA has {len(seqs)} sequences but PepNet returned {len(df)} rows.")

    # Final clean dataframe
    out = pd.DataFrame({
        "id": ids,
        "sequence": df["sequence"].astype(str),
        "AMP_prob": df["probability"].astype(float),
        "AMP": df["Binary"].astype(int),
    })

    final_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(final_csv, index=False)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, help="Path to input FASTA")
    parser.add_argument("--out_dir", default="outputs/pepnet", help="Directory for PepNet raw outputs")
    parser.add_argument("--final_csv", default="outputs/results_amp.csv", help="Final clean CSV output path")
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parents[1]
    fasta_path = (project_root / args.fasta).resolve() if not os.path.isabs(args.fasta) else Path(args.fasta)
    out_dir = (project_root / args.out_dir).resolve() if not os.path.isabs(args.out_dir) else Path(args.out_dir)
    final_csv = (project_root / args.final_csv).resolve() if not os.path.isabs(args.final_csv) else Path(args.final_csv)

    pepnet_csv = run_pepnet_amp(project_root, fasta_path, out_dir)
    out_df = build_clean_amp_csv(pepnet_csv, fasta_path, final_csv)

    print(f"\nWrote clean AMP results to: {final_csv}")
    print(out_df.head())


if __name__ == "__main__":
    main()