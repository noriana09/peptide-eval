import argparse
import subprocess
import sys
from pathlib import Path

def run(cmd):
    print("\n>>> Running:", " ".join(cmd))
    res = subprocess.run(cmd, check=False)
    if res.returncode != 0:
        print(f"\nCommand failed with code {res.returncode}")
        sys.exit(res.returncode)

def must_exist(path: Path, label: str):
    if not path.exists():
        print(f"\nExpected file not found ({label}): {path}")
        sys.exit(1)
    print(f"Found {label}: {path}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="e.g. inputs/test.fasta")
    ap.add_argument("--out_dir", default="outputs", help="e.g. outputs")
    ap.add_argument("--dsressol_model", default="tools/DSResSol/models/dsressol_models.h5")
    args = ap.parse_args()

    fasta = args.fasta
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files_to_clean = [
        out_dir / "results_amp.csv",
        out_dir / "results_soluble.csv",
        out_dir / "results_nontoxic.csv",
        out_dir / "results_stable.csv",
        out_dir / "final_results.csv"
    ]

    for f in files_to_clean:
        if f.exists():
            f.unlink()
            print(f"🧹 Deleted old file: {f}")

    amp_csv = out_dir / "results_amp.csv"
    sol_csv = out_dir / "results_soluble.csv"
    nontox_csv = out_dir / "results_nontoxic.csv"
    stab_csv = out_dir / "results_stable.csv"
    final_csv = out_dir / "final_results.csv"

    # 1) AMP
    run([sys.executable, "src/run_pepnet_amp.py",
         "--fasta", fasta,
         "--out_dir", str(out_dir)])

    # 2) Solubility
    run([sys.executable, "src/run_dsressol_sol.py",
         "--fasta", fasta,
         "--model", args.dsressol_model,
         "--out_csv", str(sol_csv)])

    # 3) Non-toxic
    run([sys.executable, "src/run_toxipep_nontoxic.py",
         "--fasta", fasta,
         "--out_csv", str(nontox_csv)])

    # 4) Stability
    run([sys.executable, "src/run_stability.py",
         "--fasta", fasta,
         "--out_csv", str(stab_csv)])

    # sanity checks before merge
    must_exist(amp_csv, "AMP CSV")
    must_exist(sol_csv, "Soluble CSV")
    must_exist(nontox_csv, "NonToxic CSV")
    must_exist(stab_csv, "Stable CSV")

    # 5) Merge
    run([sys.executable, "src/merge_all_results.py",
         "--amp_csv", str(amp_csv),
         "--sol_csv", str(sol_csv),
         "--nontox_csv", str(nontox_csv),
         "--stab_csv", str(stab_csv),
         "--out_csv", str(final_csv)])

    
    intermediate_files = [
        amp_csv,
        sol_csv,
        nontox_csv,
        stab_csv
    ]
    for f in intermediate_files:
        if f.exists():
            f.unlink()
            print(f"Removed intermediate file: {f}")
        must_exist(final_csv, "Final merged CSV")
        print("\n Done! Final file:", final_csv)

if __name__ == "__main__":
    main()