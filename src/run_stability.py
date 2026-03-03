import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="e.g. inputs/test.fasta")
    ap.add_argument("--out_csv", default="outputs/results_stable.csv")
    ap.add_argument("--threshold", type=float, default=40.0, help="Instability index <= threshold => Stability=1")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[1]
    fasta_path = (root / args.fasta).resolve()
    out_csv = (root / args.out_csv).resolve()

    rows = []
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        rid = str(rec.id)
        seq = str(rec.seq)

        try:
            instability = float(ProteinAnalysis(seq).instability_index())
        except Exception:
            instability = float("nan")

        stability_bin = 1 if (instability == instability and instability <= args.threshold) else 0  # nan-safe

        rows.append({
            "id": rid,
            "sequence": seq,
            "Stability_s": instability,
            "Stability": stability_bin,   # 1/0
        })

    df = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"Wrote: {out_csv}")
    print(df.head())

if __name__ == "__main__":
    main()