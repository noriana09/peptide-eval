import argparse
import subprocess
from pathlib import Path
import pandas as pd

def parse_fasta_ids_and_seqs(fasta_path: Path):
    ids, seqs = [], []
    cur_id, cur_seq = None, []
    for line in fasta_path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_id is not None:
                ids.append(cur_id)
                seqs.append("".join(cur_seq))
            cur_id = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line)
    if cur_id is not None:
        ids.append(cur_id)
        seqs.append("".join(cur_seq))
    return ids, seqs

def run(cmd, cwd=None, env=None):
    print("\n▶ " + " ".join(cmd))
    p = subprocess.run(cmd, cwd=cwd, env=env, text=True)
    if p.returncode != 0:
        raise SystemExit(p.returncode)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="e.g. inputs/test.fasta")
    ap.add_argument("--out_csv", default="outputs/results_nontoxic.csv")
    ap.add_argument("--toxipep_txt", default="outputs/toxipep_raw.txt")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[1]
    fasta_path = (root / args.fasta).resolve()
    out_csv = (root / args.out_csv).resolve()
    tox_txt = (root / args.toxipep_txt).resolve()

    # Read FASTA ids+seqs
    ids, seqs = parse_fasta_ids_and_seqs(fasta_path)

    # Call ToxiPep predict.py using conda env python directly
    toxipep_py = Path.home() / "miniconda3" / "envs" / "toxipep" / "bin" / "python"
    predict_py = root / "tools" / "ToxiPep" / "Code" / "predict.py"

    if not toxipep_py.exists():
        raise FileNotFoundError(f"toxipep python not found: {toxipep_py}")
    if not predict_py.exists():
        raise FileNotFoundError(f"predict.py not found: {predict_py}")

    env = dict(**{k: v for k, v in dict(**__import__("os").environ).items()})
    env["KMP_DUPLICATE_LIB_OK"] = "TRUE"

    # Ensure we run from Code/ so it finds best_model_0.9.pth
    code_dir = (root / "tools" / "ToxiPep" / "Code").resolve()

    run([str(toxipep_py), "predict.py", "-i", str(fasta_path), "-o", str(tox_txt)], cwd=str(code_dir), env=env)

    # Parse raw txt: "sequence, label, prob"
    rows = []
    lines = tox_txt.read_text().splitlines()
    if len(lines) != len(seqs):
        print(f"⚠️ Warning: raw lines ({len(lines)}) != fasta seqs ({len(seqs)})")

    for i, line in enumerate(lines):
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 3:
            raise ValueError(f"Bad line: {line}")
        seq_out = parts[0]
        tox_label = int(parts[1])
        tox_prob = float(parts[2])

        # Non-toxic is the complement of toxic
        non_tox = 1 - tox_label
        non_tox_prob = 1.0 - tox_prob

        # Prefer FASTA ids and sequences for stability
        rid = ids[i] if i < len(ids) else f"seq{i+1}"
        rseq = seqs[i] if i < len(seqs) else seq_out

        rows.append({
            "id": rid,
            "sequence": rseq,
            "TOX_prob": tox_prob,
            "TOX": tox_label,
            "NON_TOX_prob": non_tox_prob,
            "NON_TOX": non_tox,
        })

    df = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"Wrote: {out_csv}")
    print(df.head())

if __name__ == "__main__":
    main()