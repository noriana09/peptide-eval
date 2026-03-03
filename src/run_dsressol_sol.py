import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras import models

CODES = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

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

def create_char_dict(codes):
    return {aa: i + 1 for i, aa in enumerate(codes)}  # 1..20, unknown -> 0

def integer_encoding(seqs, char_dict):
    enc = []
    for s in seqs:
        enc.append(np.array([char_dict.get(ch, 0) for ch in s], dtype=np.int32))
    return enc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--model", required=True, help="Path to DSResSol .h5/.keras model")
    ap.add_argument("--out_csv", default="outputs/results_soluble.csv")
    ap.add_argument("--threshold", type=float, default=0.5)
    args = ap.parse_args()

    project_root = Path(__file__).resolve().parents[1]
    fasta_path = (project_root / args.fasta).resolve()
    model_path = (project_root / args.model).resolve() if not Path(args.model).is_absolute() else Path(args.model)
    out_csv = (project_root / args.out_csv).resolve()

    ids, seqs = parse_fasta_ids_and_seqs(fasta_path)

    model = models.load_model(str(model_path))

    # infer max_length from model input shape
    in_shape = model.input_shape[0] if isinstance(model.input_shape, list) else model.input_shape
    max_length = int(in_shape[1])

    char_dict = create_char_dict(CODES)
    enc = integer_encoding(seqs, char_dict)
    x = pad_sequences(enc, maxlen=max_length, padding="post", truncating="post")

    y_pred = model.predict([x], verbose=0)

    print("DEBUG y_pred shape:", getattr(y_pred, "shape", None))

    # Handle common output shapes:
    # (N, 1) -> probability column 0
    # (N, 2) -> take class-1 probability (usually "soluble")
    # (N,)   -> already ok
    if y_pred.ndim == 2 and y_pred.shape[1] == 1:
        sol_prob = y_pred[:, 0].astype(float)
    elif y_pred.ndim == 2 and y_pred.shape[1] == 2:
        sol_prob = y_pred[:, 1].astype(float)  # class 1 prob
    elif y_pred.ndim == 1:
        sol_prob = y_pred.astype(float)
    else:
        raise ValueError(f"Unexpected y_pred shape: {y_pred.shape}")

    df = pd.DataFrame({
        "id": ids,
        "sequence": seqs,
        "SOL_prob": sol_prob,
        "SOL": (sol_prob > args.threshold).astype(int),
    })

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"Wrote: {out_csv}")
    print(df.head())

if __name__ == "__main__":
    main()