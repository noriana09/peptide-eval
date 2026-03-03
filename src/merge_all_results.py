import argparse
from pathlib import Path
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amp_csv", default="outputs/results_amp.csv")
    ap.add_argument("--sol_csv", default="outputs/results_soluble.csv")
    ap.add_argument("--nontox_csv", default="outputs/results_nontoxic.csv")
    ap.add_argument("--stab_csv", default="outputs/results_stable.csv")
    ap.add_argument("--out_csv", default="outputs/final_results.csv")
    ap.add_argument(
        "--prob_thr",
        type=float,
        default=0.95,
        help="threshold for 'All above 0.95' (excluding stability)"
    )
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[1]

    amp_path = (root / args.amp_csv).resolve()
    sol_path = (root / args.sol_csv).resolve()
    nontox_path = (root / args.nontox_csv).resolve()
    stab_path = (root / args.stab_csv).resolve()
    out_path = (root / args.out_csv).resolve()

    amp = pd.read_csv(amp_path)
    sol = pd.read_csv(sol_path)
    nontox = pd.read_csv(nontox_path)
    stab = pd.read_csv(stab_path)

    # -------------------------
    # Normalize columns (AMP)
    # -------------------------
    if "AMP_prob" not in amp.columns and "probability" in amp.columns:
        amp = amp.rename(columns={"probability": "AMP_prob"})
    if "AMP" not in amp.columns and "Binary" in amp.columns:
        amp = amp.rename(columns={"Binary": "AMP"})
    if "sequence" not in amp.columns and "Sequence" in amp.columns:
        amp = amp.rename(columns={"Sequence": "sequence"})

    # -------------------------
    # Normalize columns (SOL)
    # -------------------------
    if "SOL_prob" not in sol.columns and "probability" in sol.columns:
        sol = sol.rename(columns={"probability": "SOL_prob"})
    if "SOL" not in sol.columns and "Binary" in sol.columns:
        sol = sol.rename(columns={"Binary": "SOL"})
    if "sequence" not in sol.columns and "Sequence" in sol.columns:
        sol = sol.rename(columns={"Sequence": "sequence"})

    # -------------------------
    # Normalize columns (NON-TOX)
    # -------------------------
    if "NON_TOX_prob" not in nontox.columns:
        raise ValueError("nontox csv must contain NON_TOX_prob")
    if "NON_TOX" not in nontox.columns:
        raise ValueError("nontox csv must contain NON_TOX")

    # -------------------------
    # Normalize columns (STABILITY)
    # expecting: id, sequence, Stability_s, Stability
    # -------------------------
    if "Stability_s" not in stab.columns and "STAB_index" in stab.columns:
        stab = stab.rename(columns={"STAB_index": "Stability_s"})
    if "Stability" not in stab.columns and "STABLE" in stab.columns:
        stab = stab.rename(columns={"STABLE": "Stability"})

    if "Stability_s" not in stab.columns or "Stability" not in stab.columns:
        raise ValueError(
            f"stability csv must contain Stability_s and Stability. Found: {list(stab.columns)}"
        )

    # -------------------------
    # Merge
    # -------------------------
    df = amp.merge(sol, on=["id", "sequence"], how="inner")
    df = df.merge(nontox, on=["id", "sequence"], how="inner")
    df = df.merge(
        stab[["id", "sequence", "Stability_s", "Stability"]],
        on=["id", "sequence"],
        how="inner"
    )

    # Rename to your final naming convention
    df = df.rename(columns={
        "NON_TOX_prob": "Tox_prob",
        "NON_TOX": "Toxicity",
        "SOL_prob": "Sol_prob",
        "SOL": "Solubility",
    })

    # -------------------------
    # Derived columns
    # -------------------------
    # ALL = 1 when all binaries are 1
    df["ALL"] = (
        (df["AMP"] == 1) &
        (df["Toxicity"] == 1) &
        (df["Solubility"] == 1) &
        (df["Stability"] == 1)
    ).astype(int)

    # All above threshold (excluding stability)
    df["All above 0.95"] = (
        (df["AMP_prob"] > args.prob_thr) &
        (df["Tox_prob"] > args.prob_thr) &
        (df["Sol_prob"] > args.prob_thr)
    ).astype(int)

    # Final column order
    keep = [
        "id", "sequence",
        "AMP_prob", "Tox_prob", "Sol_prob", "Stability_s",
        "AMP", "Toxicity", "Solubility", "Stability",
        "ALL", "All above 0.95"
    ]
    df = df[keep]

    # ---------------------------
    # Add summary row at the end
    # ---------------------------
    total_all = int((df["ALL"] == 1).sum())
    total_above = int((df["All above 0.95"] == 1).sum())

    summary_row = {
        "id": "TOTAL",
        "sequence": "",
        "AMP_prob": "",
        "Tox_prob": "",
        "Sol_prob": "",
        "Stability_s": "",
        "AMP": "",
        "Toxicity": "",
        "Solubility": "",
        "Stability": "",
        "ALL": total_all,
        "All above 0.95": total_above
    }

    df = pd.concat([df, pd.DataFrame([summary_row])], ignore_index=True)

    # Write
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)

    print(f"Wrote merged CSV: {out_path}")
    print(df.head())

    print("\nCounts:")
    print("AMP:", int((df["AMP"] == 1).sum()))
    print("Non-Toxic:", int((df["Toxicity"] == 1).sum()))
    print("Soluble:", int((df["Solubility"] == 1).sum()))
    print("Stable:", int((df["Stability"] == 1).sum()))
    print("ALL:", total_all)
    print(f"All above {args.prob_thr:.2f}:", total_above)

if __name__ == "__main__":
    main()