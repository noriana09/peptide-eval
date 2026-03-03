# Peptide Multi-Objective Evaluation Pipeline

This project is a simple evaluation pipeline for generated peptide sequences.

It takes peptide sequences (FASTA file) as input and evaluates them across four properties:

- AMP (Antimicrobial Activity)
- Non-Toxicity
- Solubility
- Stability

At the end, it produces a single clean CSV file with all results merged together.

---

## What this project does

Given a FASTA file with peptide sequences, the pipeline:

1. Runs each external predictor separately
2. Collects probabilities and binary classifications (0/1)
3. Merges everything into one file
4. Adds summary columns:
   - `ALL` → 1 if all four properties are satisfied
   - `All above 0.95` → 1 if AMP, Toxicity and Solubility probabilities > 0.95
5. Adds a final `TOTAL` row showing how many sequences satisfy:
   - ALL
   - All above threshold

After execution, only one file remains:


---

## Required External Tools

This pipeline depends on the following external tools.  
They must be downloaded and placed inside a `tools/` directory.

### 1) AMP Prediction
PepNet  
http://liulab.top/PepNet/server

### 2) Toxicity Prediction
ToxiPep  
https://github.com/GGCL7/ToxiPep

### 3) Solubility Prediction
DSResSol  
https://github.com/mahan-fcb/DSResSol

Stability is computed using BioPython (ProteinAnalysis).

---

## Input

A FASTA file containing peptide sequences.



---

## How to Run

From the root directory:

```bash
python src/run_all.py --fasta inputs/test.fasta --out_dir outputs
