# usage : python 05_fasta.py --input peptides_9mer.tsv   --wt_fasta peptides_wt.fasta --mut_fasta peptides_mut.fasta

import pandas as pd
import argparse

def sanitize(text):
    """Nettoie les chaînes pour en faire des identifiants valides pour FASTA"""
    return str(text).replace(" ", "_").replace(".", "_").replace("/", "_")

def main():
    parser = argparse.ArgumentParser(description="Génère deux fichiers FASTA pour peptides WT et mutés à partir d’un TSV")
    parser.add_argument("--input", required=True, help="Fichier TSV contenant WT_9mer et MUT_9mer")
    parser.add_argument("--wt_fasta", default="wt.fasta", help="Fichier de sortie FASTA pour les peptides normaux")
    parser.add_argument("--mut_fasta", default="mut.fasta", help="Fichier de sortie FASTA pour les peptides mutés")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    # Séparer correctement la colonne combinée "Transcript_IDHGVS_p" si nécessaire
    if "Transcript_IDHGVS_p" in df.columns and "Transcript_ID" not in df.columns:
        df[["Transcript_ID", "HGVS_p"]] = df["Transcript_IDHGVS_p"].str.extract(r"(ENST\d+)\s*(p\.\w+)")

    wt_records = []
    mut_records = []

    for _, row in df.iterrows():
        gene = sanitize(row.get("Gene_Name", "NA"))
        tx = sanitize(row.get("Transcript_ID", "NA"))
        pos = row.get("Mutant_AA_Position_in_9mer", "NA")
        chrom = row.get("CHROM", "NA")
        position = row.get("POS", "NA")
        ref = row.get("REF", "X")
        alt = row.get("ALT", "X")

        base_id = f"{gene}_{tx}_pos{pos}_chr{chrom}_{position}_{ref}>{alt}"

        wt_seq = row.get("WT_9mer", "")
        mut_seq = row.get("MUT_9mer", "")

        if len(wt_seq) == 9 and len(mut_seq) == 9:
            wt_records.append(f">{base_id}\n{wt_seq}")
            mut_records.append(f">{base_id}\n{mut_seq}")

    with open(args.wt_fasta, "w") as f:
        f.write("\n".join(wt_records) + "\n")

    with open(args.mut_fasta, "w") as f:
        f.write("\n".join(mut_records) + "\n")

    print(f"✅ FASTA WT : {args.wt_fasta}")
    print(f"✅ FASTA muté : {args.mut_fasta}")
    print(f"🔢 {len(wt_records)} peptides écrits dans chaque fichier")

if __name__ == "__main__":
    main()
