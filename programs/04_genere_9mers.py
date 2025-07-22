## python 04_genere_9mers.py --input cosmic_somatic.tsv  --cds_fasta MANE.GRCh38.v1.2.ensembl_protein.faa   --output peptides_9mer.tsv



import pandas as pd
from Bio import SeqIO
import argparse
import re

def parse_protein_fasta(fasta_path):
    tx2seq = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        desc = rec.description
        match = re.search(r"transcript:(ENST[0-9]+)\.\d+", desc)
        if match:
            tx_id = match.group(1)
            tx2seq[tx_id] = str(rec.seq)
    return tx2seq

def extract_pos_and_aa(hgvs_p):
    match = re.match(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})", hgvs_p)
    aa_3to1 = {
        'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G','His':'H','Ile':'I',
        'Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V',
        'Ter':'*'
    }
    if not match:
        return None, None
    pos = int(match.group(2))
    mut_aa = aa_3to1.get(match.group(3), None)
    return pos, mut_aa

def generate_9mers_sliding(seq, pos1, mut_aa):
    peptides = []
    for offset in range(9):
        start = pos1 - (offset + 1)
        end = start + 9
        if start < 0 or end > len(seq):
            continue
        wt = seq[start:end]
        mut = list(wt)
        mut[offset] = mut_aa
        peptides.append((offset + 1, wt, ''.join(mut)))
    return peptides

def main():
    parser = argparse.ArgumentParser(description="Génère 9 peptides 9-mer sliding window avec AA muté")
    parser.add_argument("--input", required=True, help="Fichier TSV avec mutations (cosmic_somatic.tsv)")
    parser.add_argument("--cds_fasta", required=True, help="FASTA protéique MANE avec annotation transcript")
    parser.add_argument("--output", required=True, help="Fichier TSV de sortie avec peptides")

    args = parser.parse_args()

    # Chargement des données mutationnelles
    df = pd.read_csv(args.input, sep="\t")

    # Renommer les colonnes si nécessaire pour correspondre aux attentes du script
    if "Feature_ID" in df.columns:
        df["Transcript_ID"] = df["Feature_ID"].str.split(".").str[0]
    if "HGVS.p" in df.columns:
        df["HGVS_p"] = df["HGVS.p"]

    protein_dict = parse_protein_fasta(args.cds_fasta)
    all_peptides = []

    for _, row in df.iterrows():
        tx_id = row["Transcript_ID"]
        hgvs_p = str(row["HGVS_p"])
        if tx_id not in protein_dict or not hgvs_p.startswith("p."):
            continue

        pos1, mut_aa = extract_pos_and_aa(hgvs_p)
        if pos1 is None or mut_aa is None:
            continue

        seq = protein_dict[tx_id]
        peptides = generate_9mers_sliding(seq, pos1, mut_aa)

        for mutant_position, wt_pep, mut_pep in peptides:
            record = row.to_dict()
            record["Mutant_AA_Position_in_9mer"] = mutant_position
            record["WT_9mer"] = wt_pep
            record["MUT_9mer"] = mut_pep
            all_peptides.append(record)

    out_df = pd.DataFrame(all_peptides)
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"{len(out_df)} peptides sliding 9-mer générés dans {args.output}")

if __name__ == "__main__":
    main()

