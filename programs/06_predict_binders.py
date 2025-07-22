## usage python 06_predict_binders.py --fasta peptides_mut.fasta

import argparse
import pandas as pd
import subprocess
from Bio import SeqIO
import uuid
import plotly.express as px
import os
import glob

HLA_SUPERTYPES = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01",
    "HLA-A*24:02", "HLA-A*26:01", "HLA-A*68:01",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*27:05",
    "HLA-B*39:01", "HLA-B*58:01", "HLA-B*15:01"
]

def read_peptides_from_fasta(fasta_file):
    peptides = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).strip().upper()
        if len(seq) == 9:
            peptides.append({'id': record.id, 'sequence': seq})
    return peptides

def write_input_csv(peptides, csv_file):
    rows = []
    for p in peptides:
        for allele in HLA_SUPERTYPES:
            rows.append({
                "peptide": p["sequence"],
                "allele": allele,
                "seq_id": p["id"]
            })
    df = pd.DataFrame(rows)
    df.to_csv(csv_file, index=False)
    return df, csv_file

def run_mhcflurry(input_csv, output_csv):
    subprocess.run([
        "mhcflurry-predict",
        input_csv,
        "--out", output_csv
    ], check=True)
    return output_csv

def classify_affinity(affinity):
    if affinity < 50:
        return "Strong binder"
    elif affinity < 500:
        return "Weak binder"
    else:
        return "Non-binder"

def generate_html_plot(df, html_file="binders_plot.html"):
    # Récupération de la palette qualitative Plotly par défaut
    default_colors = px.colors.qualitative.Plotly

    # Inversion des couleurs entre Strong et Weak binder
    color_map = {
        "Strong binder": default_colors[1],  # couleur initialement pour Weak binder
        "Weak binder": default_colors[2],    # couleur initialement pour Strong binder
        "Non-binder": default_colors[0]
    }

    fig = px.scatter(
        df,
        x="Peptide",
        y="Affinity_nM",
        color="Interpretation",
        color_discrete_map=color_map,
        hover_data=["Sequence_ID", "HLA"],
        title="Affinités MHC-I des peptides par allèle",
        labels={"Affinity_nM": "Affinité (nM)", "Peptide": "Peptide"}
    )
    fig.update_yaxes(type="log", autorange=True)

    # Masquer les étiquettes de l'axe des x
    fig.update_xaxes(showticklabels=False)

    fig.update_traces(
        hoverlabel=dict(
            bgcolor=None,
            bordercolor=None,
            font_size=12,
            font_family="Arial"
        )
    )
    fig.write_html(html_file)
    print(f"[INFO] Graphique interactif enregistré dans {html_file}")


def cleanup_temp_files():
    patterns = ["mhcflurry_input_*.csv", "mhcflurry_output_*.csv"]
    for pattern in patterns:
        for f in glob.glob(pattern):
            try:
                os.remove(f)
                print(f"[INFO] Supprimé le fichier temporaire : {f}")
            except Exception as e:
                print(f"[WARNING] Impossible de supprimer {f} : {e}")

def main(fasta_path):
    print("[INFO] Lecture des peptides 9-mers dans le fichier FASTA...")
    peptides = read_peptides_from_fasta(fasta_path)

    print(f"[INFO] {len(peptides)} peptides lus. Préparation du CSV d'entrée...")
    input_csv = f"mhcflurry_input_{uuid.uuid4().hex}.csv"
    output_csv = f"mhcflurry_output_{uuid.uuid4().hex}.csv"

    input_df, _ = write_input_csv(peptides, input_csv)

    print("[INFO] Exécution de mhcflurry-predict...")
    run_mhcflurry(input_csv, output_csv)

    print("[INFO] Lecture des résultats de mhcflurry...")
    prediction_df = pd.read_csv(output_csv)

    print("[INFO] Fusion des données et ajout des colonnes d'interprétation...")
    full_df = input_df.copy()
    full_df["Affinity_nM"] = prediction_df["mhcflurry_affinity"]

    full_df = full_df.rename(columns={
        "peptide": "Peptide",
        "allele": "HLA",
        "seq_id": "Sequence_ID"
    })

    full_df["Interpretation"] = full_df["Affinity_nM"].apply(classify_affinity)

    print("[INFO] Sauvegarde du fichier complet avec tous les résultats...")
    full_df.to_csv("06_binders_final.tsv", sep="\t", index=False)

    print("[INFO] Extraction des meilleurs binders par peptide...")
    best_binders_df = full_df.loc[full_df.groupby("Sequence_ID")["Affinity_nM"].idxmin()]
    best_binders_df.to_csv("06_best_binders_by_peptide.tsv", sep="\t", index=False)

    print("[INFO] Génération du graphique interactif HTML...")
    generate_html_plot(full_df, "06_binders_plot.html")

    print("[INFO] Nettoyage des fichiers temporaires...")
    cleanup_temp_files()

    print("[✔] Analyse terminée.")
    print("Fichiers générés dans : 06_binders_final.tsv, 06_best_binders_by_peptide.tsv, 06_binders_plot.html")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prédiction binders MHC-I avec mhcflurry et sortie HTML interactive")
    parser.add_argument("--fasta", required=True, help="Fichier FASTA contenant peptides 9-mers")
    args = parser.parse_args()
    main(args.fasta)
