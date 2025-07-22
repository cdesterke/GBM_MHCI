## usage python 10_scatter2.py --peptides peptides_9mer.tsv --binders 06_binders_final.tsv
import pandas as pd
import numpy as np
import plotly.express as px
import argparse

def main():
    parser = argparse.ArgumentParser(description="Analyse peptides et binders")
    parser.add_argument('--peptides', type=str, required=True, help="Fichier peptides TSV")
    parser.add_argument('--binders', type=str, required=True, help="Fichier binders TSV")
    args = parser.parse_args()

    # 1. Lecture des fichiers
    pep = pd.read_csv(args.peptides, sep="\t")
    binder = pd.read_csv(args.binders, sep="\t")

    # 2. Créer colonne 'conca'
    pep['conca'] = pep['Gene_Name'] + "_" + pep['HGVS.p']

    # 3. Garder colonnes utiles et supprimer doublons
    pep_unique = pep[['conca', 'FREQ', 'MUT_9mer', 'Mutant_AA_Position_in_9mer','LEGACY_MUTATION_ID']].drop_duplicates()

    # 4. Jointure many-to-many
    merged = pd.merge(
        pep_unique,
        binder,
        left_on='MUT_9mer',
        right_on='Peptide',
        how='left'
    )

    # 5. Filtrer les lignes sans "Non" dans Interpretation
    filtered = merged[~merged['Interpretation'].str.contains("Non", na=False)].copy()

    # --- Export du dataset filtré ---
    filtered.to_csv("10_peptides_mutations.tsv", sep="\t", index=False)
    print("✅ Dataset filtré exporté dans '10_peptides_mutations.tsv'")

    # 8. Graphique
    fig = px.scatter(
        filtered,
        x='FREQ',
        y='Affinity_nM',
        color='HLA',
        symbol='Interpretation',
        hover_data=['conca', 'FREQ', 'MUT_9mer'],
        title="Affinité en fonction de Frequence, couleur (HLA) et forme (Interpretation)",
        log_y=True,
        size_max=20
    )

    fig.update_layout(
    xaxis_title="Frequency",
    yaxis_title="Affinity (nM)",
    showlegend=False
    )



    # 9. Ligne horizontale à y=50
    fig.add_hline(
        y=50,
        line_dash="dash",
        line_color="black"
    )

    # 10. Export HTML avec autoresize
    fig.write_html("10_peptides_selection.html", full_html=True, include_plotlyjs='cdn', config={"responsive": True})

    print("✅ Graphique généré pour l'étape 10 avec redimensionnement dynamique")

if __name__ == "__main__":
    main()

