## usage python 07_barplot.py --input-tsv 06_binders_final.tsv --output-html 07_barplot.html


import argparse
import pandas as pd
import plotly.express as px

HLA_SUPERTYPES = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01",
    "HLA-A*24:02", "HLA-A*26:01", "HLA-A*68:01",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01",
    "HLA-B*27:05", "HLA-B*39:01", "HLA-B*58:01"
]

def main(input_tsv, output_html):
    print(f"[INFO] Lecture du fichier {input_tsv} ...")
    df = pd.read_csv(input_tsv, sep="\t")

    # Filtrer uniquement sur les 12 supertypes
    df = df[df['HLA'].isin(HLA_SUPERTYPES)]

    # Compter le nombre de peptides par HLA et par Interpretation
    counts = df.groupby(['HLA', 'Interpretation']).size().reset_index(name='Count')
    
    # Prendre la palette par défaut
    default_colors = px.colors.qualitative.Plotly

    # Assigner manuellement les bonnes couleurs aux catégories
    color_discrete_map = {
        'Non binder': default_colors[0],     # bleu
        'Weak binder': default_colors[2],    # vert
        'Strong binder': default_colors[1]   # orange
    }


    fig = px.bar(
        counts,
        x='HLA',
        y='Count',
        color='Interpretation',
        color_discrete_map=color_discrete_map,
        barmode='group',
        title="Nombre de peptides par type d'interprétation et HLA supertype",
        labels={'Count': 'Nombre de peptides (log10)', 'HLA': 'HLA Supertype'}
    )

    fig.update_layout(
        xaxis_tickangle=-45,
        yaxis_type="log",
        xaxis=dict(tickfont=dict(size=14), title_font=dict(size=16)),
        yaxis=dict(tickfont=dict(size=14), title_font=dict(size=16)),
        title_font=dict(size=18)
    )

    # Sauvegarder le plot interactif en HTML
    fig.write_html(output_html)
    print(f"[INFO] Graphique sauvegardé dans {output_html}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Barplot des counts par HLA supertype et type d'interprétation")
    parser.add_argument("--input-tsv", required=True, help="Fichier TSV en entrée avec colonnes HLA et Interpretation")
    parser.add_argument("--output-html", required=True, help="Fichier HTML de sortie pour le graphique")
    args = parser.parse_args()

    main(args.input_tsv, args.output_html)
