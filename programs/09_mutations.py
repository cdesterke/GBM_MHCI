import pandas as pd
import plotly.express as px

# 1. Chargement des données depuis un fichier CSV ou TSV
# Remplacez par le nom réel du fichier si nécessaire
df = pd.read_csv("cosmic_somatic.tsv", sep="\t")

# 2. Nettoyage des types
df['FREQ'] = pd.to_numeric(df['FREQ'], errors='coerce')
df['POS'] = pd.to_numeric(df['POS'], errors='coerce')

# 3. Création du scatter plot interactif
fig = px.scatter(
    df,
    x="POS",
    y="GENE",
    color="GENE",
    size="FREQ",
    hover_data=[
        "CHROM", "POS", "GENE", "FREQ", "Annotation", "HGVS.c", "HGVS.p", "COSMIC_ID" if "COSMIC_ID" in df.columns else "MUTATION_DESCRIPTION"
    ],
    title="Visualisation des mutations somatiques",
    labels={"POS": "Position sur le chromosome", "GENE": "Gène", "FREQ": "Fréquence"},
    template="plotly_white"
)

# 4. Exporter en HTML
fig.write_html("09_mutations_plot.html")

# Afficher dans un notebook si souhaité
fig.show()

