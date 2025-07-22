import pandas as pd

# Lecture du fichier Cosmic
cosmic = pd.read_csv("Cosmic_MutantCensus_v102_GRCh38.tsv", sep="\t", low_memory=False)

# Sélection des colonnes
cosmic = cosmic[[
    "CHROMOSOME", "GENOME_START", "GENOMIC_WT_ALLELE", "GENOMIC_MUT_ALLELE",
    "MUTATION_DESCRIPTION", "MUTATION_SOMATIC_STATUS", "LEGACY_MUTATION_ID",
    "MUTATION_AA", "MUTATION_CDS"
]]

# Filtrage des mutations d'intérêt
cosmic = cosmic[
    (cosmic["MUTATION_DESCRIPTION"] == "missense_variant") &
    (cosmic["MUTATION_SOMATIC_STATUS"] == "Confirmed somatic variant")
]

# Lecture du fichier d'annotation
data = pd.read_csv("gbm.ann.tsv", sep="\t")

# Filtrer sur les missense variants
data = data[data["Annotation"] == "missense_variant"]

# Supprimer les colonnes inutiles si elles existent
cols_to_drop = ["Distance", "ERRORS/WARNINGS/INFO"]
data = data.drop(columns=[col for col in cols_to_drop if col in data.columns])

# Fusionner avec Cosmic sur les clés génomiques
merged = pd.merge(
    data, cosmic,
    how="inner",
    left_on=["CHROM", "POS", "REF", "ALT"],
    right_on=["CHROMOSOME", "GENOME_START", "GENOMIC_WT_ALLELE", "GENOMIC_MUT_ALLELE"]
)

# Éliminer les doublons
merged = merged.drop_duplicates()

# Exporter le résultat
merged.to_csv("cosmic_somatic.tsv", sep="\t", index=False)

