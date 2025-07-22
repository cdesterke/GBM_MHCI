import pandas as pd

# Charger le fichier
df = pd.read_csv("df.csv", sep="\t", dtype=str)

# Nettoyer les colonnes pour éviter les erreurs
df.columns = [col.strip() for col in df.columns]

# Nettoyer les colonnes critiques
df = df[df['Reference_Allele'].notna() & df['Tumor_Seq_Allele2'].notna()]

# Assurer les types pour POS (1-based dans VCF)
df['POS'] = df['chromStart'].astype(int) + 1

# Assigner les valeurs obligatoires VCF
df['#CHROM'] = df['chrom'].astype(str)
df['ID'] = df['dbSNP_RS'].fillna('.')
df['REF'] = df['Reference_Allele']
df['ALT'] = df['Tumor_Seq_Allele2']
df['QUAL'] = '.'
df['FILTER'] = '.'


# Créer le champ INFO avec quelques métadonnées
df['INFO'] = ('GENE=' + df['Hugo_Symbol'].fillna('.') +
    ';TYPE=' + df['Variant_Classification'].fillna('.') +
    ';FREQ=' + df['freq'].fillna('.'))

# Garder uniquement les colonnes nécessaires au VCF
vcf_df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]

# Sauvegarder au format VCF
vcf_filename = "GBM.vcf"
with open(vcf_filename, 'w') as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene symbol\">\n")
    f.write("##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant classification\">\n")
    f.write("##INFO=<ID=FREQ,Number=1,Type=Float,Description=\"Allele frequency in tumor samples\">\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    vcf_df.to_csv(f, sep="\t", index=False, header=False)

print(f"✅ Fichier VCF généré avec succès : {vcf_filename}")

