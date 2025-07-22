import sys
import csv

def parse_info_field(info_str):
    """Parse le champ INFO en dict clé=valeur."""
    info_dict = {}
    if info_str == ".":
        return info_dict
    for entry in info_str.split(";"):
        if "=" in entry:
            key, val = entry.split("=", 1)
            info_dict[key] = val
        else:
            info_dict[entry] = True  # flag sans valeur
    return info_dict

def parse_ann_field(ann_str):
    """
    Parse le champ ANN (SnpEff), multiple annotations séparées par ','
    Chaque annotation est un pipe-separated fields.
    """
    ann_list = []
    if not ann_str:
        return ann_list
    for ann in ann_str.split(","):
        fields = ann.split("|")
        ann_list.append(fields)
    return ann_list

def main(vcf_path, tsv_path):
    with open(vcf_path) as vcf, open(tsv_path, "w", newline="") as tsvfile:
        writer = None

        for line in vcf:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").strip().split("\t")
                # Colonnes classiques du VCF
                base_cols = header[:8]  # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
                
                # On va élargir avec les clés INFO qu'on rencontre au fur et à mesure,
                # mais comme on ne sait pas encore lesquelles, on attend la première ligne de données
                continue

            # ligne variant
            fields = line.strip().split("\t")
            base_values = fields[:8]
            info_field = base_values[7]

            info_dict = parse_info_field(info_field)

            # On récupère les annotations ANN
            ann_list = parse_ann_field(info_dict.get("ANN", ""))

            # Au premier variant, on construit le header complet
            if writer is None:
                # Colonnes INFO à ajouter (hors ANN)
                info_keys = sorted(k for k in info_dict.keys() if k != "ANN")

                # Colonnes ANN (fixées ici d'après doc SnpEff)
                ann_columns = [
                    "Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID",
                    "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
                    "HGVS.c", "HGVS.p", "cDNA.pos/cDNA.length", "CDS.pos/CDS.length",
                    "AA.pos/AA.length", "Distance", "ERRORS/WARNINGS/INFO"
                ]

                # Header final
                full_header = base_cols[:-1] + info_keys + ann_columns
                writer = csv.DictWriter(tsvfile, fieldnames=full_header, delimiter="\t", extrasaction="ignore")
                writer.writeheader()

            # Construction du dict ligne
            row = dict(zip(base_cols[:-1], base_values[:-1]))

            # Ajouter les infos hors ANN
            for k in info_keys:
                row[k] = info_dict.get(k, "")

            # Prendre la première annotation ANN seulement (si plusieurs, tu peux modifier pour toutes)
            if ann_list:
                ann_first = ann_list[0]
                # On complète avec la taille de ann_columns, certains champs peuvent manquer
                for i, col in enumerate(ann_columns):
                    row[col] = ann_first[i] if i < len(ann_first) else ""
            else:
                for col in ann_columns:
                    row[col] = ""

            writer.writerow(row)

    print(f"Extraction terminée, résultat dans {tsv_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python vcf_to_tsv.py input.vcf output.tsv")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

