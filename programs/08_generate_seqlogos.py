import pandas as pd
import logomaker
import matplotlib.pyplot as plt
import io
import base64
from collections import defaultdict
import os

# Couleurs classiques par acide aminé
AA_COLORS = {
    'A': '#8dd3c7', 'C': '#ffffb3', 'D': '#bebada', 'E': '#fb8072',
    'F': '#80b1d3', 'G': '#fdb462', 'H': '#b3de69', 'I': '#fccde5',
    'K': '#d9d9d9', 'L': '#bc80bd', 'M': '#ccebc5', 'N': '#ffed6f',
    'P': '#d9d9d9', 'Q': '#d9d9d9', 'R': '#d9d9d9', 'S': '#d9d9d9',
    'T': '#d9d9d9', 'V': '#d9d9d9', 'W': '#d9d9d9', 'Y': '#d9d9d9'
}

def build_pwm(peptides):
    if len(peptides) == 0:
        return None
    length = len(peptides[0])
    counts = [defaultdict(int) for _ in range(length)]
    for pep in peptides:
        for i, aa in enumerate(pep):
            counts[i][aa] += 1
    pwm_data = {}
    all_aas = sorted({aa for pos in counts for aa in pos.keys()})
    for aa in all_aas:
        pwm_data[aa] = [counts[i].get(aa, 0) for i in range(length)]
    pwm_df = pd.DataFrame(pwm_data)
    pwm_df = pwm_df / pwm_df.sum(axis=1).values[:, None]
    return pwm_df

def generate_logo(pwm_df):
    plt.figure(figsize=(max(6, pwm_df.shape[0]*0.6), 3))
    ax = plt.gca()
    pwm_df.index = range(1, pwm_df.shape[0] + 1) 
    logo = logomaker.Logo(pwm_df, ax=ax, color_scheme=AA_COLORS)
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Position')
    ax.set_ylim(0, 1)
    plt.tight_layout()

    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    plt.close()
    buf.seek(0)
    img_bytes = buf.read()
    img_b64 = base64.b64encode(img_bytes).decode('utf-8')
    return img_b64

def create_output_dir(path="output"):
    if not os.path.exists(path):
        os.makedirs(path)

def create_html(hla_logos, out_html_path):
    html_header = """
    <html><head><title>Seqlogos HLA</title></head><body>
    <h1>Seqlogos par HLA</h1>
    """
    html_footer = "</body></html>"
    content = ""
    for hla, (img_b64, pep_count) in hla_logos.items():
        content += f"<h2>{hla} (n={pep_count} peptides)</h2>\n"
        content += f'<img src="data:image/png;base64,{img_b64}" alt="Logo {hla}"><br><br>\n'

    with open(out_html_path, "w") as f:
        f.write(html_header + content + html_footer)

def main():
    df = pd.read_csv("06_binders_final.tsv", sep='\t')

    df_binders = df[~df['Interpretation'].str.contains("Non-binder", case=False, na=False)]
    print(f"Nombre total de peptides binders : {len(df_binders)}")

    hla_groups = df_binders.groupby('HLA')['Peptide'].apply(list)

    hla_logos = {}
    for hla, peptides in hla_groups.items():
        if len(peptides) < 5:
            print(f"[INFO] Trop peu de peptides pour {hla} ({len(peptides)}). Ignoré.")
            continue
        lengths = set(len(p) for p in peptides)
        if len(lengths) != 1:
            print(f"[WARNING] Peptides de longueurs différentes pour {hla}, ignoré.")
            continue
        pwm_df = build_pwm(peptides)
        if pwm_df is None:
            continue
        logo_img_b64 = generate_logo(pwm_df)
        hla_logos[hla] = (logo_img_b64, len(peptides))  # <-- on stocke aussi le n

    if len(hla_logos) == 0:
        print("[WARNING] Aucun logo généré.")
        return

    create_output_dir()
    out_html = "08_seqlogos.html"
    create_html(hla_logos, out_html)
    print(f"[INFO] Logos générés et enregistrés dans {out_html}")

if __name__ == "__main__":
    main()

