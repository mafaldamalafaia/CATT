import pandas as pd
import re

df = pd.read_csv("fh_LDLR_clinvar_full.csv")


def map_label(x):
    if pd.isna(x):
        return None
    x = x.lower()
    if "pathogenic" in x:
        return 1
    if "benign" in x:
        return 0
    return None


def extract_protein(name):
    m = re.search(r"\(p\.([^)]+)\)", str(name))
    return m.group(1) if m else None


def extract_cdna(name):
    m = re.search(r":c\.([^ ]+)", str(name))
    return m.group(1) if m else None


df["protein_variant"] = df["Name"].apply(extract_protein)

df["y"] = df["ClinicalSignificance"].apply(map_label)
df = df[df["y"].notna()].copy()

print(df["y"].value_counts())

df["protein_variant"] = df["Name"].apply(extract_protein)
df["cdna_variant"] = df["Name"].apply(extract_cdna)
df["genomic_pos"] = df["PositionVCF"]
df["ref"] = df["ReferenceAlleleVCF"]
df["alt"] = df["AlternateAlleleVCF"]
df["variant_type"] = df["Type"]

final_cols = [
    "protein_variant",
    "cdna_variant",
    "genomic_pos",
    "ref",
    "alt",
    "variant_type",
    "y"
]

bio_df = df[final_cols].dropna()

bio_df.to_csv("LDLR_biology.csv", index=False)
