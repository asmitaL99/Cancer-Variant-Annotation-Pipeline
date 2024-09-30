#Cancer Variant Annotation Pipeline in Python

import pandas as pd
from Bio import SeqIO
from biomart import BiomartServer
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------
# Step 1: Read Variant Data
# ---------------------------
# Define input VCF file path 
vcf_file = "data/variants.vcf"

# Read the VCF file using Biopython's SeqIO
vcf = SeqIO.parse(vcf_file, "vcf")

# Convert the VCF records to a DataFrame
vcf_records = []
for record in vcf:
    for variant in record.features:
        vcf_records.append({
            'seqname': record.id,
            'start': variant.location.start,
            'end': variant.location.end,
            'REF': record.ref,
            'ALT': record.alts[0] if record.alts else None
        })
vcf_df = pd.DataFrame(vcf_records)

# Print a summary
print("Step 1: VCF file successfully read.")
print(vcf_df.head())

# ---------------------------
# Step 2: Annotate Variants
# ---------------------------
# Connect to Ensembl BioMart
server = BiomartServer("http://www.ensembl.org/biomart")
hsapiens = server.datasets['hsapiens_gene_ensembl']

# Prepare variant IDs
vcf_df['variant_id'] = vcf_df['seqname'] + "_" + vcf_df['start'].astype(str) + "_" + vcf_df['REF'] + "/" + vcf_df['ALT']

# Query Ensembl VEP for annotations
results = hsapiens.search({
    'attributes': [
        'ensembl_gene_id', 'external_gene_name', 'consequence_type_tv',
        'protein_start', 'protein_end', 'sift_score', 'polyphen_score'
    ],
    'filters': {'variant_id': vcf_df['variant_id'].tolist()}
})

# Convert the results into a DataFrame
vep_df = pd.DataFrame(results)

# Merge annotations with original variant data
annotated_variants = pd.merge(vcf_df, vep_df, left_on='variant_id', right_on='Uploaded_variation')

# Save the annotated variants to a CSV file
annotated_variants.to_csv("results/annotated_variants.csv", index=False)

print("Step 2: Variant annotation completed.")

# ---------------------------
# Step 3: Predict Functional Impact
# ---------------------------
# Convert SIFT and PolyPhen scores to numeric
annotated_variants['sift_score'] = pd.to_numeric(annotated_variants['sift_score'], errors='coerce')
annotated_variants['polyphen_score'] = pd.to_numeric(annotated_variants['polyphen_score'], errors='coerce')

# Save the results
annotated_variants.to_csv("results/functional_impact_variants.csv", index=False)

print("Step 3: Functional impact prediction completed.")

# ---------------------------
# Step 4: Pathway Enrichment Analysis
# ---------------------------
# Extract unique gene symbols for pathway analysis
genes = annotated_variants['external_gene_name'].unique()

# Map gene symbols to Entrez IDs (using gseapy's `id_conversion` function)
gene_mapping = gp.id_mapping(genes, 'symbol', 'entrezgene', 'hsapiens')

# Perform pathway enrichment analysis using Reactome
enr = gp.enrichr(gene_list=gene_mapping['entrezgene'].tolist(),
                 gene_sets='Reactome_2016', 
                 organism='Human', 
                 cutoff=0.05)

# Save the pathway enrichment results
enr.results.to_csv("results/pathway_enrichment_results.csv", index=False)

print("Step 4: Pathway enrichment analysis completed.")

# ---------------------------
# Step 5: Visualization
# ---------------------------
# Plot the distribution of variant consequences
plt.figure(figsize=(10, 6))
sns.countplot(data=annotated_variants, x='consequence_type_tv', color='steelblue')
plt.xticks(rotation=45)
plt.xlabel('Variant Consequence')
plt.ylabel('Count')
plt.title('Distribution of Variant Consequences')
plt.tight_layout()

# Save the consequence plot
plt.savefig("results/variant_consequences.png")
plt.close()

# Visualize the top 10 enriched pathways
enr_results = enr.results.sort_values('Adjusted P-value').head(10)
plt.figure(figsize=(8, 6))
sns.barplot(data=enr_results, y='Term', x='Adjusted P-value', color='blue')
plt.xlabel('Adjusted P-value')
plt.ylabel('Pathway')
plt.title('Top 10 Enriched Pathways')
plt.tight_layout()

# Save the pathway enrichment plot
plt.savefig("results/pathway_enrichment.png")
plt.close()

print("Step 5: Visualization completed. Check the results directory for plots and CSV files.")
