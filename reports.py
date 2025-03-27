import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#more time reqiured to tidy code!!!
# read in output csv from mapping script, add in missing uniprot codes
csv_file = "output_fungiDB_with_pLDDT.csv"  
df = pd.read_csv(csv_file)
df['uniprotID_extracted'] = df['gene_product'].str.extract(r';Acc:(\S+)')
df.loc[df['uniprotID'].isna() & df['gene_product'].str.contains('UniProtKB', na=False), 'uniprotID'] = df['uniprotID_extracted']
print(df)

# Filter data based on the presence of uniprotID
with_uniprot = df[df['uniprotID'].notna()]['mean_pLDDT']
without_uniprot = df[df['uniprotID'].isna()]['mean_pLDDT']

# Create histogram
plt.figure(figsize=(8, 5))
plt.hist([with_uniprot, without_uniprot], bins=20, color=['blue', 'red'], edgecolor='black', alpha=0.7, label=['With UniProt ID', 'No UniProt ID'])
plt.xlabel('Mean pLDDT')
plt.ylabel('Count')
plt.title('Histogram of Mean pLDDT')
plt.legend(loc='upper left')
plt.savefig("mean_plddt_histogram.png", dpi=300, bbox_inches='tight')

# Three Filters
filtered_none_df = df[(df['uniprotID'].isna()) | (df['uniprotID'] == 'None')]
filtered_df = df[(df['uniprotID'].notna())]
dubious_rows = df[df["gene_product"].str.contains("Dubious", case=False, na=False)]

# Basic Stats on the three groups
mean_protein_length_all = round(df["protein_length"].mean())
mean_protein_length_none = round(filtered_none_df["protein_length"].mean())
mean_protein_length_uniprot = round(filtered_df["protein_length"].mean())
mean_protein_length_dubious = round(dubious_rows["protein_length"].mean())

mean_pLDDT_all = round(df["mean_pLDDT"].mean())
mean_pLDDT_none = round(filtered_none_df["mean_pLDDT"].mean())
mean_pLDDT_uniprot = round(filtered_df["mean_pLDDT"].mean())
mean_pLDDT_dubious = round(dubious_rows["mean_pLDDT"].mean())

exon_count_all = df["no_exons"].value_counts()
exon_count_none = filtered_none_df["no_exons"].value_counts()
exon_count_uniprot = filtered_df["no_exons"].value_counts()
exon_count_dubious = dubious_rows["no_exons"].value_counts()

count_rows_above_70_all = (df["plddt_percentage"] > 70).sum()
count_rows_above_70_none = (filtered_none_df["plddt_percentage"] > 70).sum()
count_rows_above_70_uniprot = (filtered_df["plddt_percentage"] > 70).sum()
count_rows_above_70_dubious = (dubious_rows["plddt_percentage"] > 70).sum()

# Create a DataFrame to print all the statistics in a table
stats_data = {
    "Category": ["All", "Non-Uniprot", "Uniprot", "Dubious ORF"],
    "Mean Protein Length": [mean_protein_length_all, mean_protein_length_none, mean_protein_length_uniprot, mean_protein_length_dubious],
    "Mean pLDDT": [mean_pLDDT_all, mean_pLDDT_none, mean_pLDDT_uniprot, mean_pLDDT_dubious],
    "Exon Count Distribution": [exon_count_all.to_dict(), exon_count_none.to_dict(), exon_count_uniprot.to_dict(), exon_count_dubious.to_dict()],
    "Rows with pLDDT > 70%": [count_rows_above_70_all, count_rows_above_70_none, count_rows_above_70_uniprot, count_rows_above_70_dubious],
}

# Filter Non-Uniprot category with Rows where pLDDT > 70%
non_uniprot_above_70 = filtered_none_df[filtered_none_df["plddt_percentage"] > 70]
print("Top 4 rows for 'Non-Uniprot' category where Rows with pLDDT > 70%:")
print(non_uniprot_above_70.head(4))
#test different thesholds
non_uniprot_above_70 = filtered_none_df[filtered_none_df["plddt_percentage"] > 50]
non_uniprot_above_70 = filtered_none_df[filtered_none_df["plddt_percentage"] > 20]


# Sort the DataFrame by 'plddt_percentage' in descending order
non_uniprot_above_70_sorted = non_uniprot_above_70.sort_values(by="plddt_percentage", ascending=False)
non_uniprot_above_70_sorted = non_uniprot_above_70_sorted[['fungiDBID', 'gene_product', 'plddt_percentage', 'protein_length']]
print("All rows for 'Non-Uniprot' category where pLDDT > 50%, sorted by pLDDT percentage in descending order:")
print(non_uniprot_above_70_sorted)

dubious_rows = df.loc[df["gene_product"].str.contains("Dubious", case=False, na=False), 
                      ["fungiDBID", "gene_product", "plddt_percentage", "protein_length"]]
print(dubious_rows)

# Convert the stats data into a DataFrame for easier easier formattin in excell for final report
stats_df = pd.DataFrame(stats_data)
print("Basic Statistics Table:")
print(stats_df)
stats_df.to_csv("summary_statistics.csv", index=False)
print("Summary statistics saved as 'summary_statistics.csv'")

#********************violoin plots****************************
# Create a new DataFrame for visualization
plot_df = pd.concat([
    df.assign(Category="All"),
    filtered_none_df.assign(Category="Non-Uniprot"),
    filtered_df.assign(Category="Uniprot"),
    dubious_rows.assign(Category="Dubious ORF")
])

# Remove rows with protein length > 1000 for the protein length plot
plot_df_filtered = plot_df[plot_df["protein_length"] <= 1000]

sns.set(style="whitegrid")

# Violin plot for mean pLDDT
plt.figure(figsize=(7, 6))
sns.violinplot(x="Category", y="mean_pLDDT", data=plot_df, inner="box")
plt.title("Distribution of Mean pLDDT")
plt.ylabel("Mean pLDDT")
plt.xlabel("Category")
plt.xticks(rotation=20, ha="right")
plt.tight_layout()
plt.savefig("violin_mean_pLDDT.png", dpi=300)
plt.close()

# Violin plot for protein length (filtered - removed prot length >3000)
plt.figure(figsize=(7, 6))
sns.violinplot(x="Category", y="protein_length", data=plot_df_filtered, inner="box")
plt.title("Distribution of Mean Protein Length (Filtered ≤3000)")
plt.ylabel("Protein Length")
plt.xlabel("Category")
plt.xticks(rotation=20, ha="right")
plt.tight_layout()
plt.savefig("violin_protein_length.png", dpi=300)
plt.close()


