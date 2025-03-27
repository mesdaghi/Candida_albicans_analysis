import os
import pandas as pd

# Construct dictionary to hold UniProt mapping v50
mapping = {}
with open('FungiDB-50_CalbicansSC5314_UniProtMapping.txt', 'r') as fungi_mapping_file:
    for line in fungi_mapping_file:
        temp = line.split()
        mapping[temp[0]] = temp[1]

# Construct dictionary for exon counts
gff_file = "FungiDB-68_CalbicansSC5314.gff"
exon_counts = {}
with open(gff_file, 'r') as file:
    for line in file:
        if line.startswith("#") or not line.strip():
            continue
        columns = line.strip().split("\t")
        if len(columns) < 9:
            continue
        feature_type = columns[2]
        attributes = columns[8]
        if feature_type == "exon":
            attr_dict = {key: value for key, value in 
                         (field.split("=") for field in attributes.split(";") if "=" in field)}
            gene_id = attr_dict.get("gene_id")
            if gene_id:
                exon_counts[gene_id] = exon_counts.get(gene_id, 0) + 1

# Read in FungiDB data v68
fungi_prot_lines = []
with open('FungiDB-68_CalbicansSC5314_AnnotatedProteins.fasta', 'r') as fungi_prot_file:
    for entry in fungi_prot_file.read().split('>C'):
        if len(entry) > 1:
            fungi_prot_lines.append('C' + entry)

# Get pandas column headers
titles_list = ['fungiDBID']
header_list = fungi_prot_lines[1].split('|')
titles_list.extend([x.split('=')[0].strip() for x in header_list[1:]])
titles_list.extend(['uniprotID', 'no_exons', 'mean_pLDDT', 'CA_count_plddt_80'])
df = pd.DataFrame(columns=titles_list)

# Populate pandas DataFrame
for entry in fungi_prot_lines:
    db_entry = entry.split('|')
    row = [db_entry[0].strip()]
    row.extend(x.split('=')[1].strip() for x in db_entry[1:])
    
    # UniProt ID mapping
    try:
        uniprot_code = list(mapping.keys())[list(mapping.values()).index(db_entry[0].replace(" ", "").rstrip(db_entry[0][-6:]))]
    except ValueError:
        uniprot_code = None
    row.append(uniprot_code)

    # Exon count mapping
    try:
        exon_number = exon_counts[db_entry[0].replace(" ", "").rstrip(db_entry[0][-6:])]
    except KeyError:
        exon_number = None
    row.append(exon_number)

    # Append row to DataFrame
    # Ensure row has the same number of elements as df.columns
    while len(row) < len(df.columns):
        row.append(None)  # Fill missing values with None

    df.loc[len(df)] = row  # Now, row length matches df.columns

    df.loc[len(df)] = row

# Function to extract pLDDT from PDB files
def extract_plddt_from_pdb(pdb_file):
    ca_plddt_values = []
    with open(pdb_file, "r") as file:
        for line in file:
            if line.startswith("ATOM") and " CA " in line[12:16]:  # Select only Cα atoms
                plddt = float(line[61:66])  # pLDDT is in columns 62-66
                ca_plddt_values.append(plddt)
    return ca_plddt_values

# Process PDB files and update DataFrame
pdb_directory = "./models"
for pdb_file in os.listdir(pdb_directory):
    if pdb_file.endswith(".pdb"):
        file_path = os.path.join(pdb_directory, pdb_file)
        plddt_values = extract_plddt_from_pdb(file_path)

        if plddt_values:
            mean_plddt = sum(plddt_values) / len(plddt_values)
            high_plddt_count = sum(1 for p in plddt_values if p > 80)
        else:
            mean_plddt = None
            high_plddt_count = None

        # Determine row identifier
        if pdb_file.startswith("C"):
            row_id_column = "gene"
        else:
            row_id_column = "uniprotID"

        # Update DataFrame if the row exists
        row_index = df[df[row_id_column] == pdb_file.rstrip(".pdb")].index
        if not row_index.empty:
            df.loc[row_index, "mean_pLDDT"] = mean_plddt
            df.loc[row_index, "CA_count_plddt_80"] = high_plddt_count

df = df.drop_duplicates() # duplicate rows for some reason!!!
df["protein_length"] = pd.to_numeric(df["protein_length"], errors="coerce")
df["plddt_percentage"] = (df["CA_count_plddt_80"] / df["protein_length"]) * 100

print(df.to_string())
print(df)
df.to_csv("output_fungiDB_with_pLDDT.csv", index=False)
