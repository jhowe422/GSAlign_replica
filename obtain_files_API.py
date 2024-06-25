#replace search term with your desired term
import requests
import csv
from Bio import Entrez

#get the info from kegg: Kyoto Encyclopedia of Genes and Genomes
#for understanding functions of bio systems from genome sequencing and other high-throughput technologies
def get_kegg_gene_info(organism_code):
    url = f"http://rest.kegg.jp/list/{organism_code}" #url
    response = requests.get(url)
    if response.status_code == 200: #success
        return response.text.splitlines()
    else:
        print(f"Failed to retrieve gene information from KEGG for {organism_code}") #error message
        return []

# Function to retrieve gene annotations from KEGG
def get_kegg_gene_annotations(gene_id):
    url = f"http://rest.kegg.jp/get/{gene_id}" #kegg url
    response = requests.get(url)
    if response.status_code == 200: #success
        return response.text
    else:
        print(f"Failed to retrieve annotations from KEGG for {gene_id}") #print error message
        return ""

#retrieve gene information from NCBI Entrez
def get_ncbi_gene_info(term, db="gene", retmax=10):
    Entrez.email = ""  #your email here
    handle = Entrez.esearch(db=db, term=term, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

#retrieve gene name and annotations from NCBI Entrez
def get_ncbi_gene_details(gene_id):
    handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml")
    record = Entrez.read(handle)
    handle.close()
#initialize gene name and annotations info
    gene_name = ""
    annotations = ""
#get gene locus
    if record:
        if "Entrezgene_locus" in record[0]:
            gene_locus = record[0]["Entrezgene_locus"]
            if gene_locus:
                gene_name = gene_locus[0].get("Gene-commentary_products", [{}])[0].get("Gene-commentary_label", "")
#get gene annotations
        if "Entrezgene_commentary" in record[0]:
            commentary = record[0]["Entrezgene_commentary"]
            if commentary:
                for comment in commentary:
                    if "Gene-commentary_comment" in comment:
                        annotations += comment["Gene-commentary_comment"] + "\n"

    return gene_name, annotations.strip()


#usage
term = "" #your search term here
ncbi_gene_ids = get_ncbi_gene_info(term) #run ncbi function

# Write gene information and annotations to a file
output_file = "simulans.tsv" #output file
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile, delimiter="\t") #use csv writer
    writer.writerow(["Source", "Gene ID", "Gene Name", "Annotations"]) #format output table
    
    #kegg info
    organism_code = " "  #organism code of choice here
    kegg_gene_info = get_kegg_gene_info(organism_code) #kegg function
    for line in kegg_gene_info: #extrapolate info from kegg
        parts = line.split("\t")
        gene_id = parts[0] if parts else ""
        gene_name = parts[1] if len(parts) > 1 else ""
        annotations = get_kegg_gene_annotations(gene_id)
        writer.writerow(["KEGG", gene_id, gene_name, annotations])
    
    #get gene information and annotations from NCBI Entrez
    for gene_id in ncbi_gene_ids:
        gene_name, annotations = get_ncbi_gene_details(gene_id) #gene details function
        if annotations:
            if len(annotations) >= 3 and annotations[2] == gene_name:
                writer.writerow(["NCBI", gene_id, gene_name, annotations])
        else:
            writer.writerow(["NCBI", gene_id, gene_name, annotations])

print(f"Data has been written to {output_file}") #completion message
