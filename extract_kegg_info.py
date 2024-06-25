import argparse  
from Bio.SeqRecord import SeqRecord  
from Bio.Seq import Seq

def parse_kegg_file(kegg_file, organism):
    gene_records = []  #list to store SeqRecord objects representing gene sequences
    seen_genes = set()  # Set to store seen gene names to avoid duplicates
    
    with open(kegg_file, 'r') as f:
        gene_name = None
        gene_sequence = None
        in_ntseq_section = False  #flag to track if currently inside the NTSEQ section
        for line in f:
            if line.startswith('NAME'):
                gene_name = line.strip().split(maxsplit=1)[1]  #extract gene name
            elif line.startswith('NTSEQ'):
                in_ntseq_section = True  #mark beginning of sequence section
                gene_sequence = ''  #initialize gene sequence string
            elif line.startswith('///'):
                #end of gene entry, create SeqRecord and add to gene_records if not seen before
                if gene_name and gene_sequence and gene_name not in seen_genes:
                    #create SeqRecord with gene sequence and metadata
                    record_title = f">{gene_name} ({organism})"
                    gene_record = SeqRecord(Seq(gene_sequence), id=gene_name, description=record_title)
                    gene_records.append(gene_record)  #add SeqRecord to the list
                    seen_genes.add(gene_name)  #add gene name to the set of seen genes
                gene_name = None  #reset gene name
                gene_sequence = None  #reset gene sequence
                in_ntseq_section = False  #mark end of sequence section
            elif in_ntseq_section:
                gene_sequence += line.strip()  #append sequence lines to gene sequence string
    
    return gene_records  #return list of SeqRecord objects representing gene sequences

def write_fasta_file(gene_records, fasta_file):
    #write gene records to a FASTA file
    with open(fasta_file, 'w') as f:
        for record in gene_records:
            f.write(f"{record.description}\n{record.seq}\n")  # Write record metadata and sequence to file

if __name__ == "__main__":
    #command-line argument parsing
    parser = argparse.ArgumentParser(description="Parse a KEGG file and generate a FASTA file")
    parser.add_argument("kegg_file", help="Input KEGG file")  # Input KEGG file path
    parser.add_argument("organism_name", help="Organism name")  # Organism name
    parser.add_argument("-o", "--output", help="Output FASTA file", required=True)  # Output FASTA file path
    args = parser.parse_args()

    #parse KEGG file and generate SeqRecord objects
    gene_records = parse_kegg_file(args.kegg_file, args.organism_name)
    
    #write to a FASTA file
    write_fasta_file(gene_records, args.output)
