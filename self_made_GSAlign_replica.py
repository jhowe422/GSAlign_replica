import argparse
from Bio import SeqIO
#class to parse files since it is so much data
class GeneEntry:
    def __init__(self, gene_id=None, sequence=None, organism=None, name=None): #init
        self.gene_id = gene_id #get NCBI protein id
        self.sequence = sequence #get sequence
        self.organism = organism #organism name
        self.name = name #gene name
#parse files
def parse_fasta_files(file_paths):
    entries = [] #empty list to store entries

    for file_path in file_paths:
        print(f"Parsing file: {file_path}")
        for record in SeqIO.parse(file_path, "fasta"): #use seqio to parse the fiile
            entry = GeneEntry(gene_id=record.id, sequence=str(record.seq), organism=file_path.split('.')[0], name=record.description)
            entries.append(entry) #add to entry

    return entries
#GSAlign class to perform alignment and scoring
class GSAlign:
    def __init__(self): #init
        pass
#main to loop through
    def main(self, file_paths, output_file): 
        entries = parse_fasta_files(file_paths) #get entries from parse function

        with open(output_file, 'w') as f: #write everything to output file
            f.write("Genome Comparison Results:\n\n")

            for i, entry1 in enumerate(entries): #seq1
                for j, entry2 in enumerate(entries): #seq2
                    if i < j and entry1.organism != entry2.organism:  #compare genes only between different species
                        similarity_score = self.compute_similarity(entry1.sequence, entry2.sequence) #get sim score
                        #write to output file
                        f.write(f"Comparison between {entry1.gene_id} ({entry1.name}, {entry1.organism}) and "
                                f"{entry2.gene_id} ({entry2.name}, {entry2.organism}):\n")
                        f.write(f"Similarity Score: {similarity_score}\n\n")
#similarity score function
    def compute_similarity(self, seq1, seq2):
        #Ukkonen's edit distance as measure of similarity
        edit_distance = self.ukkonen_edit_distance(seq1, seq2)
        
        #normalize the edit distance to get sim score
        max_length = max(len(seq1), len(seq2))
        similarity_score = 100 - (edit_distance / max_length) * 100 #sim score
        
        return similarity_score
#get ukkonen dist
    def ukkonen_edit_distance(self, seq1, seq2):
        #initialize the matrix
        prev_row = list(range(len(seq2) + 1))
        current_row = [0] * (len(seq2) + 1)

        #compute the edit distance row by row
        for i in range(1, len(seq1) + 1): #for a base in the seq1
            current_row[0] = i #this is first row
            for j in range(1, len(seq2) + 1): #for base in seq2
                cost = 0 if seq1[i - 1] == seq2[j - 1] else 1 #calculate the cost to change to the corresponding base in 1
                current_row[j] = min(current_row[j - 1] + 1, prev_row[j] + 1, prev_row[j - 1] + cost) #update matrix
            prev_row = current_row.copy()

        return current_row[-1] #return matrix row
    #BWT
    def run_bwt_alignment(self, seq1, seq2):
        #uses seq1 and 2
        bwt_seq1 = self.bwt(seq1 + "$") 
        bwt_seq2 = self.bwt(seq2 + "$")
    #do this for the sequence that is best match
        matches = self.find_exact_matches(bwt_seq1, seq2)
        alignment_score = len(matches)

        return alignment_score
#bwt: create suffix array with the last base in the seq
#get indices based on the cyclic rotation with this array and the consecutive last character with shuffling
#joins characters using suffix array to gie BWT seq
    def bwt(self, sequence):
        sa = sorted(range(len(sequence)), key=lambda x: sequence[x:])
        bwt_sequence = ''.join(sequence[i - 1] for i in sa)
        return bwt_sequence
#get the matches using the bwt indices
    def find_exact_matches(self, bwt_sequence, query_sequence):
        #Find exact matches between a query sequence and a Burrows-Wheeler Transform (BWT) sequence.
        #bwt_sequence (str): The Burrows-Wheeler Transform (BWT) sequence.
        #query_sequence (str): The query sequence to search for in the BWT sequence.
        #Returns list of tuples containing the indices of exact matches found in the BWT sequence.

        matches = []  #list to store indices of exact matches
        top = 0  #initialize top index of the search range
        bottom = len(bwt_sequence) - 1  #initialize the bottom index of the search range

        #iterate through BWT sequence while search range is valid
        while top <= bottom:
            #check if query sequence is not empty
            if query_sequence:
                symbol = query_sequence[-1]  #get last symbol of the query sequence
                query_sequence = query_sequence[:-1]  #remove last symbol from the query sequence

                #check if the symbol is present in the current search range of the BWT sequence
                if symbol in bwt_sequence[top:bottom + 1]:
                    #find the index of the first occurrence of the symbol in the current search range
                    top_index = bwt_sequence[top:bottom + 1].index(symbol) + top
                    #find the index of the last occurrence of the symbol in the current search range
                    bottom_index = bwt_sequence[top:bottom + 1].rindex(symbol) + top
                    
                    #update the search range based on the occurrences of the symbol
                    top = self.first_occurrence(bwt_sequence, symbol) + self.count_occurrences(bwt_sequence, symbol, top)
                    bottom = self.first_occurrence(bwt_sequence, symbol) + self.count_occurrences(bwt_sequence, symbol, bottom + 1) - 1
                    
                    #add the indices of the exact match to the list of matches
                    matches.append((top_index, bottom_index))
                else:
                    return []  #if the symbol is not found, return an empty list
            else:
                #if the query sequence is empty, add the current search range to the list of matches
                matches.append((top, bottom))
                top += 1  #increment the top index to continue the search

        return matches  #return list of matches

#get first occurrence
    def first_occurrence(self, bwt_sequence, symbol):
        return bwt_sequence.find(symbol)
#count occurrances
    def count_occurrences(self, bwt_sequence, symbol, index):
        return bwt_sequence[:index].count(symbol)
#main to initialize parsers and run everything
def main():
    parser = argparse.ArgumentParser(description='Perform genome comparison using GSAlign.') 
    parser.add_argument('files1', metavar='F1', type=str, nargs='+',
                        help='input FASTA files for the first genome') #input file 1
    parser.add_argument('files2', metavar='F2', type=str, nargs='+',
                        help='input FASTA files for the second genome') #input file 2
    parser.add_argument('-o', '--output', metavar='output_file', type=str,
                        help='output file to save comparison results') #output file

    args = parser.parse_args() #run args

    gsalign = GSAlign()
    gsalign.main(args.files1 + args.files2, args.output)

if __name__ == "__main__":
    main()
