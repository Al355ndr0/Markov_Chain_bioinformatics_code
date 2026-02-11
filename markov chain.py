import math
import random

#ALESSANDRO ZIRUOLO 



NUCLEOTIDES = ['A', 'C', 'G', 'T']

#genome loading from the FASTA file 
def fasta(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    return ''.join(line.strip() for line in lines if not line.startswith('>')).upper()

#loading of CpG coordinates
def cpg_coordinates(filename, chromosome='chr22'):
    coordinates = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if parts[0] == chromosome:
                start = int(parts[1])
                end = int(parts[2])
                coordinates.append((start, end))
    return coordinates

def sequence_extraction(genome, coords):
    sequences = []
    for start, end in coords:
        seq = genome[start:end]
        if 'N' not in seq:
            sequences.append(seq)
    return sequences

def sequence_extraction_outside_cpg(genome, coords):
    lengths = [end - start for start, end in coords]
    genome_len = len(genome)
    sequences = []
    for length in lengths:
        for i in range(10):   #number of tries to find a valid sequence without 'N'
            start = random.randint(0, genome_len - length - 1)
            end = start + length
            seq = genome[start:end]
            if 'N' not in seq:
                sequences.append(seq)
                break
    return sequences

#Matrix initialization and dinucleotide counting
def init_matrix():
    matrix = {}  
    for x in NUCLEOTIDES:
        matrix[x] = {}
        for y in NUCLEOTIDES:
            matrix[x][y] = 0
    return matrix

def count_dinucleotides(sequences):
    matrix = init_matrix()
    for seq in sequences:
        for i in range(len(seq) - 1):
            x, y = seq[i], seq[i + 1]
            if x in NUCLEOTIDES and y in NUCLEOTIDES:
                matrix[x][y] += 1
    return matrix

def probabilities(matrix):
    probs = {}
    for x in NUCLEOTIDES:
        total = 0
        for y in matrix[x]:
            total += matrix[x][y]
        probs[x] = {}
        for y in NUCLEOTIDES:
            probs[x][y] = (matrix[x][y] + 1) / (total + 4)  # pseudocount
    return probs

#Calculate log probabilities for a sequence
def log_probability(seq, probs):
    log_prob = math.log(0.25)  # Uniform base prior
    for i in range(1, len(seq)):
        x, y = seq[i - 1], seq[i]
        if x in probs and y in probs[x]:
            log_prob += math.log(probs[x][y])
        else:
            log_prob += math.log(1e-10)
    return log_prob

def score_sequence(seq, inside_probs, outside_probs):
    return log_probability(seq, inside_probs) - log_probability(seq, outside_probs)


# Main execution
genome = fasta("chr22.fa")
cpg_coords = cpg_coordinates("model-based-cpg-islands-hg19.txt")

inside_seqs = sequence_extraction(genome, cpg_coords)
outside_seqs = sequence_extraction_outside_cpg(genome, cpg_coords)

inside_counts = count_dinucleotides(inside_seqs)
outside_counts = count_dinucleotides(outside_seqs)

inside_probs = probabilities(inside_counts)
outside_probs = probabilities(outside_counts)

query = input("Enter a DNA sequence (100 nt) to score: ").strip().upper()
if len(query) < 2:
    print("Sequence too short.")
else:
    score = score_sequence(query, inside_probs, outside_probs)
    print(f"\nLog-ratio score S(X): {score:.4f}")
    if score > 0:
        print("Likely CpG island ")
    else:
        print("Not a CpG island ")



#to run this code from the TERMINAL, open it from the directory where this file is located and copy the following command: python3 markov\ chain.py 
#make sure to have the files chr22.fa and model-based-cpg-islands-hg19.txt in the same directory

#otherwise make sure to open the code from the folder containing also chr22.fa and model-based-cpg-islands-hg19.txt