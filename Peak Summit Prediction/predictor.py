


import numpy as np
from collections import Counter
import csv
import zipfile

def create_initial_pwm(seqs, k):
    """Create an initial PWM with uniform probabilities."""
    pwm = np.ones((4, k)) / 4  # Initialize with uniform probabilities
    return pwm

def seq_to_numeric(seq):
    """Convert a DNA sequence to a numeric array."""
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return np.array([mapping[s] for s in seq])

def update_pwm(seqs, motif_starts, k):
    """Update the PWM based on the current motif positions."""
    counts = np.zeros((4, k))
    for i, seq in enumerate(seqs):
        start = motif_starts[i]
        for j in range(k):
            nucleotide = seq_to_numeric(seq[start:start+k])[j]
            counts[nucleotide, j] += 1
    pwm = (counts + 1) / (len(seqs) + 4)  # Laplace smoothing
    return pwm

def expectation(seqs, pwm, k):
    """E-step: Compute the probability of each k-mer being the motif."""
    probs = []
    for seq in seqs:
        seq_numeric = seq_to_numeric(seq)
        prob = []
        for i in range(len(seq) - k + 1):
            k_mer = seq_numeric[i:i+k]
            p = np.product([pwm[nucleotide, j] for j, nucleotide in enumerate(k_mer)])
            prob.append(p)
        probs.append(prob)
    return probs

def maximization(seqs, probs, k):
    """M-step: Find the most likely start positions of motifs."""
    motif_starts = []
    for prob in probs:
        motif_starts.append(np.argmax(prob))
    return motif_starts

def em_motif_finding(seqs, k, max_iter=200):
    """EM algorithm for motif finding."""
    # Initialize with random positions
    motif_starts = np.random.choice(len(seqs[0]) - k + 1, len(seqs))
    pwm = create_initial_pwm(seqs, k)

    for iteration in range(max_iter):
        print(iteration)
        # E-step
        probs = expectation(seqs, pwm, k)
        # M-step
        motif_starts = maximization(seqs, probs, k)
        # Update PWM
        pwm = update_pwm(seqs, motif_starts, k)

    return pwm

def calculate_pwm_score(pwm, k_mer_numeric):
    """Calculate the score of a k-mer based on the PWM."""
    score = 0
    for i, nucleotide in enumerate(k_mer_numeric):
        score += np.log(pwm[nucleotide, i])  # Using log probabilities
    return score

def find_motif_positions_in_sequences(pwm, seqs, k):
    """Find positions of motifs in a list of sequences based on the PWM."""
    motif_positions = []
    for seq in seqs:
        best_score = -np.inf
        best_position = None
        for i in range(len(seq) - k + 1):
            k_mer = seq[i:i+k]
            k_mer_numeric = seq_to_numeric(k_mer)
            score = calculate_pwm_score(pwm, k_mer_numeric)
            if score > best_score:
                best_score = score
                best_position = i
        motif_positions.append(best_position)
    return motif_positions


def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        seq_id = None
        for line in lines:
            if line.startswith('>'):
                seq_id = line.strip()[1:]
                sequences[seq_id] = ''
            else:
                sequences[seq_id] += line.strip()
    return sequences

centered_sequences = read_fasta('boundcentered.fasta')
offset_sequences = read_fasta('boundrandomoffset.fasta')

pwm = em_motif_finding(list(centered_sequences.values()), 30)
motif_positions = find_motif_positions_in_sequences(pwm, list(offset_sequences.values()), 30)
result_dict = {'seq{}'.format(i+1): value for i, value in enumerate(motif_positions)}
csv_filename = "predictions.csv"


# Write the dictionary to the CSV file
with open(csv_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')
    for key, value in result_dict.items():
        writer.writerow([key, value])

zip_filename = "love_letter_to_you_3.zip"

# Create a Zip file and add the CSV file to it
with zipfile.ZipFile(zip_filename, 'w') as zipf:
    zipf.write(csv_filename)
