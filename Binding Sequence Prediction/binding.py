# -*- coding: utf-8 -*-
"""

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1cbtbvCw1BkT9purE9WJn0Vow_tSQ-uzD
"""

from Bio import SeqIO
from sklearn.feature_extraction.text import CountVectorizer

def read_fasta(file):
    sequences = []
    for record in SeqIO.parse(file, "fasta"):
        sequences.append(str(record.seq).upper())
    return sequences

def kmer_counting(sequences, k=2):
    cv = CountVectorizer(ngram_range=(k, k), analyzer='char')
    kmer_features = cv.fit_transform(sequences)
    return kmer_features

# Read the sequences
bound_seqs = read_fasta("bound.fasta")
unbound_seqs = read_fasta("notbound.fasta")
test_seqs = read_fasta("test.fasta")


# Convert sequences to k-mer features
X_bound = kmer_counting(bound_seqs)
X_unbound = kmer_counting(unbound_seqs)
X_test = kmer_counting(test_seqs)


# Labels
y_bound = [1] * X_bound.shape[0]  # 1 for bound
y_unbound = [0] * X_unbound.shape[0]  # 0 for unbound

# Combine bound and unbound for training
import numpy as np
X_train = np.vstack((X_bound.toarray(), X_unbound.toarray()))
y_train = np.array(y_bound + y_unbound)

from sklearn.ensemble import RandomForestClassifier

# Initialize and train the classifier
clf = RandomForestClassifier(n_estimators=100, random_state=42)
clf.fit(X_train, y_train)

# Predict probabilities for the test set
predictions = clf.predict_proba(X_test)[:, 1]  # Get probability for the 'bound' class

# Get the top 6000 sequences
import csv
import zipfile


top_indices = np.argsort(predictions)[::-1][:6000]
top_seq_names = [f'seq{i+1}' for i in top_indices]

csv_filename = "predictions.csv"

# Write the sequence names to the CSV file
with open(csv_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for seq_name in top_seq_names:
        writer.writerow([seq_name])  # Pass seq_name as a single-element list

zip_filename = "electric_relaxationerssssss.zip"

# Create a Zip file and add the CSV file to it
with zipfile.ZipFile(zip_filename, 'w') as zipf:
    zipf.write(csv_filename)
