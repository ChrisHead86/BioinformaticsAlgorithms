import sys
from typing import List, Dict, Iterable
import csv
import zipfile
import os

sys.setrecursionlimit(1000000)

def de_bruijn_kmers(k_mers: List[str]) -> Dict[str, List[str]]:
    """Forms the de Bruijn graph of a collection of k-mers."""
    de_bruijn_graph = {}

    for kmer in k_mers:
        prefix = kmer[:-1]
        suffix = kmer[1:]

        if prefix not in de_bruijn_graph:
            de_bruijn_graph[prefix] = []
        de_bruijn_graph[prefix].append(suffix)

    return de_bruijn_graph

def eulerian_path(g: Dict[str, List[str]]) -> Iterable[str]:
    graph = g.copy()
    path = []
    in_degree = {}
    out_degree = {}
    edges = 0


    n = len(graph) + 10

    def setUp():
        nonlocal in_degree, out_degree, edges
        in_degree = {node: 0 for node in graph}
        out_degree = {node: 0 for node in graph}
        for node, neighbors in graph.items():
            for neighbor in neighbors:
                if neighbor not in in_degree:
                    in_degree[neighbor] = 0
                if neighbor not in out_degree:
                    out_degree[neighbor] = 0
        edges = 0

        for from_node in graph:
            for to_node in graph[from_node]:
                in_degree[to_node] += 1
                out_degree[from_node] += 1
                edges += 1

    def hasEulerianPath():
        return True
        if edges == 0:
            return False

        start_nodes = 0
        end_nodes = 0

        for node in in_degree:
            if out_degree[node] - in_degree[node] > 1 or in_degree[node] - out_degree[node] > 1:
                return False
            elif out_degree[node] - in_degree[node] == 1:
                start_nodes += 1
            elif in_degree[node] - out_degree[node] == 1:
                end_nodes += 1

        

        return (end_nodes == 0 and start_nodes == 0) or (end_nodes == 1 and start_nodes == 1)

    def findStartNode():
        start = 0

        for node in graph:
            if out_degree[node] - in_degree[node] == 1:
                return node
            if out_degree[node] > 0:
                start = node

        return start

    def dfs(at):
        nonlocal path
        while out_degree[at]:
            next_node = graph[at].pop()
            out_degree[at] -= 1
            dfs(next_node)
        path.append(at)

    def getEulerianPath():
        setUp()

        if not hasEulerianPath():
            return []

        dfs(findStartNode())

        if len(path) != edges + 1:
            return []

        return path[::-1]

    return getEulerianPath()

def PathToGenome(path: Iterable[str]) -> str:
    """Constructs a genome string from a given Eulerian path."""
    genome = path[0]
    for pattern in path[1:]:
        genome += pattern[-1]
    return genome

def string_reconstruction(patterns: List[str], k: int) -> List[str]:
    """Reconstructs the order of the reads from a list of k-mers."""
    dB = de_bruijn_kmers(patterns)
    path = eulerian_path(dB)
    return path

def read_fasta(file_name):
    sequences = {}
    with open(file_name, 'r') as file:
        current_id = None
        current_sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_sequence] = current_id
                current_id = line[1:]
                current_sequence = ''
            else:
                current_sequence += line
        if current_id is not None:
            sequences[current_sequence] = current_id
    return sequences


def count_kmers(reads, k):
    kmer_counts = {}
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    return kmer_counts

def get_node_degrees(adjacency_list):
    node_degrees = {}
    for outgoing_node, incoming_nodes in adjacency_list.items():
        
        if outgoing_node in node_degrees:
            node_degrees[outgoing_node][1] += len(incoming_nodes)
        else:
            node_degrees[outgoing_node] = [0, len(incoming_nodes)]

        
        for incoming_node in incoming_nodes:
            if incoming_node in node_degrees:
                node_degrees[incoming_node][0] += 1
            else:
                node_degrees[incoming_node] = [1, 0]

    return node_degrees

def find_paths(adjacency_list):

    
    node_degrees = get_node_degrees(adjacency_list)

    

    
    paths = []

    
    for v in list(node_degrees.keys()):
        if node_degrees[v] != [1,1]:
            if node_degrees[v][1] > 0:
                for w in adjacency_list[v]:
                    non_branching_path = [v,w]
                    while node_degrees[w] == [1,1]:
                        u = adjacency_list[w][0]
                        non_branching_path.append(u)
                        w = u
                    paths.append(non_branching_path)

    for path in paths:
        for node in path:
            if node in adjacency_list:
                del adjacency_list[node]

    
    for node in list(adjacency_list.keys()):
        if node_degrees[node] != [1,1]:
            del adjacency_list[node]


    while adjacency_list:
        start_node = list(adjacency_list.keys())[0]
        
        curr_node = start_node
        next_node = adjacency_list[start_node][0]

        cycle = [start_node]

        first_time_visiting_start_node = True
        while curr_node != start_node or first_time_visiting_start_node:
            first_time_visiting_start_node = False
            del adjacency_list[curr_node] 
            cycle.append(next_node)
            curr_node = next_node
            if next_node not in adjacency_list:
                continue
            next_node = adjacency_list[next_node][0]
    
        paths.append(cycle)

    return paths

def make_contigs(paths):
    contigs = []

    for path in paths:
        contig = path[0] 
        for i in range(1, len(path)):
            contig += path[i][-1] 
        contigs.append(contig)

    return contigs


def assemble_genome(contigs, paired_end_reads):
    # Concatenate the contigs into a single sequence
    assembled_genome = ''.join(contigs)
    
    # If paired-end reads are provided, use them to scaffold contigs
    if paired_end_reads:
        for read_pair in paired_end_reads:
            read1, read2 = read_pair
            orientation1, orientation2 = read1[1], read2[1]

            if orientation1 == "+" and orientation2 == "-":
                contig1, contig2 = read1[0], read2[0]
            elif orientation1 == "-" and orientation2 == "+":
                contig1, contig2 = read2[0], read1[0]
            else:
                continue  # Skip read pairs with inconsistent orientations
            
            # Find the index of the end of contig1 in the assembled genome
            index_contig1_end = assembled_genome.find(contig1[-10:])
            
            # Find the index of the start of contig2 in the assembled genome
            index_contig2_start = assembled_genome.find(contig2[:10])
            
            if index_contig1_end != -1 and index_contig2_start != -1:
                # Add the sequence between the end of contig1 and the start of contig2
                assembled_genome += contig2[(index_contig2_start + 10):]
                
    return assembled_genome

def pairs_to_tuples(reads_dict):
    read_pairs = []
    paired_reads = {}

    # Group reads by identifier (excluding /1 or /2 suffix)
    for sequence, name in reads_dict.items():
        identifier = name.split('/')[0]
        if identifier not in paired_reads:
            paired_reads[identifier] = {'/1': None, '/2': None}
        if '/1' in name:
            paired_reads[identifier]['/1'] = sequence
        elif '/2' in name:
            paired_reads[identifier]['/2'] = sequence

    # Construct tuples for each pair of reads
    for identifier, reads in paired_reads.items():
        if reads['/1'] is not None and reads['/2'] is not None:
            read_pairs.append(((reads['/1'], '+'), (reads['/2'], '-')))

    return read_pairs

def assemble_genome(contigs, paired_end_reads):
    # Concatenate the contigs into a single sequence
    assembled_genome = ''.join(contigs)
    
    # If paired-end reads are provided, use them to scaffold contigs
    if paired_end_reads:
        for read_pair in paired_end_reads:
            read1, read2 = read_pair
            orientation1, orientation2 = read1[1], read2[1]

            if orientation1 == "+" and orientation2 == "-":
                contig1, contig2 = read1[0], read2[0]
            elif orientation1 == "-" and orientation2 == "+":
                contig1, contig2 = read2[0], read1[0]
            else:
                continue  # Skip read pairs with inconsistent orientations
            
            # Find the index of the end of contig1 in the assembled genome
            index_contig1_end = assembled_genome.find(contig1[-10:])
            
            # Find the index of the start of contig2 in the assembled genome
            index_contig2_start = assembled_genome.find(contig2[:10])
            
            if index_contig1_end != -1 and index_contig2_start != -1:
                # Add the sequence between the end of contig1 and the start of contig2
                assembled_genome += contig2[(index_contig2_start + 10):]
                
    return assembled_genome


def sort_reads(reads_dict, assembled_genome, max_errors=3):
    read_positions = {}  # Dictionary to store positions of reads in the genome

    # Iterate through each read in the dictionary
    for read_seq, read_name in reads_dict.items():
        print(read_name)
        # Initialize the minimum edit distance and position
        min_edit_distance = float('inf')
        best_position = None

        # Iterate through all possible starting positions in the assembled genome
        for pos in range(len(assembled_genome) - len(read_seq) + 1):
            # Calculate the edit distance between the read and the substring of the assembled genome
            edit_distance = sum(1 for i in range(len(read_seq)) if read_seq[i] != assembled_genome[pos + i])

            # Update the minimum edit distance and best position if this position is better
            if edit_distance < min_edit_distance:
                min_edit_distance = edit_distance
                best_position = pos

        # If the minimum edit distance is within the specified maximum errors, store the position
        if min_edit_distance <= max_errors:
            read_positions[read_name] = [best_position]

    # Sort the reads based on their positions in the genome
    sorted_reads = sorted(read_positions.items(), key=lambda x: x[1])

    # Return the sorted reads
    return [read[0] for read in sorted_reads]
"""

def sort_reads(reads_dict, assembled_genome):
    read_positions = {}  # Dictionary to store positions of reads in the genome

    # Iterate through each read in the dictionary
    for read_seq, read_name in reads_dict.items():
        
        # Find all occurrences of the read in the assembled genome
        positions = [pos for pos in range(len(assembled_genome)) if assembled_genome.startswith(read_seq, pos)]

        # Store the positions of the read in the genome if positions exist
        if positions:
            read_positions[read_name] = positions

    # Sort the reads based on their positions in the genome
    sorted_reads = sorted(read_positions.items(), key=lambda x: x[1])

    # Return the sorted reads
    return [read for read, number in sorted_reads]
    """

def extract(sorted_reads):
    return [read_tuple[0] for read_tuple in sorted_reads]

def write_read_names_to_csv(read_names, csv_filename):
    with open(csv_filename, 'w') as csvfile:
        for read_name in read_names:
            print("Writing to CSV:", read_name)
            csvfile.write(">" + read_name + "\n")

def create_zip_file_with_csv(csv_filename, zip_filename):
    with zipfile.ZipFile(zip_filename, 'w') as zipf:
        zipf.write(csv_filename)


def main():
    file_name = "project2b_reads.fasta"
    sequences = read_fasta(file_name)
    new_file = "project2b_paired_reads.fasta"
    paired_reads = read_fasta(new_file)
    paired_reads = pairs_to_tuples(paired_reads)
    first_key = next(iter(sequences))
    kmer_size = len(first_key)
    input = [] 


    for key in sequences:
        input.append(key)

    k = 20
    kmer_count = count_kmers(input, k)
    kmer_count = {kmer: freq for kmer, freq in kmer_count.items() if freq > 3}

    de_bruijn = de_bruijn_kmers(kmer_count)

    new_paths = find_paths(de_bruijn)

    contigs = make_contigs(new_paths)

    assembled = assemble_genome(contigs, paired_reads)

    result = sort_reads(sequences, assembled)

    

    print(len(result))


    csv_filename = 'predictions.csv'
    zip_filename = 'project2b-this_one.zip'

    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for item in result:
            writer.writerow([f'>{item}'])
            print([f'>{item}'])

    with zipfile.ZipFile(zip_filename, 'w') as zipf:
        zipf.write(csv_filename)

    os.remove(csv_filename)

    print(f"CSV file '{csv_filename}' with '>' prepended has been zipped into '{zip_filename}'")




if __name__ == "__main__":
    main()
