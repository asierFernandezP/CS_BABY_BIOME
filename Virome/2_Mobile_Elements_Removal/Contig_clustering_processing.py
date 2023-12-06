import argparse

# Define command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('Clusters', help='input file with clusters of contigs')
parser.add_argument('Plasmid_sequences', help='input file with identified plasmid sequences')
args = parser.parse_args()

#Read the files
with open(args.Clusters, 'r') as clusters, open(args.Plasmid_sequences, 'r') as plasmids, open('Viral_sequences_no_plasmids.txt', 'w') as output:
    # Read the sequence names to exclude from plasmids file into a set
    exclude = set(line.strip() for line in plasmids)

# Loop over the lines in clusters file:
    for line in clusters:
        # Split the line by tabs and discard the first column (representative sequence - it is already included in the 2nd column)
        cols = line.strip().split('\t')[1:]
        # Split each sequence by comma
        seqs = [seq.strip() for col in cols for seq in col.split(',')]
        # Check if any of the sequences in the line are in the exclude set
        if any(seq in exclude for seq in seqs):
            # If so, skip the line
            continue
        # If no sequence is present in the exclude set, write the sequences to output file
        for seq in seqs:
            output.write(seq + '\n')
