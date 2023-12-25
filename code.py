import os
from Bio import SeqIO
import numpy as np

input_dir = r'C:\Users\Aseel\Desktop\four\pdb'
output_dir = r'C:\Users\Aseel\Desktop\four\pdb\out'

fasta_filename = os.path.join(output_dir, "fasta.fasta")
ss_filename = os.path.join(output_dir, "AA-SS.txt")

fasta_file = open(fasta_filename, "w")
ss_file = open(ss_filename, "w")

protein_lengths = []

for filename in os.listdir(input_dir):
    if filename.endswith('.ent'):
        pdb_path = os.path.join(input_dir, filename)

        ss_data = []
        with open(pdb_path) as handle:
            sequence = next(SeqIO.parse(handle, "pdb-atom"))

        SeqIO.write(sequence, fasta_file, "fasta")

        protein_lengths.append(len(sequence))

        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('HELIX'):
                    start_residue = line[21:25].strip()
                    end_residue = line[33:37].strip()
                    seq = sequence.seq[int(line[21:25].strip()) - 1:int(line[33:37].strip()) - 1]
                    ss_class = "h" * len(seq)
                    ss_data.append((start_residue, end_residue,ss_class, seq ))
                if line.startswith('SHEET'):
                    start_residue = line[22:26].strip()
                    end_residue = line[33:37].strip()
                    seq = sequence.seq[int(line[22:26].strip()) - 1:int(line[33:37].strip()) - 1]
                    ss_class = "b" * len(seq)
                    ss_data.append((start_residue, end_residue,ss_class, seq))

        count = 1
        for record in ss_data:
            ss_file.write(f"{count}\n{record[0]}:{record[1]}\n{record[2]}\n{record[3]}\n")
            count += 1

ent_files = {}
file_count = 0
with open("Summary.txt", "w") as file:
    for file_name in os.listdir(input_dir):
        file_path = os.path.join(input_dir, file_name)
        if os.path.isfile(file_path):
            file_ext = os.path.splitext(file_name)[1]
            if file_ext == ".ent":
                ent_files[file_name] = ent_files.get(file_name, 0) + 1
                file_count += 1

    file.write(f"Number of files : {file_count}\n")
    for file_name in ent_files.keys():
        file.write(f"{file_name}\n")



fasta_file.close()
ss_file.close()

mean_length = np.mean(protein_lengths)
std_length = np.std(protein_lengths)

with open(fasta_filename, 'a') as fasta_file:
    fasta_file.write(f"\nMean length: {mean_length}\n")
    fasta_file.write(f"Standard deviation of length: {std_length}\n")

