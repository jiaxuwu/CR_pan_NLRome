module load miniconda/3-py3.10 agat/1.2.0
source activate agat-1.2.0

agat_sp_extract_sequences.pl \
        --cfs \
        -g westar_helixer.gff \
        -f Brassica_napus.Westar.v0.genome.fa \
        -p -o westar_peptide.fasta


from Bio import SeqIO

input_file = "H1_peptide.fasta"
output_file = "H1_peptide_final.fasta"

# remove any description from the fasta file
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        record.description = ""
        SeqIO.write(record, outfile, "fa")


def modify_fasta(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                # Split the header line into parts
                parts = line.strip().split('_')
                # Modify the parts to remove "_helixer" and "_polished"
                new_parts = [part for part in parts if part not in ['helixer', 'polished']]

                # Reconstruct the modified header line
                new_header = '_'.join(new_parts)

                # Write the modified header line to the output file
                f_out.write(new_header + '\n')
            else:
                # Write the sequence line to the output file as is
                f_out.write(line)

if __name__ == "__main__":
    input_file = "input.fasta"
    output_file = "output.fasta"
    modify_fasta(input_file, output_file)
    print("Sequence names modified. Output written to", output_file)        
        
