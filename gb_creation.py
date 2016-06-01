#!/usr/bin/python

from Bio import SeqIO
import sys

strain_name = sys.argv[1]

# open the genbank file and take the interesting informations from the header to stock them in a 'header' list
with open("strains/{0}/sequence.gb".format(strain_name), 'r') as genbank_f:

    header = list()
    footer = list()
    for line in genbank_f:
        if line.find("COMMENT") != -1:
            break
        header.append(line)
    header.append("FEATURES             Location/Qualifiers\n")
# take the information from the footer to stock them in a 'footer' list. As the footer is at the end of the genbank file, the fact that the cursor is now at the line 'COMMENT' (and not at the beginning of the file) is not an issue.

    flag = False
    for line in genbank_f:
        if line.find("ORIGIN") != -1:
            flag = True
        if flag:
            footer.append(line)

# takes the information from the fna file to complete the 'header' list, especially the size of the DNA (which should be the same that indicated in the original genbank file, as both the files were taken from the NCBI database)
with open("strains/{0}/total.fna".format(strain_name), 'r') as fna_f:

    headline = fna_f.readline()
    organism_components = headline.split("|")[4].lstrip().split(" ")
    organism = " ".join(organism_components[0:3])
    strain = organism_components[2]
    nucleotid = ""
    for line in fna_f:
        nucleotid += line.strip()
    size = len(nucleotid)

source = '''     source          1..{_size_}
                     /organism="{_organism_}"
                     /mol_type="{_mol_type_}"
                     /strain="{_strain_}"\n'''.format(_size_=size, _organism_=organism, _mol_type_="genomic DNA", _strain_=strain)
header.append(source)

# print ''.join(header)

# takes the informations from the previously created fasta file and create a new fasta file with less informations to use as a input for psi-phi
output_fasta = open("./strains/{0}/input_psi_phi_{0}.faa".format(strain_name), 'w')
# takes the informations from the previously created fasta file and convert them into a genbank format to stock them in a 'section' list
with open("./strains/{0}/proteins_into_clusters_{0}".format(strain_name), 'r') as fasta_f:

    section = list()
    reader = SeqIO.parse(fasta_f, "fasta")
    for protein in reader:
        beginning = protein.description.split("|")[4].split(" # ")[1]
        end = protein.description.split("|")[4].split(" # ")[2]
        protein_sequence = str('/translation="' + protein.seq)
        protein_id = protein.description.split(" # ")[0].replace("|", "_")
        output_fasta.write(">{0}\n".format(protein_id))
        output_fasta.write("{0}\n".format(str(protein.seq)))

        string = ''
        for i in range(len(protein_sequence)):
            string = string + protein_sequence[i]
            if (i+1) % 59 == 0:
                string = string + '\n' + ' '*21

        gene_details = '''     gene            {_beginning_}..{_end_}
                     /locus_tag=""
                     /old_locus_tag=""
     CDS             {_beginning_}..{_end_}
                     /locus_tag=""
                     /old_locus_tag=""
                     /inference=""
                     /note=""
                     /codon_start=1
                     /transl_table=11
                     /product=""
                     /protein_id="{_protein_id_}"
                     /db_xref=""
                     {_protein_sequence_}"\n'''.format(_beginning_=beginning, _end_=end, _protein_id_=protein_id, _protein_sequence_=string)
        section.append(gene_details)

output_fasta.close()

# print ''.join(section)
# print ''.join(footer)

# output_fasta.close()

# open a new text file and transfer the informations into a genbank file
with open("strains/{0}/input_psi_phi_{0}.gb".format(strain_name), 'w') as f:
    f.write(''.join(header))
    f.write(''.join(section))
    f.write(''.join(footer))
