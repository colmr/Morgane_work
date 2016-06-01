#!/usr/bin/python

import csv
import sys

strain = sys.argv[1]

# Take all the positions, in the nucleotidic sequence, of the lantibiotic clusters found in the bacteria strain
list_positions_clusters = []
with open("cluster_newScript_C1_50kb_3rdscore_ggsearch_newinfo.csv", "r") as f:
    reader = csv.DictReader(f, delimiter=",")
    for line in reader:
        if line["Species"] == strain:
            list_positions_clusters.append((int(line["Start"]), int(line["End"])))

print list_positions_clusters


# open the fasta file and store it in a dictionary
with open("./strains/{0}/translate_prod.txt".format(strain), "r") as f:
    in_silico_proteome = {}
    for line in f:
        if line[0] == ">":
            description = line
            in_silico_proteome[description] = ""
        else:
            in_silico_proteome[description] += line.strip()


# fonction which looks if a protein is in the cluster or not
def in_cluster(details):
    position = ( int(details.split("|")[4].split(" # ")[1]), int(details.split("|")[4].split(" # ")[2]) )
    for cluster_position in list_positions_clusters:
        if position[0] >= cluster_position[0] and position[1] <= cluster_position[1]:
            return True
    return False

# write a new fasta file with only defined proteins that are present into clusters
with open("./strains/{0}/proteins_into_clusters_{0}".format(strain), "w") as f:
    for description in in_silico_proteome.keys():
        if in_cluster(description) is True:
            f.write("{0}{1}\n\n".format(description, in_silico_proteome[description]))
