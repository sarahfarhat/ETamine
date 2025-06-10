#!/usr/bin/python3
import argparse

def parse_fasta(fasta_file):
    """Crée un dictionnaire des séquences FASTA."""
    dic = {}
    with open(fasta_file, "r") as fil:
        for line in fil.readlines():
            if line.startswith(">"):
                fid = line.strip()
                dic.setdefault(fid, [])
            else:
                dic[fid].extend(char for char in line.strip())
    return dic

def extract_seqname(dic):
    """Extrait le nom basé sur les 4 premiers champs des séquences."""
    # On récupère le premier identifiant dans le dictionnaire
    first_id = next(iter(dic.keys()))
    # Extraction des 4 premiers champs séparés par "_"
    fields = first_id.lstrip(">").split("_")
    seqname = "_".join(fields[:4])  # Combine les 4 premiers champs
    return seqname

def calculate_consensus(dic, with_insertions):
    """Calcule le consensus en fonction de l'option avec ou sans insertions."""
    tail_seq = len(next(iter(dic.values()), []))
    consensus = ""

    for i in range(tail_seq):
        tokeep = {}
        for el in dic:
            char = dic[el][i]
            tokeep[char] = tokeep.get(char, 0) + 1

        if not with_insertions and "-" in tokeep:
            del tokeep["-"]

        if len(tokeep) > 1:
            sorted_tokeep = sorted(tokeep.items(), key=lambda x: x[1], reverse=True)
            toprint = sorted_tokeep[0][0]
            if len(sorted_tokeep) > 1 and sorted_tokeep[0][1] == sorted_tokeep[1][1]:
                toprint = sorted_tokeep[1][0] if sorted_tokeep[0][0] == "-" else sorted_tokeep[0][0]
        else:
            toprint = list(tokeep.keys())[0]

        consensus += toprint
    return consensus

def main():
    parser = argparse.ArgumentParser(description='Consensus sequence generator (with or without insertions).')
    parser.add_argument('-f', '--fasta', dest='fastaFile', type=str, help='Input FASTA file of a multiple alignment', required=True)
    parser.add_argument('-o', '--output', dest='output_file', type=str, default="consensus.fa", help="Output file")
    args = parser.parse_args()

    dic = parse_fasta(args.fastaFile)
    num_sequences = len(dic)

    # Choisir le mode en fonction du nombre de séquences
    with_insertions = num_sequences > 2

    consensus = calculate_consensus(dic, with_insertions)

    # Extrait les 4 premiers champs pour le nom de la séquence
    seqname = extract_seqname(dic)

    with open(args.output_file, "w") as f:
        f.write(f">{seqname}\n")
        f.write(consensus + "\n")

if __name__ == "__main__":
    main()
