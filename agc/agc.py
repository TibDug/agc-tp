#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Thibault Dugauquier"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Thibault Dugauquier"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Thibault Dugauquier"
__email__ = "thibault.dug@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """prend deux arguments correspondant au fichier fasta ou fasta.gz et à la longueur
    minimale des séquences et retourne un générateur de séquences de longueur l >= minseqlen:
    yield sequence"""
    if amplicon_file.endswith(".gz"):
        file_in = gzip.open(amplicon_file)
        seq = b""
        for ligne in file_in:
            if ligne.startswith(b">"):
                if len(seq) >= minseqlen:
                    yield seq.decode('ascii')
                seq = b""
            else:
                seq += ligne[:-1]
        yield seq.decode('ascii')
    else:
        file_in = open(amplicon_file)
        seq = ""
        for ligne in file_in:
            if ligne.startswith(">"):
                if len(seq) >= minseqlen:
                    yield seq
                seq = ""
            else:
                seq += ligne[:-1]
        yield seq
    file_in.close()

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """Prend trois arguments correspondant au fichier fasta,  la longueur minimale des
    séquences et leur comptage minimum. Elle fait appel au générateur fourni par read_fasta
    et retourne un générateur des séquences uniques ayant une occurrence O>=mincount ainsi
    que leur occurrence. Les séquences seront retournées par ordre décroissant d’occurrence:
    yield [sequence, count]"""
    seq_dict = {}
    for seq in list(read_fasta(amplicon_file, minseqlen)):
        print(seq)
        if seq in seq_dict:
            seq_dict[seq] += 1
        else:
            seq_dict[seq] = 1
    for seq, seq_count in sorted(seq_dict.items(), key = lambda item: item[1], reverse = True):
        if seq_count >= mincount:
            yield [seq, seq_count]

def get_chunks(sequence, chunk_size):
    """prend une séquence et un longueur de segment l: chunk_size et retourne une liste de
    sous-séquences de taille l non chevauchantes"""
    chunk_list = []
    seq_position = 0
    while True:
        if seq_position+chunk_size >= len(sequence):
            break
        chunk_list.append(sequence[seq_position:seq_position+chunk_size])
        seq_position += chunk_size
    return chunk_list

def cut_kmer(sequence, kmer_size):
    """prend une séquence, une taille de k-mer et retourne un générateur de k-mer"""
    for nuc in range(len(sequence[:-kmer_size+1])):
        yield sequence[nuc:nuc+kmer_size]

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """prend un dictionnaire ayant pour clé un index de kmer et pour valeur une liste
    d’identifiant des séquences dont ils proviennent."""
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer in kmer_dict:
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer] = [id_seq]
    return kmer_dict

def search_mates(kmer_dict, sequence, kmer_size):
    """prend un dictionnaire ayant pour clé un index de kmer et pour valeur une liste
    d’identifiant des séquences dont ils proviennent, une séquence et une longueur de kmer:
    kmer_size."""
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size) if \
    kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]

def get_identity(alignment_list):
    """prend un alignement (sous forme de liste) et calcule le pourcentage d’identité entre
    les deux séquences"""
    identical_count = 0
    for indice_nuc in range(len(alignment_list[0])):
        if alignment_list[0][indice_nuc] == alignment_list[1][indice_nuc]:
            identical_count += 1
    return identical_count * 100 / len(alignment_list[0])

def detect_chimera(perc_identity_matrix):
    """prend une matrice donnant par segment le taux d’identité entre la séquence candidate
    et deux séquences parentes et retourne un booléen indiquant si la séquence candidate est
    une chimère (True) ou ne l’est pas (False)"""
    pass

def get_unique(ids):
    """get_unique"""
    return {}.fromkeys(ids).keys()

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Fait appel au générateur fourni par dereplication_fulllength et retourne un
    générateur des séquences non chimérique au format: yield [sequence, count]"""
    pass

def common(lst1, lst2):
    """retourne les éléments communs entre deux listes"""
    return list(set(lst1) & set(lst2))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    main()
