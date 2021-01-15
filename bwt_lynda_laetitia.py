# Projet realisé par:
    # MESSAD DIT MAHTAL Lynda & HOLLIER Laetitia
    # Formation M2BI univerité de Paris

import argparse
import os
import sys
import copy
from collections import Counter



def get_arguments():
    """
    This function parses the command line arguments

    Arguments: python3 command-line

    Returns : fasta_file : string
    """

    #Declaration of expexted arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="Input fasta sequences file.", type=str, required=True, dest = 'seq_file')
    parser.add_argument("-r", "--reads", help="Input reads file.", type=str, required=True, dest = 'reads_file')
    parser.add_argument("-o", "--output", help="CSV output file.", type=str, default = "results.csv", dest = 'results_file')


    return parser.parse_args()

def parse_fasta(fasta_file):
    """
    This function reads the input fasta file to get the sequence
    
    Arguments: input fasta file  : str

    Returns : sequence :str
    """

    sequence = ""

    with open(fasta_file, "r") as f_in:
        for line in f_in:
            if line[0] == ">":
                continue
            else:
                sequence += line.strip()
    return sequence 

def parse_reads(reads_file):
    """
    This function reads the reads file in fastq format

    Arguments:
        FASTQ file : str

    Return:
        reads : in dico format , key = read's id ; value = sequence

    """

    reads_dict = {}

    with open(reads_file, "r") as f_in:
        for line in f_in:
            if line.startswith(">"):
                read_id = line.split(".")[0].split(">")[1]
                reads_dict[read_id] = ""
            else:
                reads_dict[read_id] += line.strip()

    return reads_dict

def burrows_wheeler(seq):
    """
    This function follows the Borrows Wheeler transform method.

    Parameters:
        seq: str / referenced sequence from the fasta file

    Return:
        bwt : str / B.W transformed sequence
    """
    permuts = ['$' + seq]
    for i in range(len(seq)):
        last_nucl= permuts[i][-1]
        element = last_nucl + permuts[i][0:-1]
        permuts.append(element)

    # sored result
    permuts.sort() 

    # get the last column
    bwt = ""
    for permut in permuts:
        bwt += permut[-1]

    return bwt

def bwt_reverse(bwt):
    """
    This function generates the inverse of the transformed sequence

    Parameters:
        bwt: str / the transformed sequence

    Return:
        seq : str / inverse of B.W transformed
    """
    reverse = [""] * len(bwt)

    for i in range(len(bwt)):
        elmt = [bwt[i] + reverse[i] for i in range(len(bwt))]

        reverse = sorted(elmt)

    seq = [row for row in reverse if row.startswith('$')][0]
    return seq.strip('$')

def tally_table(bwt):
    """
    Generats a tally table from the transformed b.w sequence
    
    Parameter
    ---------
    bwt: str / B.W transformed sequence
    
    Return
    ------
    tally: list/ contains dictionary of the occurences
    """

    dico= {'$': 0, 'A': 0, 'T': 0, 'C': 0, 'G': 0}
    tally =  []
    for nucl in bwt:
        dico[nucl] += 1
        tally.append(copy.deepcopy(dico))

    return tally

def count_number(bwt):
    """
    Generats the number of character occurrences
    
    Parameter
    ---------
    bwt: str / the transformed sequence 
    
    Return
    ------
    inferior: dict / contains dicitonary of the occurences
    """
    inferior, count = {'$': 0, 'A': 0, 'C': 0, 'G': 0, 'T': 0}, Counter(sorted(bwt))
    nb = 0
    for nucl in inferior.keys():
        inferior[nucl] = nb
        nb += count[nucl]
    return inferior



def fill(text, width=80):
    """
    Split text with a line return to respect fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))



def write_transform(bwt):
    """
    Whrite the BW transformed sequence in a text file

    Parameter
    ---------
    bwt: str / the transformed sequence 
    """
    with open("BW_transformed_sequence.txt", "w") as filout:
        filout.write("> Burrows-Wheeler Transform" + "\n")
        filout.write(fill(bwt))



if __name__ == '__main__':
    # récupérer les arguments
    args = get_arguments()
    
    # les séquences du fasta file
    sequences = parse_fasta(args.seq_file)
    print("la sequence trouvée(s) : {} \n".format(sequences))
    # les reads 
    reads = parse_reads(args.reads_file)
    print("les reads trouvés : {} \n".format(reads))
    # construction de la bwt
    bwt  = burrows_wheeler(sequences)
    print("la transformée de B.W : {} \n".format(bwt))

    rev =  bwt_reverse(bwt)
    print("la reverse de B.W : {} \n".format(rev))

    tally_nb = tally_table(bwt)
    print("tally table : {} \n".format(tally_nb))

    number = count_number(bwt)
    print("numeros inferieurs: {} \n".format(number))
