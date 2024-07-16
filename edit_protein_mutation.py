import os
import re
import sys
from collections import defaultdict
import argparse
def DNA_reverse_complement2(sequence):
    """Reverse complement of a DNA sequence.

    Args:
        sequence (str): The DNA sequence to be reverse complemented.

    Returns:
        str: The reverse complement of the input DNA sequence.
    """
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    return sequence.translate(trantab)[::-1]

coden2aa_dict = {'TTT': "F", 'TTC': "F", 'TTA': "L", 'TTG': "L", 'CTA': "L", 'CTT': "L", 'CTC': "L", 'CTG': "L",
                 'ATA': "I", 'ATT': "I", 'ATC': "I", 'ATG': "M", 'GTT': "V", 'GTA': "V", 'GTG': "V", 'GTC': "V",
                 'TCT': "S", 'TCG': "S", 'TCA': "S", 'TCC': "S", 'AGC': "S", 'AGT': "S", 'CCT': "P", 'CCG': "P",
                 'CCA': "P", 'CCC': "P", 'ACT': "T", 'ACG': "T", 'ACA': "T", 'ACC': "T", 'GCA': "A", 'GCC': "A",
                 'GCT': "A", 'GCG': "A", 'GAA': "E", 'GAG': "E", 'GAC': "D", 'GAT': "D", 'AAA': "K", 'AAG': "K",
                 'AAC': "N", 'AAT': "N", 'CAA': "Q", 'CAG': "Q", 'CAT': "H", 'CAC': "H", 'TAA': "?", 'TAG': "?",
                 'TGA': "?", 'TAT': "Y", 'TAC': "Y", 'TGC': "C", 'TGT': "C", 'TGG': "W", 'CGA': "R", 'CGC': "R",
                 'CGT': "R", 'CGG': "R", 'AGA': "R", 'AGG': "R", 'GGA': "G", 'GGT': "G", 'GGC': "G", 'GGG': "G"}


def get_percentage(file,ORF_start=3,ORF_end=12,cut_off=0.01,reverse=False):
    """
    Args:
        file (_type_): Path to the file to be processed.
        ORF_start (int, optional): Start position of the ORF (Open Reading Frame). Must be an integer. Defaults to 3.
        ORF_end (int, optional): End position of the ORF. Must be an integer. Defaults to 12.
        cut_off (float, optional): Threshold or cut-off value for processing. Must be a float. Defaults to 0.01.
        reverse (bool, optional): Process in reverse sequence. Defaults to False.
    Returns:
        [protein_dict:dict, frameshift_mutation:float, early_stopping:float, both:float]
    """
    data = open(file,"r")
    if (ORF_end-ORF_start)%3 !=0:
        print("length of ORF is ERROR")
        sys.exit(1)
    first_line=0
    frameshift_mutation = 0
    early_stopping = 0
    both = 0
    protein_dict=defaultdict(lambda:0)
    for line in data:
        frameshift_mutation_flag=False
        if "Aligned_Sequence" in line:
            continue
        target,single_ref,edit_or_not,n_deleted,n_inserted,n_mutated,Reads,Reads_per=line.strip().split("\t")
        if  (float(Reads_per) >= cut_off) and (int(n_deleted)+int(n_inserted))%3!=0:
            frameshift_mutation += float(Reads_per)
            frameshift_mutation_flag=True
        if (int(n_deleted)+int(n_inserted))%3==0 and float(Reads_per) >= cut_off:
            protein_rna_Seq=target.replace("-","")[ORF_start-1:ORF_end-1].upper()
            if reverse:
                protein_rna_Seq=DNA_reverse_complement2(protein_rna_Seq)
                
            protein_seq = "".join(
            coden2aa_dict[protein_rna_Seq[_:_ + 3]]
            for _ in range(0, len(protein_rna_Seq)-len(protein_rna_Seq)%3, 3))
            

            if frameshift_mutation_flag:
                if "?" in protein_seq:
                    both+=float(Reads_per)
                    frameshift_mutation -= float(Reads_per)

            else:
                if "?" in protein_seq:
                    early_stopping+=float(Reads_per)

            protein_dict[protein_seq]+=float(Reads_per)
    return protein_dict,frameshift_mutation,early_stopping,both



def parse_arguments():
    parser = argparse.ArgumentParser(description="Process a file with given parameters.")
    parser.add_argument('file', type=str, help='Path to the file to be processed.')
    parser.add_argument('ORF_start', type=int, help='Start position of the ORF (Open Reading Frame). Must be an integer .')
    parser.add_argument('ORF_end', type=int, help='End position of the ORF. Must be an integer.')
    parser.add_argument('cut_off', type=float,default=0.01 ,help='Threshold or cut-off value for processing. Default is 25.0')
    parser.add_argument('--reverse', action='store_true',default=False ,help='Process in reverse sequence. Default is False')
    return parser.parse_args()

def main():
    args = parse_arguments()
    protein_dict,frameshift_mutation,early_stopping,both = get_percentage(args.file, args.ORF_start, args.ORF_end, args.cut_off, args.reverse)
    for key, value in protein_dict.items():
        print(f"{key}: {value}")
    print(f"frameshift_mutation : {frameshift_mutation}\nearly_stopping : {early_stopping}\nboth : {both}")

        
if __name__ == "__main__":
    main()
