# edit_protein_mutation
This script is designed to analyze the output from CRISPResso2, focusing specifically on a designated Open Reading Frame (ORF). It systematically identifies the types of amino acid mutations within the specified ORF and calculates the proportions of frameshift mutations and premature stop codons occurring. This analysis helps in understanding the impact of CRISPR-Cas9 editing on the protein function by quantifying disruptive mutations that could potentially alter the protein structure and function.

## Usage:  
`python edit_protein_mutation_public.py [-h] [--reverse] file ORF_start ORF_end cut_off`  
### arguments:
  1. Input file: Provide the path to the CRISPResso2 output file
  2. ORF_start : Start position of the ORF (Open Reading Frame). Must be an integer.
  3. ORF_end : End position of the ORF. Must be an integer.
  4. cut_off : Threshold or cut-off value for processing. Must be a float.
  5. reverse : Process in reverse sequence.
### Command Line Example:
  `python edit_protein_mutation_public.py --reverse edit_protein_mutation_input.txt 15 21 0.01 `
### output
  HS: 69.20446827943142  
  HL: 30.378289345292462  
  YL: 0.019125235162330313  
  ?L: 0.010928705807045894  (? means stop coden)  
  frameshift_mutation : 0.038250470324660626  
  early_stopping : 0.010928705807045894  
  both : 0  
