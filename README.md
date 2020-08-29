# get_ORF_aa

 This function takes three files as input as follows:
1. sequences is the fasta format of DNA sequences of one or more chromosom fragments
 2. intervals is a table containing four columns named chr, start, stop, and id. the intervals of the 
 genes to be translated are defined in this file.
 3. codons is a table containing the genetic codons for RNA translation to amino acids. 

 if no codons file provided the function will go through the R built-in translation function
 and finish the task. 

 The output of the function will be file of standard fasta format of translated amino acids
 for defined intervals, which will be saved in working directory.
 The function will also print the results on the console.

 The required packages are defined in the function, which need to be installed before running the function.
