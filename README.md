### **Humanin Protein Motif Identification Script**

#### Background Information

This is a repository containing scripts and resources for a computational tool that takes a FASTA
file containing the DNA sequence of the mitochondrial 16s gene, translates the sequence (using the EMBOSS Transeq tool) and identifies
the most likely _Humanin_ motif in each reading frame of the translated peptide sequence, using a probability 
matrix formed from identified _Humanin_ sequences in a variety of species (more details below).

#### Using the tool
Required libraries are contained in the requirements.txt file and can be installed as follows:

`pip install -r requirements.txt`

Navigate to the Project directory and run:
`python humanin.py -h`

This will produce a help message listing all the options for the tool:

_Usage: Humanin.py [options]_

_Options:_

_-h, --help            show this help message and exit_
  
_-i INPUT, --input=INPUT
                        Read data from FASTA file; can also be a directory
                        path if the -d option is used_
                        
_-o OUTPUT, --output=OUTPUT
                        Prefix of the output file_
                        
_-m MATRIX, --matrix=MATRIX
                        File containing the motif sequences used to create the
                        probability matrix; each line should hold a separate
                        sequence; all sequences should be the same length.
                        Default file is humanin.txt_
                        
_-d, --dir             Read data from FASTA files in this directory_
  
_-f FRAMES, --frames=FRAMES
                        Determines for which frames the most probable _Humanin_
                        sequence is returned; options: 1-6, All, or Best_
                        
_-p, --probs           Includes relative probability of the motif (as
                        compared to probability of consensus sequence)_

Example Usage:

`python humanin.py -i path_to_fasta -o Test\test`

This will output a file called "test_humanin_results.txt" to the Test directory that contains
 _Humanin_ motifs identified in the FASTA at the input path.

###### **Tool Options**

INPUT: Can be a path to a FASTA file or to a directory containing FASTA files. If the path is
a directory, all files in the directory will be processed and the results for each will be
combined into a single output file. The input file(s) should be a FASTA file containing the DNA sequence of the 16s gene.
This script does not support multi-organism FASTA files.

OUTPUT: Defines the prefix of the output file. The output file will be named as follows: prefix + '_humanin_results.txt'.

MATRIX: Path to the file used to create the probability matrix. The file should be a text file
containing known _Humanin_ peptide sequences - one on each line and 25 aa's long to account for potential stop codons. A default file is provided that contains _Humanin_ sequences
from the following species (in this order): 
* _Homo sapiens_
* _Pan troglodytes_
* _Gorilla gorilla_
* _Macaca mulatta_
* _Mus musculus_
* _Cavia porcellus_
* _Tupaia belangeri_
* _Felis catus_
* _Capra hircus_
* _Danio rerio_
* _Lepidotrigla microptera_
* _Dopasia gracilis_
* _Corvus corax_
* _Myxine glutinosa_
* _Ichthyomyzon unicuspis_
* _Neoceratodus forsteri_

These sequences were obtained from Figure 1 of the article: 

_Ian S. Logan. 2017. Pseudogenization of the _Humanin_ gene is common in the mitochondrial DNA of many vertebrates. Zoological Research, 38(4): 198-202_

dir: Indicates that the input is a directory containinf FASTA files

FRAMES: Indicates which frames will have their motif recorded in the output. A value of 'Best'
will include only the frames with the most likely _Humanin_ motif in the output; a value of 'All'
will include the most likely _Humanin_ motifs for each frame in the output; a number value (1-6) will include
the _Humanin_ motif for that specified frame in the output.

probs: Indicates that the output should include the relative probability for the _Humanin_ motif; this value is found by 
calculating the probability of the consensus sequence for the probability matrix (the max possible probability)

#### Translation

The translation of the 16s sequences was accomplished using the EMBOSS Transeq tool
(_Madeira F, Park YM, Lee J, et al. The EMBL-EBI search and sequence analysis tools APIs in 2019. Nucleic Acids Research. 2019 Jul;47(W1):W636-W641. DOI: 10.1093/nar/gkz268._)
via the REST API client script. Jobs are initialized using the following options:
* --frame 6 (all frames)
* --codontable 2 (Vertebrate Mitochondrial)


#### Probability Matrix

The probability is formed as follows:

Each sequence in the designated matrix file is lined up and for each position the count
unique amino acid value is calculated. This count is then used to generate a probability associated
with that particular amino acid at that particular position. These probabilities are stored in a 
dataframe which is then referenced for identifying the most probable _Humanin_ sequence within a peptide sequence.

#### Identifying _Humanin_

Once a sequence has been translated, each frame sequence is iterated through in 25 aa windows (step of 1) and the probability associated with
each amino acid at that particular position is obtained from the probability matrix; if the amino
acid has no associated probability at that that position (meaning none of the matrix sequences had that aa at that position)
the amino acid is given the minimum possible non-zero probability. These probabilities are multiplied
together to identify the probability of the subsequence being a _Humanin_ motif; the subsequence with the highest
probability is returned as the 25 aa-long _Humanin_ motif for that frame.


