# Alpha-Helix-Prediction
Bayesian probabilistic model for predicting alpha helix secondary structure

Assignment completed in fulfillment of BIEN 410 - Computational Methods in Bioengineering at McGill University. 
File directory and names were provided by the Instructor and had to be closely adhered to.


## File Breakdown

**input_file/infile.txt** - *text file containing ~5000 amino acid sequences in FASTA format*

**output_file/outfile.txt** - *text file containing the alpha helix predictions for the sequences in infile.txt*

**src/parameters.txt** - *text file containing probabilities generated during the training of the model and used when making predictions*

**src/SSPred.py** - *script for making predictions*

**src_example/** - *example files provided by the instructor*

**training_data/labels.txt** - *text file containing the same amino acid sequences as infile.txt as well as the correct alpha helix sequence*

**testing/testing.py** - *script for comparing the alpha helix sequence predicted by SSPred.py to the correct sequence in labels.txt*

**praccy.py** - *code originally used for training the model and generating the probabilities in parameters.txt, this was not intended to be run by itself and was originally a part of SSPred.py*


## Run Explanation

Running SSPred.py will make alpha helix predictions based on the contents of infile.txt, writing predictions to outfile.txt.
The accuracy of the predictions made can the be checked by running testing.py.
Changes to the probabilites used when making predictions can be made by editing parameters.txt but please refer to the document below for how the file is formatted.
Amino acid sequences can be edited or added in infile.txt, however, identical changes will have to be made to labels.txt.
Since testing.py compares outfile.txt to labels.txt and outfile.txt is based on infile.txt, not making these changes will result in an incorrect accuracy being calculated by testing.py.

## Assignment Writeup

Below is a write up done for the assignment explaining the approach, discussing advantages/disadvantages, and briefly exploring future improvements.

[Git Writeup.pdf](https://github.com/graememckinley/Alpha-Helix-Prediction/files/10491780/Git.Writeup.pdf)
