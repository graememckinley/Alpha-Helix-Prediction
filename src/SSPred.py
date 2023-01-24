# Python 3 script to edit for this project. 
# Note 1: Do not change the name of this file
# Note 2: Do not change the location of this file within the BIEN410_FinalProject package
# Note 3: This file can only read in "../input_file/input_file.txt" and "parameters.txt" as input
# Note 4: This file should write output to "../output_file/outfile.txt"
# Note 5: See example of a working SSPred.py file in ../scr_example folder

import numpy as np

infile = "../input_file/infile.txt"
parametersTxtFile = "parameters.txt"
predictionFile = "../output_file/outfile.txt"


# Read in FASTA file and return library
# This function is reused from src_example\SSPred
def readInput(file):
    inputData = {}

    with open(file, 'r') as f:
        while True:

            name = f.readline()
            seq = f.readline()

            if not seq:
                break

            inputData.update({name.rstrip(): seq.rstrip()})

    return inputData


# Write to the output file
# This function is reused from src_example\SSPred
def writeOutput(inputData, predictions, outfile):
    with open(outfile, 'w') as f:
        for name in inputData:
            f.write(name + "\n")
            f.write(inputData[name] + "\n")
            f.write(predictions[name] + "\n")


# Function for making predictions
def prediction(inputData, parametersFile):

    # Create a library to store predictions
    predictions = {}

    # Create libraries for storing probabilities
    # probH stores probabilities given there was a helix in the middle of the window during training
    # probSpace stores probabilities given there was a helix in the middle of the window during training
    probH = {}
    probSpace = {}

    # Read in the parameters text file
    # Parameters stores probabilities for each amino acid given helix and given space
    with open(parametersFile, 'r') as f:
        for line in f:

            # Alter each line in the parameters file in order to update dictionaries
            # The start of each line is the key to each of our dictionaries
            # The following numbers are our probabilities for each position in our window
            line = line.split(", ")
            aa = line[0]
            line = line[1:]
            line = [float(line[i]) for i in range(len(line))]

            # Define a window size we will use for making predictions
            windowSize = len(line)

            # The parameters file is split up into two sets of data
            # The first set corresponds to probabilities given helix
            # The second set corresponds to probabilities gives space
            if aa not in probH:
                probH.update({aa: line})
            else:
                probSpace.update({aa: line})

    # Define the prior log odds determined during training
    # It's easier to define the variable here than having a lone variable in our parameters file
    priorLogOdds = -0.4164735928449742

    # Iterate through the sequences of the input data and make predictions
    for name in inputData:

        # Define the sequence we will be making predictions on
        sequence = inputData[name]

        # Create a base prediction with no "H"
        # This gives us a list the same length as the sequence
        # We can then overwrite this list with our predictions
        pred = ["-" for aa in sequence]

        # Create a window that iterates through the sequence
        for i in range(len(sequence) - windowSize + 1):
            window = sequence[i: i + windowSize]

            # Calculate the value that allows us to make predictions
            # Iterate through the window dividing the probability given a helix for each position in
            # our window divided by the probability given a space
            predict = priorLogOdds + sum(
                [np.log(probH[window[j]][j] / probSpace[window[j]][j]) for j in range(len(window))])

            # Determine if whether to predict helix or space
            if predict > 0:
                pred[i + windowSize // 2] = "H"

        # Join the predictions into a string
        pred = "".join(pred)

        # Update the predictions library
        predictions.update({name: pred})

    return predictions


data = readInput(infile)

writeOutput(data, prediction(data, parametersTxtFile), predictionFile)
