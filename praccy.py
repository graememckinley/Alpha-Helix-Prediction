import numpy as np

"""
Original code for training the model and generating probabilities
This code was not intended to be run by itself, it was originally in SSPred.py
"""


def prediction(inputData, parametersFile):
    # Create a library to store predictions
    predictions = {}

    # Open labels file
    labels = {}

    with open("../training_data/labels.txt", 'r') as f:
        while True:

            name = f.readline()
            seq = f.readline()
            hel = f.readline()

            if not hel:
                break

            labels.update({name.rstrip(): (seq.rstrip(), hel.rstrip())})

    # Define all the abbreviations for amino acids
    abbrev = ["A", "R", "D", "N", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    # Define a window size
    windowSize = 9

    # Define and create dictionaries
    probH = {aa: np.zeros(windowSize) for aa in abbrev}
    probSpace = {aa: np.zeros(windowSize) for aa in abbrev}
    totalResidues = 0
    totalHelix = 0
    totalSpace = 0

    # Training
    for name in labels:
        sequence = labels[name][0]
        helixSequence = labels[name][1]

        for i in range(len(sequence)):
            totalResidues += 1

            # Count the total number of Helix and Space
            if helixSequence[i] == "H":
                totalHelix += 1
            else:
                totalSpace += 1

        # Create a window and update the dictionary of conditional probability
        for i in range(len(sequence) - windowSize + 1):

            # Create windows in each sequence
            windowSequence = sequence[i: i + windowSize]
            windowHelix = helixSequence[i: i + windowSize]

            # Determine if there is an H at middleAmino
            if windowHelix[windowSize // 2] == "H":

                # Iterate through the sequence window
                for j in range(windowSize):
                    probH[windowSequence[j]][j] += 1

            else:

                # Iterate through the sequence window
                for j in range(windowSize):
                    probSpace[windowSequence[j]][j] += 1

    # Calculate probabilities
    for aa in probH:
        for i in range(windowSize):
            probH[aa][i] /= totalHelix
            probSpace[aa][i] /= totalSpace

    # Print
    for aa in probH:
        print(aa, [probH[aa][i] for i in range(windowSize)])

    print("")

    for aa in probSpace:
        print(aa, [probSpace[aa][i] for i in range(windowSize)])

    for aa in probH:
        print(aa + ":", [int(probH[aa][i] * 100) for i in range(windowSize)])

    print("")

    for aa in probSpace:
        print(aa + ":", [int(probSpace[aa][i] * 100) for i in range(windowSize)])

    # Iterate through the sequences and make predictions
    for name in inputData:
        sequence = inputData[name]

        # Create a base prediction with no "H"
        pred = ["-" for aa in sequence]

        # Create a window that iterates through the sequence
        for i in range(len(sequence) - windowSize + 1):
            window = sequence[i: i + windowSize]

            predict = sum([np.log(probH[window[j]][j] / probSpace[window[j]][j]) for j in range(len(window))])

            predict += np.log((totalHelix / totalResidues) / (totalSpace / totalResidues))

            if predict > 0:
                pred[i + windowSize // 2] = "H"

        pred = "".join(pred)
        predictions.update({name: pred})

    return predictions

