import pickle

def readCases(file):
    with open(file, 'rb') as f:
        cases = pickle.load(f)
    return cases