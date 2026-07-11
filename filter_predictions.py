import glob
import os

def filter_predictions(csv, threshold):
    '''Takes an oepn ML prediciton csv file (csv) and a confidence threshold (numeric value between 0-1)
    and returns a filtered list containing only ancient reads with confidence value above
    the user defined threshold.'''
    predictions = []
    for line in csv:
        if not line.startswith('Read_ID'):
            confidence = float(line.split(',')[1])
            if confidence >= threshold:
                predictions.append(line)
    return predictions

t
directory = '/path/to/predictions/'
threshold = 0.7

for file in glob.glob(f'{directory}*predictions.csv'):
    print(file)
    with open(f'{os.path.splitext(file)[0]}_{str(threshold)}_filtered.csv', 'w') as f_out:
        with open(file, 'r') as f_in:
            for item in filter_predictions(f_in, threshold):
                f_out.write(item)
