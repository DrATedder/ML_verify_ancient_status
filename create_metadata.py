import os
import glob

directory = '/storage02/or-microbio/meta_sim_test_25/ML_data_sets'

data_dict = {}
last_set = ''

for r1 in glob.glob(directory + '/dataset_*_R1.fasta'):
    sample = os.path.basename(r1).strip('_R1.fasta')
    last_set = sample
    data_dict[sample] = []
    with open(r1, 'r') as f_in:
        for line in f_in:
            if line.startswith('>'):
                if line.endswith('_a\n'):
                    data_dict[sample].append('1')
                else:
                    data_dict[sample].append('0')
        with open('{0}/{1}_metadata.csv'.format(directory,sample), 'w') as f_out:
            data_length = len(data_dict[last_set])
            f_out.write('dataset,{0}\n'.format(','.join([str(x) for x in range(0,data_length,1)])))
            for k,v in data_dict.items():
                f_out.write('{0},{1}\n'.format(k, ','.join(v)))
    
