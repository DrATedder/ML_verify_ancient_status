import glob
import os

directory = '/storage02/or-microbio/meta_sim_test_25'

for sample in glob.glob(directory + '/aUoBsim*.fastq'):
    with open('{0}/{1}_RNT_{2}'.format(directory, os.path.basename(sample).split('_')[0]), os.path.basename(sample.split('_')[1])), 'w') as f_out:
        with open(sample, 'r') as f_in:
            for line in f_in:
                if line.startswith('@'):
                    f_out.write('{0}_a\n'.format(line.strip()))
                else:
                    f_out.write(line)
