import re
import argparse
import os

def main(args):

    with open(args.pipeline) as f:
        file = f.read()
    
    file = re.sub('{GENOME}', args.genome, file)
    file = re.sub('{DATA_PATH}', os.path.dirname(args.genome)+"/", file)
    file = re.sub('{SCRIPTS_PATH}', args.scripts_path.rstrip('/')+'/', file)
    file = re.sub('{PREFIX}', args.prefix, file)
    file = re.sub('{MIN_ORF_LENGTH}', args.min_orf, file)
    file = re.sub('{MIN_INT_LENGTH}', args.min_int, file)
    file = re.sub('{SEQ_IDENTITY}', args.seq_I, file)
    file = re.sub('{SEQ_LEN_CUTOFF}', args.seq_L, file)
    
    print(file)
    with open('Snakefile', 'w') as w:
        w.write(file)
    print('==========')
    print("Snakefile successfully generated...")
    print('==========')
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-pipeline', type=str, help='The path to the input genome.')
    parser.add_argument('-genome', type=str, help='The path to the input genome.')
    parser.add_argument('-scripts_path', type=str, help="The path to scripts dir.")
    parser.add_argument('-prefix', type=str, help='Additional prefix for your files with delimeter f.e. "hg19_".')
    parser.add_argument('-min_orf', type=str, default='300', help='Min ORF length.')
    parser.add_argument('-min_int', type=str, default='100', help='Min interval length.')
    parser.add_argument('-seq_I', type=str, default='0.75', help="This is the default cd-hit's 'global sequence identity' calculated as:\nnumber of identical amino acids in alignment\ndivided by the full length of the shorter sequence.")
    parser.add_argument('-seq_L', type=str, default='0.75', help="this is the default cd-hit's 'length identity cutoff',\nf.e., if set to 0.9, the shorter sequences need to be\nat least 0.9 length of the representative of the cluster.")
    
    args = parser.parse_args()
    
    # class Args():
    #     def __init__(self):
    #         self.pipeline = '../pipeline/Snakefile.template'
    #         self.genome = '../pipeline/data/Sclerotinia_sclerotiorum.ASM14694v1.dna.toplevel.fa'
    #         self.scripts_path = '../pipeline/scripts/'
    #         self.prefix = 'Scler_'
    #         self.min_orf = '300'
    #         self.min_int = '100'
    #         self.seq_I = '0.75'
    #         self.seq_L = '0.75'
    # args = Args()

    main(args)
    