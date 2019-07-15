import re
import pandas as pd
from Bio import SeqIO
import argparse


def main(args):
    
    pdb_id, pdb_seq = [], []
    for record in SeqIO.parse(args.input_pdb_fasta, 'fasta'):
        pdb_id.append(record.id)
        pdb_seq.append(str(record.seq))

    pdb_df = pd.DataFrame({'pdb_id':pdb_id, 'pdb_seq':pdb_seq})

    blastp_df = pd.read_csv(args.input_blastp_tsv, sep='\t', header=None)
    blastp_df[1] = blastp_df[1].str.replace('.*pdb', 'pdb')
    blastp_df.columns = ['query', 'pdb_id', 'evalue', 'query_algn', 'match_algn']

    matched_df = pd.merge(blastp_df, pdb_df, on='pdb_id')

    df_intersected = pd.read_csv(args.input_intersected_tsv, sep='\t', header=None)
    df_intersected.columns = ['query', 'query_seq']

    intersected_subset_data = pd.merge(df_intersected, matched_df, on='query', how='outer')

    intersected_subset_data['id'] = intersected_subset_data['query'].str.split(':', expand=True)[0]
    intersected_subset_data['ALT/MAIN'] = intersected_subset_data['query'].str.split(':', expand=True)[1]

    intersected_subset_data.drop('query', axis=1, inplace=True)

    intersected_subset_data = intersected_subset_data[['id', 'ALT/MAIN', 'query_seq', 'pdb_id', 'evalue', 'query_algn', 'match_algn', 'pdb_seq']]
    intersected_subset_data = intersected_subset_data.groupby(['id', 'ALT/MAIN']).first()

    intersected_subset_data_omit_NA = intersected_subset_data.groupby('id').filter(lambda g: g.isnull().sum().sum()==0)



    intersected_subset_data_omit_NA.to_csv(args.output_omit_NA, sep='\t')

    intersected_subset_data.fillna('.', inplace=True)
    intersected_subset_data.to_csv(args.output, sep='\t')
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-input_pdb_fasta', type=str)
    parser.add_argument('-input_blastp_tsv', type=str)
    parser.add_argument('-input_intersected_tsv', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-output_omit_NA', type=str)


    args = parser.parse_args()
    
    main(args)