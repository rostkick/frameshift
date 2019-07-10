import pandas as pd
import numpy as np
import argparse


def main(args):
    
    df = pd.read_csv(args.input, sep='\t')
    df['id'] = np.arange(len(df))
    df['id'] = 'fs_id' + df['id'].astype(str)
    cols = list(df.columns)
    cols = [cols[-1]] + cols[:-1]
    df = df[cols]
    
    df.to_csv(args.output_full_tsv, sep='\t', header=None, index=None)

    def gather( df, key, value, cols ):
        id_vars = [ col for col in df.columns if col not in cols ]
        id_values = cols
        var_name = key
        value_name = value
        return pd.melt( df, id_vars, id_values, var_name, value_name )

    subset = df.loc[:, ['id', 'PROT_INTERSECTED_MAIN_SEQ', 'PROT_INTERSECTED_ALT_SEQ']]
    subset = gather(subset, 'type', 'protein', ['PROT_INTERSECTED_MAIN_SEQ', 'PROT_INTERSECTED_ALT_SEQ'])
    subset = subset.sort_values(by='id', ascending=True)
    subset['id'] = subset['id'].astype(str) + ':' + subset['type']
    subset = subset.drop('type', axis=1)
    
    subset.to_csv(args.output_subset_tsv, sep="\t", header=None, index=None)

    with open(args.output_subset_fasta, 'w') as w:
        for i, j in zip(subset['id'], subset['protein']):
            w.write(f'>{i}\n{j}\n')
            
            
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-input', type=str)
    parser.add_argument('-output_full_tsv', type=str)
    parser.add_argument('-output_subset_fasta', type=str) 
    parser.add_argument('-output_subset_tsv', type=str) 

    args = parser.parse_args()
    
    main(args)