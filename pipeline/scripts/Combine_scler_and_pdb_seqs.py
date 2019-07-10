import os
import re
import argparse
from collections import defaultdict
from Bio import SeqIO
import tempfile


def parse_data(args):

    pdb_fasta = dict()
    translate_dict = defaultdict(list)

    for record in SeqIO.parse(args.pdb_fasta, 'fasta'):
        pdb_fasta[record.id] = str(record.seq)

    with open(args.pdb_tsv) as f:
        for line in f:
            id_pred = line.split()[0]
            id_pdb = '|'.join(re.split('\|', (line.split())[1])[2:5])
            translate_dict[id_pred].append(id_pdb)

    with open(args.predicted_fasta) as f:
        predicted_fasta = f.read()
    
    return pdb_fasta, translate_dict, predicted_fasta
        
        
class Batch:
    
    
    def __init__(self, batch):
        
        self.batch = batch
        
        self.global_orf_name, self.global_orf = '', ''
        self.dna_main_name, self.dna_main = '', ''
        self.dna_alt_name, self.dna_alt = '', ''
        self.prot_main_name, self.prot_main = '', ''
        self.prot_alt_name, self.prot_alt = '', ''
        self.inter_dna_main_name, self.inter_dna_main = '', ''
        self.inter_dna_alt_name, self.inter_dna_alt = '', ''
        self.prot_main_intersected_name, self.prot_main_intersected = '', ''
        self.prot_alt_intersected_name, self.prot_alt_intersected = '', ''
        
        self.pdb_main = ''
        self.pdb_alt = ''
        
        self.full_batch = ''
        
    def get_items(self):
        
        batch_list = iter(self.batch.split('\n'))
        
        for line in batch_list:
            if '>GLOBAL_ORF' in line:
                self.global_orf_name = line.lstrip('>')
                self.global_orf = next(batch_list)
                
            elif '>DNA_MAIN' in line:
                self.dna_main_name = line.lstrip('>')
                self.dna_main = next(batch_list)
                
            elif '>DNA_ALT' in line:
                self.dna_alt_name = line.lstrip('>')
                self.dna_alt = next(batch_list)
                
            elif '>INTERSECTION_DNA_MAIN' in line:
                self.inter_dna_main_name = line.lstrip('>')
                self.inter_dna_main = next(batch_list)
            
            elif '>INTERSECTION_DNA_ALT' in line:
                self.inter_dna_alt_name = line.lstrip('>')
                self.inter_dna_alt = next(batch_list)
                
            elif '>PROT_MAIN:' in line:
                self.prot_main_name = line.lstrip('>')
                self.prot_main = next(batch_list)
                
            elif '>PROT_ALT:' in line:
                self.prot_alt_name = line.lstrip('>')
                self.prot_alt = next(batch_list)
            elif '>PROT_INTERSECTED_MAIN' in line:
                self.prot_main_intersected_name = line.lstrip('>')
                self.prot_main_intersected = next(batch_list)
                
            elif '>PROT_INTERSECTED_ALT' in line:
                self.prot_alt_intersected_name = line.lstrip('>')
                self.prot_alt_intersected = next(batch_list)

    def get_pdb(self, pdb_dict, translate_dict):

        main_pdb_seqs, alt_pdb_seqs = [], []
        
        main_pdb_names = translate_dict[self.prot_main_intersected_name]
        main_pdb_seqs = [pdb_dict[k] for k in main_pdb_names]
        
        alt_pdb_names = translate_dict[self.prot_alt_intersected_name]
        alt_pdb_seqs = [pdb_dict[k] for k in alt_pdb_names]
        
        main_pdb_dict = dict(zip(main_pdb_names, main_pdb_seqs))
        alt_pdb_dict = dict(zip(alt_pdb_names, alt_pdb_seqs))

        if len(main_pdb_dict.keys())>0:
            for k, v in main_pdb_dict.items():
                self.pdb_main += f'>{k}\n{v}\n'
        else:
            self.pdb_main = '>NOT_INFO\n...\n'
        
        if len(alt_pdb_dict.keys())>0:
            for k, v in alt_pdb_dict.items():
                self.pdb_alt += f'>{k}\n{v}\n'
        else:
            self.pdb_alt = '>NOT_INFO\n...\n'
            
    def combine_batch(self):
        self.full_batch = f'>{self.global_orf_name}\n'\
                            f'{self.global_orf}\n'\
                            f'>{self.dna_main_name}\n'\
                            f'{self.dna_main}\n'\
                            f'>{self.dna_alt_name}\n'\
                            f'{self.dna_alt}\n'\
                            f'>{self.inter_dna_main_name}\n'\
                            f'{self.inter_dna_main}\n'\
                            f'>{self.inter_dna_alt_name}\n'\
                            f'{self.inter_dna_alt}\n'\
                            f'>{self.prot_main_name}\n'\
                            f'{self.prot_main}\n'\
                            f'>{self.prot_alt_name}\n'\
                            f'{self.prot_alt}\n'\
                            f'>{self.prot_main_intersected_name}\n'\
                            f'{self.prot_main_intersected}\n'\
                            f'{self.pdb_main}\n'\
                            f'>{self.prot_alt_intersected_name}\n'\
                            f'{self.prot_alt_intersected}\n'\
                            f'{self.pdb_alt}\n'\
                            f'======================\n'
                            
                            
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        
    parser.add_argument('-predicted_fasta', type=str,
                        help='The path to the input predicted proteins from the EMBOSS getorf find 1.')
    parser.add_argument('-pdb_fasta', type=str,
                        help='The path to the blastp fasta.')
    parser.add_argument('-pdb_tsv', type=str, help='A Path to the input blastp tsv file.')
    parser.add_argument('-output', type=str, default='output_annotated', help='The path to the output dir.')
    
    args = parser.parse_args()
    
    
#     class Args():

#         def __init__(self):

#             self.predicted_fasta = 'predicted_frameshifts_Scler/Scler_600orf_300intrsct.fasta'
#             self.pdb_fasta = 'blastp/pdb_fasta/Scler_600_300_001.fasta'
#             self.pdb_tsv = 'blastp/tsv_data/Scler_600_300_001.tsv'
#             self.output = 'output'
            
#     args = Args()
    
    pdb_fasta, translate_dict, predicted_fasta = parse_data(args)
    
#     try:
#         os.makedir(args.ouput+'_dir')
#     except:
#         pass
            
    with open(args.output+'_both_pdb.fasta', 'w') as w1,\
    open(args.output+'_at_least_1_pdb.fasta', 'w') as w2:
        for batch in predicted_fasta.split('\n=========================================================\n'):
            my_Batch = Batch(batch)
            my_Batch.get_items()
            my_Batch.get_pdb(pdb_fasta, translate_dict)
            
            if ('NOT_INFO' not in my_Batch.pdb_main and 'NOT_INFO' not in my_Batch.pdb_alt):
                my_Batch.combine_batch()
                w1.write(my_Batch.full_batch)
                
            if ('NOT_INFO' not in my_Batch.pdb_main or 'NOT_INFO' not in my_Batch.pdb_alt):
                my_Batch.combine_batch()
                w2.write(my_Batch.full_batch)

