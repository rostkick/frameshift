import itertools
import os
import re
import subprocess
import tempfile
from pybedtools import bedtool
import argparse

from collections import defaultdict
from Bio import SeqIO, Seq


def get_data(args):
    
    fasta_dna_dict = dict()
    orfs_dict = defaultdict(list)

    for record in SeqIO.parse(args.input_dna, 'fasta'):

        spl = re.split(' |_|\[|\]', record.description)

        if int(spl[3]) < int(spl[5]):
            fasta_dna_dict[f'{spl[0]}:ORF{spl[1]}:{spl[3]}:{spl[5]}:+'] = record.seq
            orfs_dict[spl[0]].append(f"ORF{spl[1]}:{spl[3]}:{spl[5]}:+")
        else:
            fasta_dna_dict[f'{spl[0]}:ORF{spl[1]}:{spl[3]}:{spl[5]}:-'] = record.seq
            orfs_dict[spl[0]].append(f"ORF{spl[1]}:{spl[3]}:{spl[5]}:-")
    
    return fasta_dna_dict, orfs_dict

class Intersection:
    
    def __init__(self, args, dna_1_name, dna_1, dna_2_name, dna_2):
        
        self.genome = args.input_genome
        self.min_orf, self.min_dna_inter = args.min_orf, args.interval_length
        
        self.chrom = dna_1_name.split(':')[0]
        
        self.global_orf_name, self.global_orf = '', ''
        
        self.dna_1_name, self.dna_1 = dna_1_name, dna_1
        self.dna_2_name, self.dna_2 = dna_2_name, dna_2  
        
        self.prot_1_name, self.prot_1 = dna_1_name, Seq.translate(dna_1)
        self.prot_2_name, self.prot_2 = dna_2_name, Seq.translate(dna_2)
        
        
        self.inter_dna_1_name, self.inter_dna_1 = '', ''
        self.inter_dna_2_name, self.inter_dna_2 = '', ''
        
        self.inter_prot_1_name, self.inter_prot_1 = '', ''
        self.inter_prot_2_name, self.inter_prot_2 = '', ''
        
        self.status = False
        
    def _get_interval(self, x, y, strand):
            
            if strand == "+":
                a = bedtool.BedTool(f'{self.chrom}\t{x-1}\t{y}\t.\t0\t{strand}',
                                    from_string=True)
                a = a.sequence(fi=args.input_genome, name=True, s=True)
                a = open(a.seqfn).read().lstrip('>').split('\n')
                
            elif strand == '-':
                a = bedtool.BedTool(f'{self.chrom}\t{y-1}\t{x}\t.\t0\t{strand}',
                                    from_string=True)
                a = a.sequence(fi=args.input_genome, s=True)
                a = open(a.seqfn).read().lstrip('>').split('\n')
            
            return a[1]
        
    def process_intersection(self):
        
        d1, d2 = self.dna_1_name.split(':'), self.dna_2_name.split(':')
        orf1, orf2 = d1[1], d2[1]
        x1, y1, x2, y2 = int(d1[2]), int(d1[3]), int(d2[2]), int(d2[3])
        strand1, strand2 = d1[4], d2[4]
        
        if all([strand1=='+', strand2=='+', x1<y2, y1>x2,\
                all([abs(x1-y2)>=self.min_dna_inter, abs(x2-y1)>=self.min_dna_inter])]):
            
            if all([x1>x2, y1>y2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':LEFT_INTERSECTION'
                self.dna_2_name = self.dna_2_name + ':RIGHT_INTERSECTION'
                
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y2}:+'
                self.inter_dna_1 = self._get_interval(x1, y2, '+')
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y2}:+'
                self.inter_dna_2 = self.inter_dna_1
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y2}:+'
                self.inter_prot_1 = self.prot_1[0:abs(y2-x1)//3]
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y2}:+'
                self.inter_prot_2 = self.prot_2[abs(x1-x2)//3:(abs(x1-x2)+abs(y2-x1))//3]
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y1}:+'
                self.global_orf = self._get_interval(x2, y1, '+')
                
                self.status = True
                    
            elif all([x1<x2, y1<y2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':RIGHT_INTERSECTION'
                self.dna_2_name = self.dna_2_name + ':LEFT_INTERSECTION'
                
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y1}:+'
                self.inter_dna_1 = self._get_interval(x2, y1, '+')
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y1}:+'
                self.inter_dna_2 = self.inter_dna_1
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y1}:+'
                self.inter_prot_1 = self.prot_1[abs(x2-x1)//3:(abs(x2-x1)+abs(y1-x2))//3]
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y1}:+'
                self.inter_prot_2 = self.prot_2[0:abs(y1-x2)//3]
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y2}:+'
                self.global_orf = self._get_interval(x1, y2, '+')
                
                self.status = True
                    
            elif all([x1<x2, y1>y2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':OUTSIDE'
                self.dna_2_name = self.dna_2_name + ':INSIDE'
                
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:+'
                self.inter_dna_1 = self.dna_2
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:+'
                self.inter_dna_2 = self.dna_2
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:+'
                self.inter_prot_1 = self.prot_1[abs(x2-x1)//3:(abs(x2-x1)+abs(y2-x2))//3]
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:+'
                self.inter_prot_2 = self.prot_2
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:+'
                self.global_orf = self.dna_1
                
                self.status = True
                
            elif all([x1<x2, y1<y2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':INSIDE'
                self.dna_2_name = self.dna_2_name + ':OUTSIDE'

                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:+'
                self.inter_dna_1 = self.dna_1
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:+'
                self.inter_dna_2 = self.dna_1
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{x1}:+'
                self.inter_prot_1 = self.prot_1
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:+'
                self.inter_prot_2 = self.prot_2[abs(x1-x2)//3:(abs(x1-x2)+abs(y1-x1))//3]
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:+'
                self.global_orf = self.dna_2
                
                self.status = True

        elif all([strand1=='-', strand2=='-', x1>y2, y1<x2,\
                all([abs(x1-y2)>=self.min_dna_inter, abs(x2-y1)>=self.min_dna_inter])]):
            
            if all([x1>x2, y1>y2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':LEFT_INTERSECTION'
                self.dna_2_name = self.dna_2_name + ':RIGHT_INTERSECTION'

                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y1}:-'
                self.inter_dna_1 = self._get_interval(x2, y1, '-')
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y1}:-'
                self.inter_dna_2 = self.inter_dna_1
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y1}:-'
                self.inter_prot_1 = self.prot_1[abs(x2-x1)//3:(abs(x2-x1)+abs(y1-x2))//3]
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y2}:{y1}:-'
                self.inter_prot_2 = self.prot_2[0:abs(y1-x2)//3]
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y2}:{x1}:+'
                self.global_orf = self._get_interval(y2, x1, '+')
                
                self.status = True
                    
            elif all([x1<x2, y1<y2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')

                self.dna_1_name = self.dna_1_name + ':RIGHT_INTERSECTION'
                self.dna_2_name = self.dna_2_name + ':LEFT_INTERSECTION'
                
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y2}:-'
                self.inter_dna_1 = self._get_interval(x1, y2, '-')
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y2}:-'
                self.inter_dna_2 = self.inter_dna_1
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y2}:-'
                self.inter_prot_1 = self.prot_1[0:abs(y2-x1)//3]
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y2}:-'
                self.inter_prot_2 = self.prot_2[abs(x1-x2)//3:(abs(x1-x2)+abs(y2-x1))//3]
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y1}:{x2}:+'
                self.global_orf = self._get_interval(y1, x2, '+')
                
                self.status = True   
                
            elif all([x1<x2, y1>y2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':INSIDE'
                self.dna_2_name = self.dna_2_name + ':OUTSIDE'
                
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:-'
                self.inter_dna_1 = self.dna_1
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:-'
                self.inter_dna_2 = self.dna_1
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:-'
                self.inter_prot_1 = self.prot_1
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:-'
                self.inter_prot_2 = self.prot_2[abs(x1-x2)//3:(abs(x1-x2)+abs(y1-x1))//3]
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:-'
                self.global_orf = self._get_interval(x2, y2, '-')
                
                self.status = True    
                
            elif all([x1>x2, y1<y2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':OUTSIDE'
                self.dna_2_name = self.dna_2_name + ':INSIDE'
                    
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:-'
                self.inter_dna_1 = self.dna_2
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:-'
                self.inter_dna_2 = self.dna_2
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:-'
                self.inter_prot_1 = self.prot_1[abs(x2-x1)//3:(abs(x2-x1)+abs(y2-x2))//3]
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:-'
                self.inter_prot_2 = self.prot_2 
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:-'
                self.global_orf = self.dna_1
                
                self.status = True
                
        elif all([strand1=='+', strand2=='-', x1<x2, y1>y2,\
                all([abs(x1-x2)>=self.min_dna_inter, abs(y1-y2)>=self.min_dna_inter])]):
            
            if all([x1>y2, y1>x2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':LEFT_INTERSECTION'
                self.dna_2_name = self.dna_2_name + ':RIGHT_INTERSECTION'
                
                
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{x2}:+'
                self.inter_dna_1 = self._get_interval(x1, x2, '+')
                
                self.inter_dna_2 = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{x1}:-'
                self.inter_dna_2 = self._get_interval(x2, x1, '-')
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{x2}:+'
                self.inter_prot_1 = self.prot_1[0:abs(x2-x1)//3]
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{x1}:-'
                self.inter_prot_2 = self.prot_2[0:abs(x1-x2)//3]
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y2}:{y1}:+'
                self.global_orf = self._get_interval(y2, y1, "+")
                
                self.status = True  
                
            elif all([x1<y2, y1<x2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')

                self.dna_1_name = self.dna_1_name + ':RIGHT_INTERSECTION'
                self.dna_2_name = self.dna_2_name + ':LEFT_INTERSECTION'
                
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y2}:{y1}:+'
                self.inter_dna_1 = self._get_interval(y2, y1, '+')
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y1}:{y2}:-'
                self.inter_dna_2 = self._get_interval(y1, y2, '-')
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y2}:{y1}:+'
                self.inter_prot_1 = self.prot_1[abs(y2-x1)//3:(abs(y2-x1)+abs(y1-y2))//3]
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y1}:{y2}:-'
                self.inter_prot_2 = self.prot_2[abs(y1-x2)//3:(abs(y1-x2)+abs(y2-y1))//3]
                             
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{x2}:+'
                self.global_orf = self._get_interval(x1, x2, "+")
                
                self.status = True
                
            elif all([x1<y2, y1>x2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':OUTSIDE'
                self.dna_2_name = self.dna_2_name + ':INSIDE'
                
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y2}:{x2}:+'
                self.inter_dna_1 = self._get_interval(y2, x2, '+')
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:-'
                self.inter_dna_2 = self.dna_2
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y2}:{x2}:+'
                self.inter_prot_1 = self.prot_1[abs(y2-x1)//3:(abs(y2-x1)+abs(x2-y2))//3]
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:-'
                self.inter_prot_2 = self.prot_2
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:+'
                self.global_orf = self._get_interval(x1, y1, "+")
                
                self.status = True    
                
            elif all([x1>y2, y1<x2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':INSIDE'
                self.dna_2_name = self.dna_2_name + ':OUTSIDE'

                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:+'
                self.inter_dna_1 = self.dna_1
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y1}:{x1}:-'
                self.inter_dna_2 = self._get_interval(y1, x1, '-')
                
                self.inter_prot_1_name =f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:+'
                self.inter_prot_1 = self.prot_1
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y1}:{x1}:-'
                self.inter_prot_2 = self.prot_2[abs(y1-x2)//3:(abs(y1-x2)+abs(x1-y1))//3]
                    
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:-'
                self.global_orf = self._get_interval(x2 ,y2, "-")
                
                self.status = True
                
        elif all([strand1=='-', strand2=='+', x1>x2, y1<y2,\
                all([abs(x1-x2)>=self.min_dna_inter, abs(y1-y2)>=self.min_dna_inter])]):
            
            if all([x1<y2, y1<x2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':RIGHT_INTERSECTION'
                self.dna_2_name = self.dna_2_name + ':LEFT_INTERSECTION'
                
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{x2}:-'
                self.inter_dna_1 = self._get_interval(x1, x2, "-")
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{x1}:+'
                self.inter_dna_2 = self._get_interval(x2, x1, '+')
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{x2}:-'
                self.inter_prot_1 = self.prot_1[0:abs(x2-x1)//3]
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{x1}:+'
                self.inter_prot_2 = self.prot_2[0:abs(x1-x2)//3]
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y1}:{y2}:+'
                self.global_orf = self._get_interval(y1 ,y2, "+") 
                
                self.status = True
                
            elif all([x1>y2, y1>x2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':LEFT_INTERSECTION'
                self.dna_2_name = self.dna_2_name + ':RIGHT_INTERSECTION'
                
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y2}:{y1}:-'
                self.inter_dna_1 = self._get_interval(y2, y1, "-")
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y1}:{y2}:+'
                self.inter_dna_2 = self._get_interval(y1, y2, '+')
                
                self.inter_prot_1 = f'{self.chrom}:{orf1}_vs_{orf2}:{y2}:{y1}:-'
                self.inter_prot_1 = self.prot_1[abs(y2-x1)//3:(abs(y2-x1)+abs(x2-y2))//3]
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y1}:{y2}:+'
                self.inter_prot_2 = self.prot_2[abs(y1-x2)//3:(abs(y1-x2)+abs(y2-y1))//3]
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{x1}:+'
                self.global_orf = self._get_interval(x2 ,x1, "+")
                
                self.status = True
                
            elif all([x1>y2, y1<x2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':OUTSIDE'
                self.dna_2_name = self.dna_2_name + ':INSIDE'
                
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y2}:{x2}:-'
                self.inter_dna_1 = self._get_interval(y2, x2, "-")
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:+'
                self.inter_dna_2 = self.dna_2
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y2}:{x2}:-'
                self.inter_prot_1 = self.prot_1[(y2-x1)//3:(y2-x1+x2-y2)//3]
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:+'
                self.inter_prot_2 = self.prot_2
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:-'
                self.global_orf = self._get_interval(x1 ,y1, "-")
                
                self.status = True
                
            elif all([x1<y2, y1>x2]):
                
                print(f'process: {self.chrom} {orf1} vs {orf2}')
                
                self.dna_1_name = self.dna_1_name + ':INSIDE'
                self.dna_2_name = self.dna_2_name + ':OUTSIDE'
                
                self.inter_dna_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:-'
                self.inter_dna_1 = self.dna_1
                
                self.inter_dna_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y1}:{x1}:+'
                self.inter_dna_2 = self._get_interval(y1, x1, "+")
                
                self.inter_prot_1_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x1}:{y1}:-'
                self.inter_prot_1 = self.prot_1
                
                self.inter_prot_2_name = f'{self.chrom}:{orf1}_vs_{orf2}:{y1-x2}:{y1-x2+x1-y1}:+'
                self.inter_prot_2 = self.prot_2[abs(y1-x2)//3:(abs(y1-x2)+abs(x1-y1))//3]
                
                self.global_orf_name = f'{self.chrom}:{orf1}_vs_{orf2}:{x2}:{y2}:+'
                self.global_orf = self.dna_2
                
                self.status = True
        
    def get_full_batch(self):
        
        return f">GLOBAL_ORF:{self.global_orf_name}\n{self.global_orf}\n"\
                f">DNA_MAIN:{self.dna_1_name}\n{self.dna_1}\n"\
                f">DNA_ALT:{self.dna_2_name}\n{self.dna_2}\n"\
                f">INTERSECTION_DNA_MAIN:{self.inter_dna_1_name}\n{self.inter_dna_1}\n"\
                f">INTERSECTION_DNA_ALT:{self.inter_dna_2_name}\n{self.inter_dna_2}\n"\
                f">PROT_MAIN:{self.prot_1_name}\n{self.prot_1}\n"\
                f">PROT_ALT:{self.prot_2_name}\n{self.prot_2}\n"\
                f">PROT_INTERSECTED_MAIN:{self.inter_prot_1_name}\n{self.inter_prot_1}\n"\
                f">PROT_INTERSECTED_ALT:{self.inter_prot_2_name}\n{self.inter_prot_2}\n"\
                f"=========================================================\n"
     
        
def main(args):
    
    fasta_dna_dict, orfs_dict = get_data(args)
    
    with open(args.output, 'w') as w:
        for chrom, orfs in orfs_dict.items():
            for dna1_name, dna2_name in itertools.combinations(orfs, 2):
                if all([len(fasta_dna_dict[f"{chrom}:{dna1_name}"])>=args.min_orf,
                        len(fasta_dna_dict[f"{chrom}:{dna2_name}"])>=args.min_orf]):

                    Inter = Intersection(args,
                                         f"{chrom}:{dna1_name}", fasta_dna_dict[f"{chrom}:{dna1_name}"],\
                                         f"{chrom}:{dna2_name}", fasta_dna_dict[f"{chrom}:{dna2_name}"])
                    Inter.process_intersection()
                    if Inter.status:
                        w.write(Inter.get_full_batch())

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-input_dna', type=str,
                        help='The address to the input predicted dna from the EMBOSS getorf find 3.')
    parser.add_argument('-input_genome', type=str, help='A Path to the input genome.')
    parser.add_argument('-min_orf', type=int, default=300, help='Minimum ORF length.')
    parser.add_argument('-interval_length', type=int, default=100, help='Minimum length of interval of intersection between 2 ORF.')
    parser.add_argument('-output', type=str, default='output.fasta', help='The address to the output.')

    args = parser.parse_args()

#     class Args():
#         def __init__(self):
#             self.input_dna = './input_predicted_prots_and_dna/Scler_dna_pred'
#             self.input_genome = './Sclerotinia_sclerotiorum_orgin_data/Sclerotinia_sclerotiorum.ASM14694v1.dna.toplevel.fa'
#             self.min_orf = 300
#             self.output = 'output'
#             self.interval_length = 300

#    args = Args()

    main(args)
