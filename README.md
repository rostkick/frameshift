# frameshift

## Usage
1.
```
python3 RUN.py -pipeline Snakefile.template -genome data/Sclerotinia_sclerotiorum.ASM14694v1.dna.toplevel.fa -scripts_path scripts/ -prefix Scler_ -min_orf 300 -min_int 100 -seq_I 0.75 -seq_L 0.75 -evalue 0.01
```
2.
```
snakemake efetch_sequences
```
