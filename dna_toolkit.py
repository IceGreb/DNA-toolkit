# %%
NucleoBases=['A','C','G','T']

codons={
'UUU': 'F',      'CUU': 'L',      'AUU': 'I',      'GUU': 'V',
'UUC': 'F',      'CUC': 'L',      'AUC': 'I',      'GUC': 'V',
'UUA': 'L',      'CUA': 'L',      'AUA': 'I',      'GUA': 'V',
'UUG': 'L',      'CUG': 'L',      'AUG': 'M',      'GUG': 'V',
'UCU': 'S',      'CCU': 'P',      'ACU': 'T',      'GCU': 'A',
'UCC': 'S',      'CCC': 'P',      'ACC': 'T',      'GCC': 'A',
'UCA': 'S',      'CCA': 'P',      'ACA': 'T',      'GCA': 'A',
'UCG': 'S',      'CCG': 'P',      'ACG': 'T',      'GCG': 'A',
'UAU': 'Y',      'CAU': 'H',      'AAU': 'N',      'GAU': 'D',
'UAC': 'Y',      'CAC': 'H',      'AAC': 'N',      'GAC': 'D',
'UAA': 'Stop',   'CAA': 'Q',      'AAA': 'K',      'GAA': 'E',
'UAG': 'Stop',   'CAG': 'Q',      'AAG': 'K',      'GAG': 'E',
'UGU': 'C',      'CGU': 'R',      'AGU': 'S',      'GGU': 'G',
'UGC': 'C',      'CGC': 'R',      'AGC': 'S',      'GGC': 'G',
'UGA': 'Stop',   'CGA': 'R',      'AGA': 'R',      'GGA': 'G',
'UGG': 'W',      'CGG': 'R',      'AGG': 'R',      'GGG': 'G'}

# %%
def Validate_Sequence(seq):
    tmp_seq=seq.upper()
    for nuc in tmp_seq:
        if nuc not in NucleoBases:
            raise ValueError("Wrong Sequence Input!")
    return tmp_seq

# %%
from Bio.Seq import Seq
def Base_Count(seq):
    seq=Seq(seq)
    base_count={'A':seq.count('A'), 'C':seq.count("C"), 'G':seq.count("G"), 'T':seq.count("T")}
    return base_count

#%%
def GC_ratio(bases):
    gc_ratio = round((bases['G']+bases['C'])/sum(bases.values())*100)
    return gc_ratio

#%%
def reversed_complement(seq):
    t=seq.replace("A", '%temp').replace("T", "A").replace('%temp', "T").replace("G",'%temp').replace("C","G").replace("%temp", "C")
    complementary_s=t[::-1]
    return complementary_s

#%%
def DNA_transcription(seq):
    rna=seq.replace("T","U")
    return rna

#%%
def DNA_translation(seq):
    for i in range(0,len(seq)+3,3):
        triple = seq[i-3:i]
        if triple in codons:
            if codons[triple]!='Stop':
                print(codons[triple], end='')