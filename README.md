
from Bio.Seq import Seq
my_seq = Seq("GATCG")
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
0 G
1 A
2 T
3 C
4 G
# We can also print the length of each sequence
print(len(my_seq))
5
# We can print specific positions in the sequence
print(my_seq[0])
G
print(my_seq[4])
G
print(my_seq[2])
T
#We can use a count function to find number of patterns in sequence
Seq("AAAA").count("AA")
2
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
len(my_seq)
32
my_seq.count("G")
9
#We can use count function for mathematical calculations, such as finding GC content
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
46.875
# Biopython has common mathematical calculations built in as functions, for example, GC content
from Bio.SeqUtils import gc_fraction
gc_fraction(my_seq)
0.46875
# We can also use biopython to slice sequences into parts
my_seq[4:12]
Seq('GATGGGCC')
my_seq[0::3]
Seq('GCTGTAGTAAG')
my_seq[1::3]
Seq('AGGCATGCATC')
my_seq[2:3]
Seq('T')
# Like strings, we can also use negatives to start from the end
my_seq[::-1]
Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG')
# We can turn sequences back into strings
str(my_seq)
'GATCGATGGGCCTATATAGGATCGAAAATCGC'
# You can have a verb placeholder string for the sequence
fasta_format_string = ">Name\n%s\n" % my_seq
print(fasta_format_string)
>Name
GATCGATGGGCCTATATAGGATCGAAAATCGC

# We can also add sequences together
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
seq1 + seq2
Seq('ACGTAACCGG')
seq2 + seq1
Seq('AACCGGACGT')
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
# N = unable to determine the exact nucleotide
spacer = Seq("N" *10)
# We can join the contigs together with a spacer inbetween
spacer.join(contigs)
Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')
# Some sequences may appear case sensitive
dna_seq = Seq("acgtACGT")
# Using commands we can make entire sequences upper or lower case
dna_seq.upper()
Seq('ACGTACGT')
dna_seq.lower()
Seq('acgtacgt')
dna_seq = dna_seq.upper()
dna_seq
Seq('ACGTACGT')
"gtac" in dna_seq
False
"GTAC" in dna_seq
True
# We can the complement and reverse complement sequences for DNA sequences
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
my_seq.complement()
Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG')
my_seq.reverse_complement()
Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC')
# We can also use biopython for protein sequences and their complements (complements are determined using ambiguity codes)
protein_seq = Seq("EVRNAK")
protein_seq.complement()
Seq('EBYNTM')
# We can also use biopython for DNA/RNA sequences
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
# We can get a template strand by finding the reverse complement
template_dna = coding_dna.reverse_complement()
print(template_dna)
CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
# Biopython also has a function for DNA to RNA transcription
messenger_rna = coding_dna.transcribe()
messenger_rna
Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
# Biological transcription with template strand becomes a 2-step process
template_dna.reverse_complement().transcribe()
Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')
# We can also do reverse transcription
messenger_rna.back_transcribe()
Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')
# Biopython can also translate RNA sequences into proteins
# * = stop codons
messenger_rna.translate()
Seq('MAIVMGR*KGAR*')
# We can also specify codon tables with biopython (nuclear/mitochondrial)
coding_dna.translate("Vertebrate Mitochondrial")
​
Seq('MAIVMGRWKGAR*')
# Also works with NCBI table numbers
coding_dna.translate(table = 2)
Seq('MAIVMGRWKGAR*')
# We can choose to translate just to the first stop codon
coding_dna.translate(to_stop = True)
Seq('MAIVMGR')
coding_dna.translate(table=2, to_stop=True)
Seq('MAIVMGRWKGAR')
# We can specify our stop codons with different symbols
coding_dna.translate(table=2, stop_symbol="!")
Seq('MAIVMGRWKGAR!')
# What if sequence uses a nonstandard codon?
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
gene.translate(table = "Bacterial")
Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')
gene.translate(table = "Bacterial", to_stop = True)
Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')
# We need to tell biopython that it is a complete CDS to start with methionine
gene.translate(table = "Bacterial", cds = True)
Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')
# We can look at codon tables through biopython
from Bio.Data import CodonTable
# Standard = Table 1
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
# Vertebrate Mitochondrial = Table 2
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
print(standard_table)
Table 1 Standard, SGC0

  |  T      |  C      |  A      |  G      |
--+---------+---------+---------+---------+--
T | TTT F   | TCT S   | TAT Y   | TGT C   | T
T | TTC F   | TCC S   | TAC Y   | TGC C   | C
T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
--+---------+---------+---------+---------+--
C | CTT L   | CCT P   | CAT H   | CGT R   | T
C | CTC L   | CCC P   | CAC H   | CGC R   | C
C | CTA L   | CCA P   | CAA Q   | CGA R   | A
C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
--+---------+---------+---------+---------+--
A | ATT I   | ACT T   | AAT N   | AGT S   | T
A | ATC I   | ACC T   | AAC N   | AGC S   | C
A | ATA I   | ACA T   | AAA K   | AGA R   | A
A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
--+---------+---------+---------+---------+--
G | GTT V   | GCT A   | GAT D   | GGT G   | T
G | GTC V   | GCC A   | GAC D   | GGC G   | C
G | GTA V   | GCA A   | GAA E   | GGA G   | A
G | GTG V   | GCG A   | GAG E   | GGG G   | G
--+---------+---------+---------+---------+--
print(mito_table)
Table 2 Vertebrate Mitochondrial, SGC1

  |  T      |  C      |  A      |  G      |
--+---------+---------+---------+---------+--
T | TTT F   | TCT S   | TAT Y   | TGT C   | T
T | TTC F   | TCC S   | TAC Y   | TGC C   | C
T | TTA L   | TCA S   | TAA Stop| TGA W   | A
T | TTG L   | TCG S   | TAG Stop| TGG W   | G
--+---------+---------+---------+---------+--
C | CTT L   | CCT P   | CAT H   | CGT R   | T
C | CTC L   | CCC P   | CAC H   | CGC R   | C
C | CTA L   | CCA P   | CAA Q   | CGA R   | A
C | CTG L   | CCG P   | CAG Q   | CGG R   | G
--+---------+---------+---------+---------+--
A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
--+---------+---------+---------+---------+--
G | GTT V   | GCT A   | GAT D   | GGT G   | T
G | GTC V   | GCC A   | GAC D   | GGC G   | C
G | GTA V   | GCA A   | GAA E   | GGA G   | A
G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
--+---------+---------+---------+---------+--
mito_table.stop_codons
['TAA', 'TAG', 'AGA', 'AGG']
mito_table.start_codons
['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
# We can compare sequences
seq = Seq("ACGT")
seq1 == "ACGT"
True
"ACGT" == seq1
True
# Sometimes the length of a sequence is known but not the contents
unknown_seq = Seq(None, 10)
#Python tells you "None" as in the contents are not known but the length is
unknown_seq
Seq(None, length=10)
len(unknown_seq)
10
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
seq[1000:1020]
Seq(None, length=20)
#Python can show you a portion of the sequence that is defined
seq[117512690:117512700]
Seq('CCTGAATGTG')
#You can have partial chromosomal information combined with length
seq[117512670:]
Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)
# We can also create these partially unknown sequences by uphending
seq = Seq("ACGT")
undefined_seq = Seq(None, length = 10)
seq + undefined_seq + seq
Seq({0: 'ACGT', 14: 'ACGT'}, length=18)
# This is crucial to ensure we are not changing our sequence data
my_seq = Seq("GCCATTGTAATGTGGCCGCTGAAAGGGTGCCCGA")
# We can create an equal sequence that allows to change data without affecting original sequence (mutable)
from Bio.Seq import MutableSeq
mutable_seq = MutableSeq(my_seq)
mutable_seq
MutableSeq('GCCATTGTAATGTGGCCGCTGAAAGGGTGCCCGA')
mutable_seq[5] = "C"
# 5th position will now be changed from T to C
mutable_seq
MutableSeq('GCCATCGTAATGTGGCCGCTGAAAGGGTGCCCGA')
# We can also remove data 
mutable_seq.remove("T")
# Removes the first T
mutable_seq
MutableSeq('GCCACGAATGTGGCCGCTGAAAGGGTGCCCGA')
mutable_seq.reverse()
mutable_seq
MutableSeq('AGCCCGTGGGAAAGTCGCCGGTGTAAGCACCG')
new_seq = Seq(mutable_seq)
# It is now back to being protected from editing
mutable_seq
MutableSeq('AGCCCGTGGGAAAGTCGCCGGTGTAAGCACCG')
# We can import packages that work with strings if we needed to use strings 
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
my_string = "CGTTTAGGCGCAAAGTCGAGATCGGATAGGCTAATTCTCGGGATA"
reverse_complement(my_string)
'ATATCCCGAGAATTAGCCTATCCGATCTCGACTTTGCGCCTAAACG'
transcribe(my_string)
'CGUUUAGGCGCAAAGUCGAGAUCGGAUAGGCUAAUUCUCGGGAUAU'
back_transcribe(my_string)
'CGTTTAGGCGCAAAGTCGAGATCGGATAGGCTAATTCTCGGGATAT'
translate(my_string)
'RLGAKSRSDRLILGI'
