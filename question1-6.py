def read_fasta_file(fastafile):#receives a path to a file as a string, reads a fasta file with a single DNA sequence and returns a tuple with the id and the sequence
    fasta=open(fastafile, "r")
    fasta.seek(0,0)
    header=""
    seq=""
    for line in fasta:
        if line.startswith(">"):
            header=line.strip()
        else:
            seq=seq+line.strip()
    idandseq=(header,seq)
    return idandseq
    
fasta=raw_input("enter path to fasta file ")
print read_fasta_file(fasta)


def reverse_complement(seq):#receives a DNA sequence as a string and returns its reverse complement.
    revcomp=''
    for nucleotide in seq.upper():
            if nucleotide=='A':
                nuc='t'
            elif nucleotide=='C':
                nuc='g'
            elif nucleotide=='G':
                nuc='c'
            elif nucleotide=='T':
                nuc='a'
            else:
                print nucleotide, "is not a nucleotide base"
            revcomp=nuc+revcomp#reversion
    return revcomp
            
seq=raw_input("enter seq ")
print reverse_complement(seq)


def get_ORF(sequence):#receives a DNA sequence as a string and returns the first open reading frame(0,1 or 2) where both start and stop codons are present; or -1 if it can’t be found, making sure the start codon appears first
    sequence=sequence.upper()
    
    for i in range(0,3):
        rangelist=range(i,len(sequence),3)
        codons=[]
        for index in rangelist:
            codon=sequence[index:index+3]
            ORF=i
            codons=codons+[codon]
            start=""
            for aa in codons:
                if aa=="ATG":
                    start=aa        
                if start!="" and (aa=="TAA" or aa=="TAG" or aa=="TGA"):
                    print "start and terminal codon found in ORF 0"
                    return ORF 
    else:
        print -1
        print "no ORF found"
        
sequence=raw_input("enter DNA sequence ")
print get_ORF(sequence)

def get_ORF(sequence,ORF):#receives a DNA sequence as a string and a number representing an open reading frame and returns the substring between the first start codon in it and the last stop codon.
    sequence=sequence.upper()
    print sequence
    rangelist=range(ORF,len(sequence),3)
    gene=""
    start=""
    for index in rangelist:
        aa=sequence[index:index+3]
        if aa=="ATG":
            gene=gene+aa
        
        if gene!="":
            if aa!="ATG":             
                gene=gene+aa
                if aa=="TAA" or aa=="TGA" or aa=="TAG":
                    return gene
    
sequence=raw_input("enter seq ")
ORF=input("enter ORF ")
print get_ORF(sequence,ORF)



