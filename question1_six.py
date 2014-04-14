#QUESTION1
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

#QUESTION2
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

#QUESTION3
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
                    return ORF 
    else:
        print -1
        print "no ORF found"
        
sequence=raw_input("enter DNA sequence ")
print get_ORF(sequence)

#QUESTION4
def get_gene_by_ORF(sequence,ORF):#receives a DNA sequence as a string and a number representing an open reading frame and returns the substring between the first start codon in it and the last stop codon.
    sequence=sequence.upper()
    rangelist=range(ORF,len(sequence),3)
    gene=""
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
print get_gene_by_ORF(sequence,ORF)

#QUESTION5
def translate(sequence):#translates a whole DNA sequence(string) to a protein sequence
    sequence=sequence.upper()
    rangelist=range(0,len(sequence),3)
    protseq=""
    for index in rangelist:
            aa=sequence[index:index+3]
            if aa=="ATG":
                codon="M"                           
            if aa=="ATT" or aa=="ATC" or aa=="ATA":
                codon="I"                
            if aa=="CTT" or aa=="CTC" or aa=="CTA" or aa=="CTG" or aa=="TTA" or aa=="TTG":
                codon="L"                
            if aa=="GTT" or aa=="GTC" or aa=="GTA" or aa=="GTG":
                codon="V"                 
            if aa=="TTT" or aa=="TTC":
                codon="F"                
            if aa=="TGT" or aa=="TGC":
                codon="C"                
            if aa=="GCT" or aa=="GCG" or aa=="GCA" or aa=="GCC":
                codon="A"               
            if aa=="GGT" or aa=="GGC" or aa=="GGA" or aa=="GGG":
                codon="G"
                
            if aa=="CCT" or aa=="CCC" or aa=="CCA" or aa=="CCG":
                codon="P"
                
            if aa=="ACT" or aa=="ACC" or aa=="ACA" or aa=="ACG":
                codon="T"
                
            if aa=="TCT" or aa=="TCC" or aa=="TCA" or aa=="TCG" or aa=="AGT" or aa=="AGC":
                codon="S"
                
            if aa=="TAT" or aa=="TAC":
                codon="Y"
                
            if aa=="TGG":
                codon="W"
                
            if aa=="CAG" or aa=="CAA":
                codon="Q"
                
            if aa=="AAT" or aa=="AAC":
                codon="N"
                
            if aa=="CAT" or aa=="CAC":
                codon="H"
                
            if aa=="GAA" or aa=="GAG":
                codon="E"
                
            if aa=="GAT" or aa=="GAG":
                codon="D"
                
            if aa=="AAA" or aa=="AAG":
                codon="K"
                
            if aa=="CGT" or aa=="CGC" or aa=="CGA" or aa=="CGG" or aa=="AGA" or aa=="AGG":
                codon="R"
                  
            if aa=="TAG" or aa=="TAA" or aa=="TGA":
                codon="*stop*"                             
            protseq=protseq+codon
    return protseq           
    
sequence=raw_input("enter sequence ")    
print translate(sequence)

#QUESTION6
def createfasta(filename):#receives a sequence(string) and an ID(string) and creates a string in fasta format where the ID is the content of the first line and the sequence is split into lines of maximum 60 characters
    filecreate=open("newfile.txt","w+")
    filecreate.write(">")
    filecreate.write(seqid,"\n")
    for index in range(0,len(sequence),60):
        filecreate.writelines(sequence[index],"\n")
    filecreate.close()
    
    
    
sequence=raw_input("enter sequence ")
seqid=raw_input("enter sequence identity ")
print createfasta(sequence,seqid) 