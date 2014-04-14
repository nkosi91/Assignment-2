#QUESTION7
#receives 2 arguments from the command line: the path of fasta file that contains a single DNA sequence, and the path of a new file where the program will save the translate protein. Make sure to verify the correct ORF and start & stop of the sequence
import question1_six
def multifun(source, save):
    question1_six.read_fasta_file(source)
    question1_six.get_ORF(idandseq[2])
    question1_six.get_gene_by_ORF(ORF)
    question1_six.translate(gene)
    question1_six.createfasta(save)

print multifun()