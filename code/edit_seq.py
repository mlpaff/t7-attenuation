from Bio import Entrez, SeqIO

Entrez.email = 'Matthew.paff@utexas.edu'

# Download wt T7 genebank record
handle = Entrez.efetch(db='nucleotide', id='NC_001604', rettype='gb', retmode='text')

record = SeqIO.read(handle, 'genbank')
handle.close()

newSeq=record.seq.tomutable()

#point mutations found in the T7Hi background
newSeq[15093] = 'T'
newSeq[24087] = 'A'
newSeq[30860] = 'C'
newSeq[30944] = 'G'
newSeq[32859] = 'A'
newSeq[34974] = 'G'
newSeq[36491] = 'G'

# Replace phi 9 with nonsense promoter along with the EcoRI 
#newSeq[21850:21870] = 'GAATTCagagattacaataa'

# Replace phi9 with nonsense promoter -- include the deletion of gene 8 stop codon -- name file '8st_*'
newSeq[21847:21870] = 'GAATTCcgaagagattacaataa'

# Replace phi 10
newSeq[22886:22909] = 'AAGCTTcgaagagattacaataa'

record.seq=newSeq

SeqIO.write(record, '8st_phi910.gb', 'genbank')