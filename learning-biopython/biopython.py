from Bio import SeqIO
from Bio.Seq import Seq

some_seq = Seq('AAATTCGACTG')


for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

print(some_seq.reverse_complement())