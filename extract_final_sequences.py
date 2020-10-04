#/usr/bin/python
  
from Bio import SeqIO
import sys
 
arg1=sys.argv[1]
arg2=sys.argv[2]
arg3=sys.argv[3]
 
filename = arg1
reads_dict=SeqIO.to_dict(SeqIO.parse(arg2, "fasta"))
with open(arg3, "w") as output_file:
        n=0
        for record in SeqIO.parse(filename, "fasta"):
                if record.id in reads_dict:
                        #n=n+1
                        SeqIO.write(record, output_file, "fasta")
                else:
                        #SeqIO.write(record, output_file, "fasta")
                        n=n+1
 
print(n)
