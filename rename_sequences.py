#!/usr/bin/python
  
import sys
 
arg1=sys.argv[1]
arg2=sys.argv[2]
arg3=sys.argv[3]
 
fasta=open(arg1)
names=open(arg2)
newfasta=open(arg3,'w')
 
for line in fasta:
    if line.startswith('>'):
        newname=names.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)
 
fasta.close()
names.close()
newfasta.close()
