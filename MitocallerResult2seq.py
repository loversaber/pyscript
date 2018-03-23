#time python3 MitocallerResult2seq_2.py refDG.fasta fils.txt
from sys import argv
from itertools import islice
import re
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord


seq_record =SeqIO.read(argv[1],"fasta")

def cpr(genetype,frq):
 genetypes,frqs=re.split(':|/',genetype)[1:],re.split(':|/',frq)[1:]
 d=dict(zip(genetypes,frqs))
 d_sort=OrderedDict(sorted(d.items(),key = lambda t:t[1],reverse=True))
 genefrq=list(d_sort.items())[0]
 g,f=genefrq[0],genefrq[1]
 return g,f

l0=list(range(1,19525))
def op(f):
 dpos,lhave={},[]
 with open(f,"r") as mitocall:
  for line in islice(mitocall,1,None):
   #strline=line.strip()#get substitution line
   line=line.strip().split("\t")
   pos,ref,genetype,frq=int(line[1]),line[2],line[30],line[31]
   lhave.append(pos)
   depth0,depth1=line[3],line[11]
   depth=re.split(':',depth1)[1]

   alt,altfrq=cpr(genetype,frq)
   refbase=re.split(':|/',ref)[1]
   if alt!="N" and refbase != alt and int(depth)>=10 and float(altfrq)==1:
    dpos[pos]=[refbase,alt,altfrq,line]
  dpos_sort=OrderedDict(sorted(dpos.items(),key=lambda t:t[0],reverse=False))
  lnot=list(set(l0)-set(lhave))
  lnot.sort()
 return dpos_sort,lnot

def ref2seq(d,l):
 mutseq1=MutableSeq(str(seq_record.seq))
 for k in d.keys():
  if mutseq1[k-1]==d[k][0]:
   #print("###REFSEQ GENE MATCH XLS REF###")
   mutseq1[k-1]=d[k][1]
 #for i in l:#change not get pos to ?
   #mutseq1[i-1]="?"
 return mutseq1

with open(argv[2],"r") as allfils:
 for i in allfils.readlines():
  i=i.strip()
  name=i.split(".")[0]
  otfil=open(name+".ref2seqfromMitocall_1","w")
  subfil=open(name+".sub_1","w")
  daltpos,lnot=op(i)
  for k in daltpos.keys():
   subfil.write("\t".join(daltpos[k][3])+"\n")
  #print(len(daltpos.keys()),daltpos)
  nalt,nnot=len(daltpos.keys()),len(lnot)
  changeseq=ref2seq(daltpos,lnot)
  rec=SeqRecord(changeseq,id=name,description=str(nalt)+" subs "+str(nnot)+" posnot")
  SeqIO.write(rec,otfil,"fasta")
  otfil.close()
  subfil.close()
