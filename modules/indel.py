#!/usr/bin/python3
import argparse

parser = argparse.ArgumentParser(description='Remove insertions in multiple alignment')
parser.add_argument('-f', '--fasta', dest='fastaFile', type=str, help='input fasta file of a multiple alignment', required=True)
parser.add_argument('-p','--percent', dest='per', type=int, default=80, help="percentage of gaps to remove a position - default 80")
parser.add_argument('-o','--output', dest='output_file', type=str, default="remove_inserts.fa", help="output file")
args = parser.parse_args()

n=0
dic={}
tail_seq=0

fil=open(args.fastaFile, "r")
for line in fil.readlines() :
    if str(line)[0] == '>':
        fid=str(line.splitlines()[0])
        dic.setdefault(fid, [])
        tail_seq=0
        n=n+1
    else:
        for car in str(line):
            if car!='\n' :
                dic[fid].append(car)
tail_seq=len(dic[fid])
fil.close()

#Adjust percentage for low number of sequences
tail_clust=n
per=args.per
if tail_clust < 5 :
    per=int(int(tail_clust)-1)/int(tail_clust)*100
i=1
#Remove positions with gaps
i=0
while i < tail_seq :
    count=0
    for el in dic:
        if str(dic[el][i]) == "-" :
            count=count+1
    if int(count/tail_clust*100) >= int(per) :
        for el in dic:
            del dic[el][i]
        tail_seq=tail_seq-1
    else :
        i=i+1

f = open(args.output_file, "w")
for ele in dic :
    f.write(ele + "\n")
    seq=''.join(dic[ele])
    f.write(seq + "\n")
f.close()
