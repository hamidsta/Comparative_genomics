import os
import pandas as pd
import time
import argparse
import pathlib
parser = argparse.ArgumentParser(description="calculate the Ks values for all the pairs that belong to the same family")
parser.add_argument("-p","--pairs",type = pathlib.Path,help="specify the path of the file that contains all the pairs",default="paires")
parser.add_argument("-ps","--protseq",type = pathlib.Path,help="specify the path of the file that contains all the protein sequences",default="Prunus_persica.Prunus_persica_NCBIv2.pep.all.fa")
parser.add_argument("-gs","--geneseq",type = pathlib.Path,help="specify the path of the file that contains all the cds sequences",default="Prunus_persica.Prunus_persica_NCBIv2.cds.all.fa")

args = parser.parse_args()
pair_path = args.pairs
protseq_path = args.protseq
geneseq_path = args.geneseq

# create a file: Ks to put the protein id and the Ks values
os.system('rm Pair_Ks')
f_w = open('Pair_Ks','x')

# set the column name
f_w = open('Pair_Ks','a')
f_w.write('Protein1'+'\t'+'Protein2'+'\t'+'Ks'+'\n')
f_w.close()

# read the protein pairs file
df = pd.read_csv(pair_path,sep='\t',header=None)

# a function for select the protein sequences of one pair
def select_protein(p1,p2):
    f_p = open('prot.fst','x')
    with open(protseq_path,'r') as pro_sequence:
        lines = pro_sequence.readlines()
        for i in range(0,len(lines)):
            if lines[i][:1]=='>' and (lines[i][1:9] == p1 or lines[i][1:9] ==p2):
                count = 0
                f_p = open("prot.fst", "a")
                f_p.write(lines[i+count])
                f_p.close()
                count = count+1
                while (i+count)<len(lines) and lines[i+count][:1]!='>':
                    f_p = open("prot.fst", "a")
                    f_p.write(lines[i+count])
                    f_p.close()
                    count = count+1       
    

# # a function for select the gene sequences of one pair
def select_gene(p1,p2):
    f_p = open('cds.fst','x')
    with open(geneseq_path,'r') as pro_sequence:
        lines = pro_sequence.readlines()
        for i in range(0,len(lines)):
            if lines[i][:1]=='>' and (lines[i][1:9] == p1 or lines[i][1:9] ==p2):
                count = 0
                f_p = open("cds.fst", "a")
                f_p.write(lines[i+count])
                f_p.close()
                count = count+1
                while (i+count)<len(lines) and lines[i+count][:1]!='>':
                    f_p = open("cds.fst", "a")
                    f_p.write(lines[i+count])
                    f_p.close()
                    count = count+1

# for each pair, calculate the Ks value
start_time = time.time()
for index,row in df.iterrows():
    p1 = row[0]
    p2 = row[1]
    os.system('rm prot.fst')
    select_protein(p1,p2)
    os.system('rm cds.fst')
    select_gene(p1,p2)
    # step 1: alignment of the proteic sequences
    os.system('clustalw2 -quiet -align -infile=prot.fst -outfile=prot.ali.aln')
    # step 2: cds sequences aligned by corresponding
    os.system('./pal2nal.pl prot.ali.aln cds.fst -output paml > cds.ali.phy')
    # step 3: modification of the control file
    os.system("awk -v file=cds.ali.phy '{gsub(\"XXXXX\",file); print $0}' yn00.ctl_master.txt > yn00.ctl")
    # step 4: running yn00
    path = '/Users/meimei/Programs/paml4.9h/bin/yn00'
    os.system(path)
    # get Ks value
    with open('2YN.dS','r') as f_ds:
        lines = f_ds.readlines()
        Ks = lines[2][23:29]
    f_w = open('Pair_Ks','a')
    f_w.write(p1 + '\t'+ p2 +'\t'+ Ks +'\n')
    f_w.close()
end_time = time.time()
print('total excution time:{}'.format(end_time-start_time))


