import os,sys
from statistics import mean 
#This is to get aggragate the rnaseq data...
#by taking the mean of the gene expression data between the three replicates..
#since the data was already normalized...
try:
    if len(sys.argv) < 2:
        print('USAGE: python path_to_rnaseq')
        exit(0)
except:
    pass

def calcluate_mean_per_replicate(path,files):
    gex = {}
    for f in files:
        fH = open (path+f,'r')
        lines = fH.readlines()
        for line in lines:
            tokens = line.replace('\n','').split('\t')
            if tokens[0] == 'gene_id':
                continue
            elif tokens[0] in gex:
                    gex[tokens[0]].append(float(tokens[1]))
            else:
                gex[tokens[0]] = []
                gex[tokens[0]].append(float(tokens[1]))
    names =  files[0].split('-')
    ofile = names[0]+'-'+names[1]+'.aggregated.txt'
    fOut  = open (ofile,'w')
    fOut.write('gene_id\t{}\n'.format(names[1]))
    for gene in gex:
        M = mean(gex[gene])
        fOut.write('{}\t{}\n'.format(gene,M))
    print('NOTICE: DONE Processing {} files'.format(names[1]))
    return
if __name__ == "__main__":
    files = os.listdir(sys.argv[1])
    files.sort()
    calcluate_mean_per_replicate(sys.argv[1],files[:3])
    calcluate_mean_per_replicate(sys.argv[1],files[3:6])
    calcluate_mean_per_replicate(sys.argv[1],files[6:9])
