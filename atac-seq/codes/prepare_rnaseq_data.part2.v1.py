import os, sys
#This is to create bed file ...
try:
    if len(sys.argv) < 3:
        print('USAGE:python mart_export.txt cell_aggregated.txt')
        exit(0)
except:
    pass

def process_mart_file(mart_file):
    fH = open (mart_file,'r')
    lines = fH.readlines()
    mart_data = {}
    del lines[0]
    for line in lines:
        tokens = line.replace('\n','').split('\t')
        mart_data[tokens[0]] = 'chr'+tokens[-1] +'\t'+tokens[1]+'\t'+tokens[2]
    return mart_data

def process_aggregated_file(aggregated_file):
    fH = open (aggregated_file,'r')
    lines = fH.readlines()
    del lines[0]
    genes  ={}
    for line in lines:
        tokens = line.split('\t')
        parts = tokens[0].split('.')
        genes[parts[0]] = tokens[1]

    return genes

def intersect(mart_data,genes):
    common_genes = list(mart_data.keys() & genes.keys())
    print('NOTICE: Found {} annotated genes out-of {}...'.format(len(common_genes),len(genes)))
    ofile1 = 'all_genes.bed'
    ofile2 = 'all_genes.txt'
    fOut1 = open (ofile1,'w')
    fOut2 = open (ofile2,'w')
    for gene in common_genes:
        fOut1.write('{}\n'.format(mart_data[gene]))
        fOut2.write('{}\t{}\n'.format(gene,mart_data[gene]))
    return
#
if __name__ == "__main__":
    mart_data = process_mart_file(sys.argv[1])
    genes = process_aggregated_file(sys.argv[2])
    intersect(mart_data,genes)
