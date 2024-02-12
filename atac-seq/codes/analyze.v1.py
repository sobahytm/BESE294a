import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#This is to analyze genes expressions that share genimic intervals with accessibilty data...
try:
    if len(sys.argv) <3:
        print('USAGE: python path_to_analysis path_rnaseq_aggregate')
        exit(0)
except:
    pass

def process_analysis_folder(path):
    files = os.listdir(path)
    genes_intervals = {}
    sig_peaks_intervals = {}
    for file in files:
        if file == 'all_genes.txt':
            fH = open (path+file,'r')
            lines = fH.readlines()
            for line in lines:
                tokens = line.replace('\n','').split('\t')
                interval = tokens[1] +':'+tokens[2]+'-'+tokens[3]
                genes_intervals[interval] = tokens[0]
        elif 'SigPeaks.intersected.bed' in file:#load data from bed files first...
            fH = open (path+file,'r')
            lines = fH.readlines()
            cell_type = file.replace('_as_reference.SigPeaks.intersected.bed','').replace('differential_analysis_res.','')
            if cell_type not in sig_peaks_intervals:
                sig_peaks_intervals[cell_type] = {}
            for line in lines:
                tokens = line.replace('\n','').split('\t')
                sig_peak_interval = tokens[3]+':'+tokens[4]+'-'+tokens[5]
                overlapped_gene_interval = tokens[0]+':'+tokens[1]+'-'+tokens[2]
                sig_peaks_intervals[cell_type][sig_peak_interval] = {'overlapped_gene_interval':overlapped_gene_interval,'log2FoldChange':0}
    #loop again!
    for file in files:
        if 'SigPeaks.csv' in file:
            fH = open (path+file,'r')
            lines = fH.readlines()
            cell_type = file.replace('_as_reference.SigPeaks.csv','').replace('differential_analysis_res.','')
            if cell_type not in sig_peaks_intervals:
                sig_peaks_intervals[cell_type] = {}
            for line in lines:
                tokens = line.split(',')
                parts = tokens[0].split('_')
                interval = parts[0]+':'+parts[1]+'-'+parts[2]
                interval = interval.replace('"','')
                if interval in sig_peaks_intervals[cell_type]:
                    sig_peaks_intervals[cell_type][interval]['log2FoldChange'] = float(tokens[2])
                #    print(sig_peaks_intervals[cell_type])
 

    print('NOTICE: Done Processing analysis folder...')
    return genes_intervals, sig_peaks_intervals

def process_aggregated_folder(path):
    files = os.listdir(path)
    genes_expressions = {}
    for file in files:
        fH = open (path+file,'r')
        lines = fH.readlines()
        del lines[0]
        file_name = file.split('-')
        cell_type = file_name[1].replace('.aggregated.txt','')
        if cell_type not in genes_expressions:
            genes_expressions[cell_type] = {}
        for line in lines:
            tokens = line.replace('\n','').split('\t')
            gene_name = tokens[0].split('.')
            genes_expressions[cell_type][gene_name[0]] = float(tokens[1])
    print('NOTICE: Done Processing aggragated counts folder...')
    return genes_expressions
def analyze_data(genes_intervals, sig_peaks_intervals,genes_expressions):
    genes_expressions_with_less_accessibility = {}
    genes_expressions_with_more_accessibility = {}
    for cell_type in sig_peaks_intervals:
        genes_expressions_with_less_accessibility[cell_type] = {}
        genes_expressions_with_more_accessibility[cell_type] = {}
        for interval in sig_peaks_intervals[cell_type]:
            if sig_peaks_intervals[cell_type][interval]['log2FoldChange'] > 0:
                associated_gene = genes_intervals[sig_peaks_intervals[cell_type][interval]['overlapped_gene_interval']]
                #associated_gene_count =  genes_expressions[genes_intervals[sig_peaks_intervals[cell_type][interval]['overlapped_gene_interval']]]
                associated_gene_count = genes_expressions[cell_type][associated_gene]
                genes_expressions_with_less_accessibility[cell_type][associated_gene] =associated_gene_count
            elif sig_peaks_intervals[cell_type][interval]['log2FoldChange'] < 0:
                associated_gene = genes_intervals[sig_peaks_intervals[cell_type][interval]['overlapped_gene_interval']]
                associated_gene_count = genes_expressions[cell_type][associated_gene]
                genes_expressions_with_more_accessibility[cell_type][associated_gene] =associated_gene_count
    print('NOTICE: Done data intersection...')
    
    print('NOTICE: Mac has {} associated genes with less chromatin accessbility'.format(len(genes_expressions_with_less_accessibility['Mac'])))
    print('NOTICE: Mac has {} associated genes with more chromatin accessbility'.format(len(genes_expressions_with_more_accessibility['Mac'])))
    print()
    print('NOTICE: Neu has {} associated genes with less chromatin accessbility'.format(len(genes_expressions_with_less_accessibility['Neu'])))
    print('NOTICE: Neu has {} associated genes with more chromatin accessbility'.format(len(genes_expressions_with_more_accessibility['Neu'])))
    print()
    print('NOTICE: Mon has {} associated genes with less chromatin accessbility'.format(len(genes_expressions_with_less_accessibility['Mon'])))
    print('NOTICE: Mon has {} associated genes with more chromatin accessbility'.format(len(genes_expressions_with_more_accessibility['Mon'])))

    #check for shared genes...
    print(len(list( genes_expressions_with_more_accessibility['Mon'].keys() & genes_expressions_with_less_accessibility['Neu'].keys())))
    #generate heatmap using rnaseq counts for intersected genes between Mon & Neu with opposite chromatin accessibilty...
    Mon_genes = []
    Neu_genes = []
    gene_names = []
    for gene in genes_expressions_with_more_accessibility['Neu']:
        Mon_genes.append(float( genes_expressions_with_more_accessibility['Neu'][gene]))
    #for gene in genes_expressions_with_less_accessibility['Neu']:
        Neu_genes.append(float( genes_expressions_with_less_accessibility['Mon'][gene]))
        gene_names.append(gene)
    #log transform...
    Mon_genes = np.log2(pd.DataFrame(Mon_genes) + 1)
    Neu_genes = np.log2(pd.DataFrame(Neu_genes) + 1)

    #add gene names...
    Mon_genes.index = gene_names
    Neu_genes.index = gene_names

    combined_genes= pd.concat([Mon_genes, Neu_genes], axis=1)
    plt.figure(figsize=(10, 6))
    x_axis_labels = ['Monocyte','Neutrophil']
    sns.heatmap(combined_genes,annot=False,xticklabels=x_axis_labels)  # Change annot=True if you want to display values
    plt.title('Expression of genes associated different accessbility')
    plt.xlabel('Cell Type')
    plt.ylabel('Genes')
    plt.show()

    return
#
if __name__ == "__main__":
    genes_intervals, sig_peaks_intervals = process_analysis_folder(sys.argv[1])
    genes_expressions = process_aggregated_folder(sys.argv[2])
    analyze_data(genes_intervals, sig_peaks_intervals,genes_expressions)