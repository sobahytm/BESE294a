import os, sys
#This is to filter peaks based on lof & padj ...
#and write the peaks to a bed file...
try:
    if len(sys.argv) < 2:
        print('USAGE: python differential_analysis_res.csv')
        exit(0)
except:
    pass
def process_DEx_atact_data(de_file):
    fH = open (de_file,'r')
    lines = fH.readlines()
    del lines[0]
    output= []
    output.append(lines[0]);del lines[0]
    for line in lines:
        tokens = line.replace('\n','').split(',')
        try:
            if float(tokens[-1]) < 0.01:
                if float(tokens[2]) > 1:
                    output.append(line)
                elif float(tokens[2]) < -1:
                    output.append(line)
        except Exception as error:
            print(error)
    print('NOTICE: Found {} significant peaks.'.format(len(output)-1))
    return output

def writing(file_name,output):
    sig_peaks_file = file_name.replace('.csv','.SigPeaks.csv')
    bed_file = file_name.replace('.csv','.SigPeaks.bed')
    fOut1 = open (sig_peaks_file,'w')
    fOut2 = open (bed_file,'w')
    for line in output:
        tokens = line.split(',')
        if tokens[0] == 'id':
            fOut1.write(line)
        else:
            fOut1.write(line)
            interval = tokens[0].replace('"','').replace('_','\t') + '\n'
            fOut2.write(interval)
    return
if __name__ == "__main__":
    output = process_DEx_atact_data(sys.argv[1])
    writing(sys.argv[1],output)