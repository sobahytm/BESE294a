import os, sys

try:
    if len(sys.argv) <2:
        print('USAGE: python count_matrix.csv')
        exit(0)
except:
    pass

def read_counts_table(counts_file):
    fH = open(counts_file,'r')
    lines = fH.readlines()
    header  =  lines[0]
    del lines[0]
    negeative_counts =[];clean_counts = []
    for line in lines:
        tokens = line.replace('\n','').split(',')
        check =0
        for i in range(1,len(tokens),++1):
            try:
                if float(tokens[i]) < 0:
                    negeative_counts.append(line)
                    check += 1
            except  Exception as error:
                print(error)
                check +=1
        if check ==0:
            clean_counts.append(line)
                #print(line)
    if negeative_counts == 0:
        print('NOTICE: No negative counts!')
    else:
        print('WARNING: Negative counts were found!')
    fOut = open (counts_file.replace('.csv','.filtered.csv'),'w')
    fOut.write(header)
    for line in clean_counts:
        fOut.write(line)

    return

if __name__ == "__main__":
    read_counts_table(sys.argv[1])