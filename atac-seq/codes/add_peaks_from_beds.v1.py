import os, sys
from pybedtools import BedTool
#This is add peaks from bed files into one BED file...

def adding_peaks_per_sample (path):
    #load files...
    files = os.listdir(path)
    bed_1 =  BedTool(path+files[0])
    bed_2 =  BedTool(path+files[1])
    bed_3 =  BedTool(path+files[2])

    #sort files...
    bed_1_sorted = bed_1.sort()
    bed_2_sorted = bed_2.sort()
    bed_3_sorted = bed_3.sort()

    #merge files...
    bed_1_merged = bed_1_sorted.merge()
    bed_2_merged = bed_2_sorted.merge()
    bed_3_merged = bed_3_sorted.merge()

    #intersect...
    intersection_1_2 = bed_1_merged.intersect(bed_2_merged, wo=True, f=0.3, r=True)

    intersection_1_2_3 = intersection_1_2.intersect(bed_3_merged,wo=True,f=0.3,r=True)
    intersection_1_2_3_merged = intersection_1_2_3.merge()


    #reduce files...

    return intersection_1_2_3_merged
def read_black_list(black_list_bed):
    blacklist_bed  = BedTool(black_list_bed)
    return blacklist_bed
def conditions_concatenation(c1_bed,c2_bed, c3_bed,black_list_bed):
    print(len(c1_bed));print(len(c2_bed));print(len(c3_bed))
    print()
    reference_bed = c1_bed.cat(c2_bed)
    reference_bed = reference_bed.cat(c3_bed)
    print(len(reference_bed))
    #reduce...
    reference_bed = reference_bed.merge()
    print(len(reference_bed))
    #filter...
    reference_bed = reference_bed.subtract(black_list_bed)
    print(len(reference_bed))
    print(type(reference_bed))
    reference_bed.saveas('reference.bed')
    return
if __name__ == "__main__":
    condition_1_bed = adding_peaks_per_sample(sys.argv[1])
    condition_2_bed = adding_peaks_per_sample(sys.argv[2])
    condition_3_bed = adding_peaks_per_sample(sys.argv[3])
    black_list_bed = read_black_list(sys.argv[4])
    conditions_concatenation(condition_1_bed,condition_2_bed, condition_3_bed,black_list_bed)

