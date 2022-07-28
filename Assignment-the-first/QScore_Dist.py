#!/usr/bin/env python

#Use code from PS9 (PS4 but simplified w/ arrays) and argparse

#Import modules
import bioinfo
import gzip
import numpy as np
import argparse

#Define arguments to use on command line
def get_args():
    parser = argparse.ArgumentParser(description="A program to produce a mean QScore distribution for each base position for specified file.")
    parser.add_argument("-f","--filename",help="File name (FASTQ)",type=str)
    parser.add_argument("-l","--readlength",help="Read Length for file",type=int)
    parser.add_argument("-o","--output",help="Output file name",type=str)
    return parser.parse_args()

#Initializing global variables from argparse
args = get_args()
f = args.filename #filename - "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" 
                  #filename - "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"             
                  #filename - "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" 
                  #filename - "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" 
                  #filename - test - "/projects/bgmp/ssoriano/bioinfo/Bi622/Demultiplex/Assignment-the-first/test_P1.fq"
l = args.readlength
o = args.output

#Populate numpy array of arrays
qs_array = np.zeros((l),dtype=float)
#print(qs_array)

#Open file and fill array
with gzip.open(f,"rt") as fh:
    line_count = 0
    num_records = 0
    for line in fh:
        strip_line = line.strip()
        if line_count % 4 == 3:
            for char in range(0,len(strip_line)):
                conv_char = bioinfo.convert_phred(strip_line[char])
                qs_array[char] += conv_char
            num_records += 1
        line_count += 1   
#print(qs_array)

#Make new array containing the mean of each inner array in qs_array
mean = np.zeros((l),dtype=float)
i = 0
for score in qs_array:
    print("Mean calc base number: "+str(i))
    mean[i] = score/num_records
    i += 1
print(mean)

import matplotlib.pyplot as plt

#Create positions for x values
x = np.arange(l)
plt.xlabel("Position")

#Set mean qscores to y values
y = mean
plt.ylabel("QScore Mean")

#Set graph title
plt.title("QScore Distribution: Mean QScore per Pos")

#Plot bar graph and save as .png file
plt.bar(x,y)
plt.savefig(o)