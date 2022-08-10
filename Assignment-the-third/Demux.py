#!/usr/bin/env python

import bioinfo
import gzip
import argparse
import numpy as np
#import itertools


#Function: Argparse
def get_args():
    parser = argparse.ArgumentParser(description="A program to produce a mean QScore distribution for each base position for specified file.")
    parser.add_argument("-R1","--inputfilenameR1",help="Read 1 File name (FASTQ)",type=str)
    parser.add_argument("-R2","--inputfilenameR2",help="Read 2 File name (FASTQ)",type=str)
    parser.add_argument("-R3","--inputfilenameR3",help="Read 3 File name (FASTQ)",type=str)
    parser.add_argument("-R4","--inputfilenameR4",help="Read 4 File name (FASTQ)",type=str)
    parser.add_argument("-KI","--knownindexfilename",help="Output file name",type=str)
    return parser.parse_args()

#Initializing global variables from argparse - uncomment for final files
args = get_args()
R1file = args.inputfilenameR1 #input file R1
R2file = args.inputfilenameR2 #input file R2
R3file = args.inputfilenameR3 #input file R3
R4file = args.inputfilenameR4 #input file R4
knownInd_file = args.knownindexfilename #input file known indexes

#Function: Reverse complement of input sequence - in bioinfo.py v0.4

#Function: Create record of 4 lines (stripping newline chars), returning record as a list
def create_record(filehandle):
    head_line: str = filehandle.readline().strip()
    seq_line: str = filehandle.readline().strip()
    plus_line: str = filehandle.readline().strip()
    qscore_line: str = filehandle.readline().strip()
    record: list = [head_line,seq_line,plus_line,qscore_line]
    #print(record)
    return record

#Function: Add index-pair of current record to header lines for R1 and R4. Inputs: old header line (typ:str), index1, rev. comp. index2. Output: updated header line (typ:string).
def add_IndexPairHeader(old_header,ind1,ind2):
    new_header = old_header + " " + ind1 + "-" + ind2
    ind_pair = ind1 + "-" + ind2
    return new_header,ind_pair
    
#Function: Opening output files - input set of known indexes, return dictionary of file handles (dict lengeth = 26, each key has a list of two file handles, R1 and R2)
def open_outputs(ki_set):
    #Initialize file handle (fh) dictionary
    fh_dict: dict = {}
    for index in ki_set:
        #Dual matched case: For each of the 24 known indexes, open R1 and R2 files for writing (48 files total) and add to fh_dict
        output_filename_R1: str = index + "_R1.fq"
        output_filename_R2: str = index + "_R2.fq"
        fhDMRead1 = open(output_filename_R1,"w")
        fhDMRead2 = open(output_filename_R2,"w")
        fh_dict[index] = [fhDMRead1,fhDMRead2]

    #Index hopped case: open R1 and R2 output files for writing (2 files total) and add to fh_dict
    fhIHRead1 = open("IH_R1.fq","w")
    fhIHRead2 = open("IH_R2.fq","w")
    fh_dict["indexhop"] = [fhIHRead1,fhIHRead2]

    #Unknown case: open R1 and R2 output files for writing (2 files total) and add to fh_dict
    fhURead1 = open("U_R1.fq","w")
    fhURead2 = open("U_R2.fq","w")
    fh_dict["unknown"] = [fhURead1,fhURead2]

    return fh_dict

#Function: Writing to output file. Input: fh, biological record. Output: None.
def write_file(output_fh,output_record):
    output_fh.write(output_record[0]+"\n"+output_record[1]+"\n"+output_record[2]+"\n"+output_record[3]+"\n")
    return None

#Function: Close all output files created in open_outputs() function (52 total). Input: fh dict. Output: None.
def close_outputs(out_fh_dict):
    for key in out_fh_dict:
        out_fh_dict[key][0].close()
        out_fh_dict[key][1].close()
    return None

#Function: Calculate mean quality score (qscore) of given index read. Input: index record. Output: mean qscore of index record (type:int)
def mean_qscore(index_qseq):
    qs_sum: int = 0
    for char in range(0,len(index_qseq)):
        conv_char = bioinfo.convert_phred(index_qseq[char])
        qs_sum += conv_char
    qs_mean: int = qs_sum/len(index_qseq)
    return qs_mean

#MAIN CODE

#4 input data files - use argparse before submit
# R1file = "/projects/bgmp/ssoriano/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test1_P2.fq"
# R2file = "/projects/bgmp/ssoriano/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test2_P2.fq"
# R3file = "/projects/bgmp/ssoriano/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test3_P2.fq"
# R4file = "/projects/bgmp/ssoriano/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test4_P2.fq"

#Known indexes file - use argparse before submit
# knownInd_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"

#Initialize 3 counter variables for total dual matched, index hopped, and unknown records
count_DM: int = 0
count_IH: int = 0
count_U: int = 0

#Initialize 3 dictionaries to count index-pair occurances (key=index-pair, value=count)
dict_DM: dict = {}
dict_IH: dict = {}
#dict_U: dict = {}

#Create set of known indexes from known indexes input file
known_ind_set: set = set()
with open(knownInd_file,"r") as fhKI:
    line_counter: int = 0 #Line counter used to omit header line of indexes.txt file
    for line in fhKI: #Go through each line
        if line_counter != 0: #For every line but the first (header) line
            k_index: str = line.strip().split("\t")[4] #Separate columns by tab char and get 4th column of that line
            known_ind_set.add(k_index) #Add index str to set
        line_counter += 1 #Increment line counter

#Open 52 output files - create dictionary of file handles (fh)
output_fh_dict: dict = open_outputs(known_ind_set)

#Open all four data files and separate records
i: int = 0 #Record counter
with gzip.open(R1file,"rt") as fhR1, gzip.open(R2file,"rt") as fhR2, gzip.open(R3file,"rt") as fhR3, gzip.open(R4file,"rt") as fhR4: #Swap with gzip, rt for final files
    while True:
        i += 1 #Increment record counter (FASTQ starts w/ 1)

        #Use create_record() to create records as LISTS from the current 4 lines of each input file
        R1_record: list = create_record(fhR1)
        R2_record: list = create_record(fhR2)
        R3_record: list = create_record(fhR3)
        R4_record: list = create_record(fhR4)

        if R1_record[0] == "":
            #EOF case - Only need to evaluate one file - all files are same line length
            break

        # if i>1: #FOR TESTING
        #     break #FOR TESTING
        
        #Make reverse complement of index2
        rc_index2: str = bioinfo.rev_comp(R3_record[1])

        #Edit header lines to add index-pair for current record to R1 and R4 header lines
        R1_record[0],index_pair = add_IndexPairHeader(R1_record[0],R2_record[1],rc_index2)
        R4_record[0],index_pair = add_IndexPairHeader(R4_record[0],R2_record[1],rc_index2)

        if "N" in R2_record[1] or "N" in R3_record[1]:
            #Unknown case: If Index1 or Index2 contain "N", write R1 + R4 records to unknown R1 + R2 output files
            count_U += 1 #Increment total U counter
            # if index_pair in dict_U:
            #     #If current index-pair is in dict_IH as key, increment value
            #     dict_U[index_pair] += 1
            # else:
            #     #Else add current index-pair to dict_IH as key, value = 1
            #     dict_U[index_pair] = 1
            write_file(output_fh_dict["unknown"][0],R1_record)
            write_file(output_fh_dict["unknown"][1],R4_record)

        else:
            #Calculate mean quality score of both index sequences
            I1_qs_mean: float = mean_qscore(R2_record[3])
            I2_qs_mean: float = mean_qscore(R3_record[3])

            if I1_qs_mean < 35 or I2_qs_mean < 35:
                #Unknown case: If mean quality score of either index is less than 35 cutoff score - see Demux LB - write R1 + R4 records to unknown R1 + R2 output files
                count_U += 1 #Increment total U counter
                # if index_pair in dict_U:
                #     #If current index-pair is in dict_IH as key, increment value
                #     dict_U[index_pair] += 1
                # else:
                #     #Else add current index-pair to dict_IH as key, value = 1
                #     dict_U[index_pair] = 1
                write_file(output_fh_dict["unknown"][0],R1_record)
                write_file(output_fh_dict["unknown"][1],R4_record)
            
            else:
                #Else indexes meet or exceed qscore cutoff
                if R2_record[1] not in known_ind_set or rc_index2 not in known_ind_set:
                    #Unknown Case: If either index1 or index2 is not in list of known indexes, write R1 + R4 records to unknown R1 + R2 output files
                    count_U += 1 #Increment total U counter
                    # if index_pair in dict_U:
                    #     #If current index-pair is in dict_IH as key, increment value
                    #     dict_U[index_pair] += 1
                    # else:
                    #     #Else add current index-pair to dict_IH as key, value = 1
                    #     dict_U[index_pair] = 1
                    write_file(output_fh_dict["unknown"][0],R1_record)
                    write_file(output_fh_dict["unknown"][1],R4_record)

                elif R2_record[1] != rc_index2:
                    #Index Hopped Case: If index1 is not equal to the rev. comp. of index2, write R1 + R4 records to index hopped R1 + R2 output files
                    count_IH += 1 #Increment total IH counter
                    if index_pair in dict_IH:
                        #If current index-pair is in dict_IH as key, increment value
                        dict_IH[index_pair] += 1
                    else:
                        #Else add current index-pair to dict_IH as key, value = 1
                        dict_IH[index_pair] = 1
                    write_file(output_fh_dict["indexhop"][0],R1_record)
                    write_file(output_fh_dict["indexhop"][1],R4_record)

                elif R2_record[1] == rc_index2:
                    #Dual-Matched Case: If index1 is equal to the rev. comp. of index2, write R1 + R4 records to dual-matched R1 + R2 output files (index1 used as fh_dict key)
                    count_DM += 1 #Increment total DM counter
                    if index_pair in dict_DM:
                        #If current index-pair is in dict_DM as key, increment value
                        dict_DM[index_pair] += 1
                    else:
                        #Else add current index-pair to dict_DM as key, value = 1
                        dict_DM[index_pair] = 1
                    write_file(output_fh_dict[R2_record[1]][0],R1_record)
                    write_file(output_fh_dict[R2_record[1]][1],R4_record)

#FINAL OUTPUTS

#Total counts of dual-matched (DM), index hopped (IH), unknown (U), and total records
count_total: int = count_DM + count_IH + count_U
print("Total Dual-Matched (DM) count: " + str(count_DM))
print("Total Index Hopped (IH) count: " + str(count_IH))
print("Total Unknown (U) count: " + str(count_U))
print("Total records: " + str(count_total))

#Percentages of DM, IH, U records
percent_DM: float = (count_DM/count_total)*100
percent_IH: float = (count_IH/count_total)*100
percent_U: float = (count_U/count_total)*100
print("Percent DM: " + str(percent_DM))
print("Percent IH: " + str(percent_IH))
print("Percent U: " + str(percent_U))
print("")

#Counts of each index pair combo (out of total counts - same as count_total above)
print("Index Pair Results: DM")
for DMpair in dict_DM:
    print(DMpair + "\t" + str(dict_DM[DMpair]) + "/" + str(count_total) + "\t" + str((dict_DM[DMpair]/count_total)*100))
print("")

print("Index Pair Results: IH")
for IHpair in dict_IH:
    print(IHpair + "\t" + str(dict_IH[IHpair]) + "/" + str(count_total) + "\t" + str((dict_IH[IHpair]/count_total)*100))

#DO NOT PRINT UNKNOWNS OR COLLECT IN DICT - save runtime
# print("Index Pair Results: U")
# for Upair in dict_U:
#     print(Upair + "\t" + str(dict_U[Upair]) + "/" + str(count_total) + "\t" + str((dict_U[Upair]/count_total)*100))

#Create distributions of dual matched index-pairs
import matplotlib.pyplot as plt

fig = plt.figure()
#Create positions for x values
x = dict_DM.keys()
plt.xlabel("Index-Pair")

#Set mean qscores to y values
y = dict_DM.values()
plt.ylabel("Number Occurances")

#Set graph title
plt.title("Demux Distribution: Dual Matched Index-Pairs")

#Plot bar graph and save as .png file
plt.bar(x,y)
plt.xticks(rotation=30, ha="right")
plt.savefig("DM_Dist.png")

plt.close()

#Close all output files
close_outputs(output_fh_dict)
