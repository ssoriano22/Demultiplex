#!/usr/bin/env python

import bioinfo
import gzip
import matplotlib
import argparse
#import itertools


#Function: Argparse
def get_args():
    parser = argparse.ArgumentParser(description="A program to produce a mean QScore distribution for each base position for specified file.")
    parser.add_argument("-IF","--inputfilename",help="File name (FASTQ)",type=str)
    parser.add_argument("-KI","--knownindexfilename",help="Output file name",type=str)
    return parser.parse_args()

#Initializing global variables from argparse - uncomment for final files
# args = get_args()
# f = args.inputfilename #input file
# ki = args.knownindexfilename #output file

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
    return new_header
    
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

#Function: Writing to output files
def write_file(outputfh):
    pass

#Function: Close all output files created in open_outputs() function (52 total). Input: fh dict. Output: None.
def close_outputs(out_fh_dict):
    for key in out_fh_dict:
        out_fh_dict[key][0].close()
        out_fh_dict[key][1].close()
    return None

#MAIN CODE

#4 input data files - use argparse before submit
R1file = "/projects/bgmp/ssoriano/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test1_P2.fq"
R2file = "/projects/bgmp/ssoriano/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test2_P2.fq"
R3file = "/projects/bgmp/ssoriano/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test3_P2.fq"
R4file = "/projects/bgmp/ssoriano/bioinfo/Bi622/Demultiplex/TEST-input_FASTQ/test4_P2.fq"

#Known indexes file - use argparse before submit
knownInd_file = "/projects/bgmp/shared/2017_sequencing/indexes.txt"

#Create set of known indexes from indexes.txt file on Talapas
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
with open(R1file,"r") as fhR1, open(R2file,"r") as fhR2, open(R3file,"r") as fhR3, open(R4file,"r") as fhR4: #Swap with gzip, rt for final files
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
        # if i>1: #FOR TESTING w/ full file
        #     break #FOR TESTING w/ full file
        
        #Make reverse complement of index2
        rc_index2: str = bioinfo.rev_comp(R3_record[1])

        #Edit header lines to add index-pair for current record to R1 and R4 header lines
        R1_record[0] = add_IndexPairHeader(R1_record[0],R2_record[1],rc_index2)
        R4_record[0] = add_IndexPairHeader(R4_record[0],R2_record[1],rc_index2)

        if "N" in R2_record[1] or "N" in R3_record[1]:
            #If Index1 or Index2 contain "N", write R1 + R4 records to unknown R1 + R2 output files
            output_fh_dict["unknown"][0].write(R1_record[0]+"\n"+R1_record[1]+"\n"+R1_record[2]+"\n"+R1_record[3]+"\n")
            output_fh_dict["unknown"][1].write(R4_record[0]+"\n"+R4_record[1]+"\n"+R4_record[2]+"\n"+R4_record[3]+"\n")

        else:

            pass

close_outputs(output_fh_dict)

#Quality = check for Ns in INDEXES (look per base, not ave) to filter records before checking quality scores

#Summarize/explain code output in md file
#Output: % reads from each sample, amount of index swapping, data figures/numbers