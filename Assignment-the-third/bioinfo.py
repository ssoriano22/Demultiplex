#!/usr/bin/env python

# Author: Sophia Soriano ssoriano@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This bioinfo.py module contians several useful bioinformatics functions including:
convert_phred, qual_score, validate_DNA_seq, validate_base_seq, gc_content, and
oneline_fasta. Use --help on any of these functions to learn more about each. '''

__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = ""
RNA_bases = ""

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    phred_num: int = (ord(letter)) - 33
    return phred_num

def qual_score(phred_score: str) -> float:
    """Takes a phred score string and returns a single average quality score (float)
    for the whole sequence."""
    score_sum = 0
    for char in phred_score:
        num_score = convert_phred(char)
        score_sum += num_score
    return score_sum/len(phred_score)

def validate_DNA_seq(DNA):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts, Gs, and Cs. False otherwise. Case insensitive.'''
    DNA = DNA.upper()
    return len(DNA) == DNA.count("A") + DNA.count("T") + DNA.count("G") + DNA.count("C") + DNA.count("N")

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    seq = seq.upper()
    return len(seq) == seq.count("A") + seq.count("U" if RNAflag else "T") + seq.count("G") + seq.count("C") + seq.count("N")

def gc_content(DNA):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_DNA_seq(DNA), "String contains invalid characters"
    DNA = DNA.upper()         #Make sure sequence is all uppercase
    Gs = DNA.count("G")       #count the number of Gs
    Cs = DNA.count("C")       #count the number of Cs
    return (Gs+Cs)/len(DNA)

def oneline_fasta(inputfile,outputfile):
    '''Given input FASTA file and output file name, returns new FASTA file with 2-line
    per record (1 header line, 1 seq line) format'''
    #Open input fasta file and intermediate output fasta file (seq on one line) and separate records
    with open(inputfile,"r") as fh:
        with open(outputfile,"w") as fh2:
            #For each record in fasta file
            for i,line in enumerate(fh):
                line = line.strip()
                if line.startswith(">"):
                    if i == 0:
                        #For first line
                        #print(line)
                        fh2.write(f"{line}\n")
                    else:
                        #For all other header lines besides first line in file
                        #print(f"\n{line}")
                        fh2.write(f"\n{line}\n")
                else:
                    #For seq lines
                    #print(line, end="")
                    fh2.write(f"{line}")
            #For end of file
            #print("")
            fh2.write("\n")

def rev_comp(inputDNA):
    '''Given an input DNA index sequence, returns the reverse complement of that sequence'''
    rev_compDNA = ""
    #Reverse input DNA string
    rev_inputDNA = inputDNA[::-1]
    for base in rev_inputDNA:
        #For each base, add the complementary base to new rev comp string
        if base == "A":
            rev_compDNA += "T"
        elif base == "C":
            rev_compDNA += "G"
        elif base == "T":
            rev_compDNA += "A"
        elif base == "G":
            rev_compDNA += "C"
        elif base == "N":
            #Unknown "N" base is still unknown
            rev_compDNA += "N"
    #print(rev_compDNA)
    return rev_compDNA
    

if __name__ == "__main__":
    # write tests for functions above
    # Tests for convert_phred()
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    # Tests for qual_score()
    phred_score: str = "FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@"
    assert qual_score(phred_score) == 37.62105263157895, "wrong average phred score"
    print("You calcluated the correct average phred score")
    # Tests for gc_content
    assert gc_content("GCGCGC") == 1, "messed up calc when all GC"
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("correctly calculated GC content")
    # Tests for validate_DNA_seq()
    assert validate_DNA_seq("aaaaa") == True, "DNA string not recognized"
    print("Correctly identified a DNA string")
    assert validate_DNA_seq("Hi there!") == False, "Non-DNA identified as DNA"
    print("Correctly determined non-DNA")
    # Tests for validate_base_seq()
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")