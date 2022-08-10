# Demultiplexing (Demux) Summary 
### Sophia Soriano

## Demux Code Output: Read Count Data

```
Total Dual-Matched (DM) count: 239750646
Total Index Hopped (IH) count: 352378
Total Unknown (U) count: 123143711
Total records: 363246735
Percent DM: 66.00214754855264
Percent IH: 0.09700789189474751
Percent U: 33.90084455955261
```

Summary: Approximately 2/3 of all reads had correctly dual-matched indexes, while approximately 1/3 were either low quality or had an unknown index sequence. Only a very small percentage (0.1%) of reads had hopped indexes, which is a positive sign that the correct adapter concentration was likely used in library prep. The percentage of unknown reads seems a little high, so the sequencing run that produced this data seems to have  not produced very high quality index sequences (<35 qscore cutoff or unknown index, likely caused by a base calling error).

## Demux Code Output: Dual-Matched Reads 
Note: See "slurm-21916045.out" for complete dual-matched and index hopped Demux.py outputs.

Bash command line used to sort Demux slurm output (DM) by percentage column:

```
$ head -33 slurm-21916045.out | tail -24 | sort -rnk 3
```

Bash-sorted Demux Output (DM reads):

```
TACCGGAT-TACCGGAT       51581889/363246735      14.200234724752583
TCTTCGAC-TCTTCGAC       31203016/363246735      8.590033438290918
CTCTGGAT-CTCTGGAT       25933763/363246735      7.139434577436739
CTAGCTCA-CTAGCTCA       13661286/363246735      3.760883356597823
TGTTCCGT-TGTTCCGT       12247310/363246735      3.3716228722606414
AGAGTCCA-AGAGTCCA       8234276/363246735       2.2668547867333206
TAGCCATG-TAGCCATG       7951357/363246735       2.188968608348262
TATGGCAC-TATGGCAC       7880872/363246735       2.169564442196569
TCGAGAGT-TCGAGAGT       7791551/363246735       2.1449748199388496
ATCATGCG-ATCATGCG       7448228/363246735       2.0504597240220206
GTCCTAAG-GTCCTAAG       6658709/363246735       1.8331091124604328
AACAGCGA-AACAGCGA       6606498/363246735       1.8187356866401014
AGGATAGC-AGGATAGC       6571536/363246735       1.8091108238041012
ACGATCAG-ACGATCAG       6234505/363246735       1.7163278838555835
GTAGCGTA-GTAGCGTA       6033486/363246735       1.66098836373574
ATCGTGGT-ATCGTGGT       5031584/363246735       1.3851697799843954
GATCAAGG-GATCAAGG       4900440/363246735       1.3490664960828898
GCTACTCT-GCTACTCT       4720342/363246735       1.29948642208718
CGATCGAT-CGATCGAT       4398491/363246735       1.2108824598244496
TCGGATTC-TCGGATTC       3149983/363246735       0.8671744840321827
CGGTAATC-CGGTAATC       3097661/363246735       0.8527705004698803
GATCTTGC-GATCTTGC       2861868/363246735       0.7878578729689063
TCGACAAG-TCGACAAG       2786134/363246735       0.7670086835054416
CACTTCAC-CACTTCAC       2765861/363246735       0.761427628523626
```

DM Index Distribution:

![Histogram of DM Indexes](/projects/bgmp/ssoriano/bioinfo/Bi622/Demultiplex/Assignment-the-third/DM_Dist.png)

Summary: The TACCGGAT index accounted for the highest percentage of dual-matched indexes (14.2%), followed by the TCTTCGAC index (8.6%) and the CTCTGGAT index (7.1%). 

## Demux Code Output: Index Hopped Reads
Note: See "slurm-21916045.out" for complete dual-matched and index hopped Demux.py outputs.

Bash command line used to sort Demux slurm output (IH) by fraction column:

```
$ head -587  slurm-21916045.out | tail -552 | sort -rnk 2
```

Bash-sorted Demux Output (IH reads): Ten highest percentage index hopped index-pairs shown here

```
TATGGCAC-TGTTCCGT       63153/363246735 0.017385703411759504
TGTTCCGT-TATGGCAC       61181/363246735 0.016842821725568984
CTAGCTCA-TCGACAAG       9926/363246735  0.0027325778991516606
GATCAAGG-TCTTCGAC       9552/363246735  0.0026296175793569074
TCGACAAG-ATCATGCG       5896/363246735  0.001623139159117287
CTCTGGAT-TACCGGAT       5821/363246735  0.0016024920361637936
TACCGGAT-CTCTGGAT       5796/363246735  0.0015956096618459626
TACCGGAT-TCTTCGAC       5626/363246735  0.001548809516484711
GTCCTAAG-TATGGCAC       4869/363246735  0.0013404112221407853
TCTTCGAC-TACCGGAT       4550/363246735  0.0012525921258452605
...
```

Summary: Two index hopped index-pairs were more frequent than any other observed swapped index-pairs: TATGGCAC-TGTTCCGT (63153 occurances, 0.017%) and TGTTCCGT-TATGGCAC (61181 occurances, 0.017%). However, given the total number of reads (363246735), even these highest occurances of index hopping are very small percentages of the total number of reads.