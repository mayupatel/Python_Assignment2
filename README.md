# Python_Assignment2
Parsing the FASTQ file and accomplishing the tasks



Given a FASTQ file, read in the records and store in an appropriate data structure.

Using this data structure, add functions that will accomplish the following tasks:

 

1. Search record header lines, and return records that match the search criteria.

2. Select a record, and generate a count of each amino acid coded for by the codons of this sequence record.  Keep in mind that because these records are not necessarily in the proper reading frame, so the user should be prompted to select a reading frame (0, +1, +2).  Keep in mind that if you build your data structure of the FASTQ file using a dictionary, that there is no order to the records in dictionaries.  How you plan to present the choice of records to the user during this step may require some cleverness.  Presenting only a subset of the total 25 records might be a good idea.

3. Allow the user to trim the nucleotide sequence record based on the quality score. User should be prompted for a score cutoff, after which a new file ("assignment2_trimmed.fastq") should be generated of the FASTQ records with the sequence and quality score lines trimmed based on that selection.  Two things of note: first, the scoring line is scored using ASCII characters. The score, from 0 - 42 (with 42 being the highest) can be found by using the built in function ord() on the character, then subtracting 64.  For instance:

ord("g") -64 = 39

Secondly, these quality scores are subject to a degree of randomness.  Generally speaking, the first part of the read is high quality, and that tends to drop off the longer the read continues.  It is best to use a sliding window approach, similar to how you identify codons, to find the area at which the score drops off.  Once the average score in your window drops below the threshold, then trim at that location.  How big you make your window is up to you.

These three functions should be called from a main method, that like the previous assignment, collects user inputs and directions, then executes whatever the user selects. The three functions should use the various user inputs and return a result.  The functions should not ask for user input themselves, they should only interpret arguments and execute a task.  


Explanation of FASTQ files:
 
FASTQ files follow this pattern

1. @header
2. "sequence"
3. +header
4. "quality scores"

The header line is repeated twice, before both the sequence and the confidence scores.  It has a different leading symbol though, to help differentiate between the subsequent line.  This predictable 1,2,3,4 pattern makes reading through files relatively easy.  Line 1 of a record will always be a header, line 2 of a record will always be the sequence, etc etc.  Use this to your advantage. 