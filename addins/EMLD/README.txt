A quick and handy program to calculate pairwise linkage disequilibrium
It takes in genotype data from unrelated individuals and uses EM algorithm 
to estimate pair-wise haplotype frequencies, then calculates LD statistics.

Input file: genotype data(save as "file-name.dat")
First line is the number of individuals and the number of SNP markers.
Start from the second line: ID, 1st alle of 1st marker, 2nd allele of 1st marker, 
1st alle of 2nd marker, 2nd allele of 2nd marker,.....delimited by tab or space.
Alleles are coded as 1's and 2's, missing values are coded as 0.
Example:
230	7
300203	2	2	1	2	
300401	1	2	1	1	
300601	0	0	1	2	
300802	2	2	1	1	
301001	2	2	2	2	
301103	1	1	1	1	
301201	1	2	1	1	
301304	1	2	0	0	
301503	2	2	1	1	
301601	1	1	1	1	
......
......

Output files:
1) LD.xt: a file in the format of input file for GOLD (http://www.sph.umich.edu/csg/abecasis/GOLD/)
with pair-wise LD(D, D' and r2). So you can plug it in directly and get LD plot.
2) HapFreq.txt: Pair-wise haplotype frequecies estimated by EM algorithm.

Source Code: java
Usage:
If you have JDK1.4.2 or above installed and have the path correctly settled, 
a simple click "EMLD.bat" will give you all the results. 
Otherwise:
1) Download: JDK1.4.2 (http://java.sun.com/j2se/1.4.2/download.html)
2) Install it in your computer and set the path as said in the installation instruction.
3) click "EMLD.bat"

In case it does not work, compile the program yourself and run it.
>javac *.java
>java EMLD inputfile