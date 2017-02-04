# PGEToolbox

Population Genetics & Evolution Toolbox
(C) 2007-2017 James Cai

1. Introduction
================
PGEToolbox is a Matlab-based software package for analysis of polymorphism and 
divergence data for population genetics and evolution. It estimates several 
basic statistics of DNA sequence variation and carries out statistical tests 
of selective neutrality under the infinite alleles model, such as Tajima's D 
test, Fu & Li's tests and Fay & Wu's H test. The significance of tests is 
determined from the distribution of the statistics obtained by coalescent 
simulation. The toolbox performs McDonald-Kreitman test (and several extensions). 
PGEToolbox also contains functions for handling SNP (Single Nucleotide 
Polymorphism) genotype data. PGEToolbox is open-sourced, can be easily 
extended or tailored for specific tasks, and scaled up for large data sets.

2. Availability
================
For academic uses, PGEToolbox is available free of charge at
http://bioinformatics.org/pgetoolbox

3. Installation
================
1) Requirements: 
	 i) MATLAB, version 6.0 and higher; 
	ii) MBEToolbox (on which PGEToolbox depends)

2) Download MBEToolbox from http://bioinformatics.org/mbetoolbox.
3) Create a directory called 'mbetoolbox', then unzip the file and copy
   all files there.
4) Adds directory of mbetoolbox to the MATLAB path.
5) Download PGEToolbox from http://bioinformatics.org/pgetoolbox.
6) Create a directory called 'pgetoolbox', then unzip the file and copy
   all files there.
7) Adds directory of pgetoolbox to the MATLAB path.

4. How to start?
==================
1) Start up MATLAB. 
2) Run "PGEDEMO" for demos.


If you have any comments or suggestions then, please, send me an email at
jamescai[at]genomezoo.net
