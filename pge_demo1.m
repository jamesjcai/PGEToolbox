%%PGEToolbox DEMO - Neutrality Tests
% Welcome to PGEToolbox.  This is a demonstration of
% PGEToolbox's functions for neutrality tests.
%
% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


cd(fileparts(which('pge_demo1')))
% Current working directory has been changed.

%%
% First, let's just show how to read a FASTA formatted file
% into a MATLAB structure, ALN, representing the alignment.
% Click OK to assign ALN as coding sequences.

filename=fullfile('seq_examples','mktest.fas');
formatid=1; % 1 - FASTA format; 2 - Phylip format
[aln] = pge_openfile(filename,formatid)


%%
% Now let's view the sequences information in this alignment.
% Information includes name, locus, population, count and sequence.

aln.population([1:6])=ones(1,6);
aln.population([7:11])=ones(1,5)*2;

pge_viewdata(aln)

%%
% Extract sequences for one of populations.

alnHsa=aln;
alnHsa.seq=aln.seq([1:6],:);
alnHsa.seqnames=aln.seqnames([1:6]);
alnHsa.locus=aln.locus([1:6]);
alnHsa.population=aln.population([1:6]);
alnHsa.count=aln.count([1:6]);

%%
% REPORTPOLYSITES reports polymorphic sites of sequences

reportpolysites(alnHsa);

%%
% Nucleotide diversity in a population

nucdiv(alnHsa);

%%
% Estimates theta_pi and theta_w

estimatetheta(alnHsa);

%%
% Tajima (1989) proposed to use the two different estimates theta_pi
% and theta_w to detect selection.

tajima89d_test(alnHsa);

%%
% Statistical Tests of Neutrality of Mutations --
% Fu and Li's D* and F* tests

fuli93dsfs_test(alnHsa);

%%
% R2 test, Fu (97) F test and Wall's B and Q tests

r2_test(alnHsa);
Fu97Fs(alnHsa);
wall99BQ(alnHsa);

%%
% Command-line McDonald-Kreitman test. You may have to assign sequences into
% two populations.

mktestcmd(aln);

%%
% GUI for McDonald-Kreitman test. You may have to assign sequences into
% two populations.

mktestgui(aln);

%%
% Coalescent simulations dialog

coalsimdlg

%%
% Sliding window analysis for Tajima's D statisitic

slidingfun(@tajima89d_test,alnHsa,30,5);

%%
% Thank you for viewing this introduction to PGEToolbox neutrality test funcations.