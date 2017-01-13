function [aln] = snp_hap2aln(haplodata,hmarkinfo)
%SNP_HAP2ALN - converts HAPLODATA to ALN

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[n,m]=size(haplodata);
[seq] = snp_hap2seq(haplodata,hmarkinfo);

names={};
for (k=1:n),
   names{k}=['Seq_',num2str(k)];
end

m=size(seq,2);
aln.seqnames=names;
aln.seqtype=1;       % non-coding nucleotide
aln.seq=seq;
aln.locus=ones(n,1);
aln.population=ones(n,1);
aln.count=ones(n,1);
aln.pos=[1:m];

