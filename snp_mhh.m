function [mhh,hh]=snp_mhh(genodata)
%SNP_MHH - MHH, most frequent haploytpe homozygosity 
% Syntax: [mhh,hh]=snp_mhh(genodata)
%
%MHH - Most frequent haploytpe homozygosity 
%HH  - haplotype homozygosity 
 	
%Kimura R, Fujimoto A, Tokunaga K, Ohashi J.	
%A practical genome scan for population-specific strong selective sweeps that have reached fixation.
%PLoS ONE. 2007 Mar 14;2:e286.
%PMID: 17356696

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


G=snp_hhgeno(genodata);
%G=G(:,[4:10]);
[n,m]=size(G);

%if m>4    

[~,b,c]=counthaplotype(G);
mhh=b(sum(c,2)==m)/n;

G(G==2)=1;
[~,b,c]=counthaplotype(G);
hh=b(sum(c,2)==m)/n;

if isempty(mhh), mhh=0; end
if isempty(hh), hh=0; end

fprintf('MHH = %f\t',mhh);
fprintf(' HH = %f\n',hh);

