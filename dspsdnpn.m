function [Ds,Ps,Dn,Pn,Ls,Ln] = dspsdnpn(aln)
%DSPSDNPN - calculates Ds, Ps, Dn, and Pn

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

Ds=0;Ps=0;Dn=0;Pn=0;Ls=0;Ln=0;
seq=aln.seq;
if (length(unique(aln.population))==2)
    taxset1=find(aln.population==1);
    taxset2=find(aln.population==2);
else
    [s,v] = choosebox('Name','Interspecies Sequences','PromptString',...
        'Sequences from species 1:','SelectString','Sequences from species 2:',...
        'ListString',aln.seqnames);
        if (v~=1),	return; else
            [n]=size(aln.seq,1);
            taxset1=s;
            taxset2=setdiff(1:n,s);
        end
end
seq1=seq(taxset1,:);
seq2=seq(taxset2,:);

[Ps1,Pn1]=pspn(seq1);
[Ps2,Pn2]=pspn(seq2);
Ps=Ps1+Ps2;
Pn=Pn1+Pn2;
[Ds,Dn]=dsdn(seq1,seq2);


if nargout>4
[Ls,Ln] = lsln(aln);
end
