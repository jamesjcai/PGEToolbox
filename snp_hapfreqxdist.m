function [s,px]=snp_hapfreqxdist(hap,winsize,steplen,pos)
%SNP_HAPFREQXDIST - product of haplotype frequencies and distance
%
%Syntax: s=snp_hapfreqxdist(hap,winsize,steplen)
%        figure; snp_hapfreqxdist(hap,winsize)
%
%See also: NUCDIV

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<3
    steplen=1;
end

if nargin<2
    [numHap,sizHap,seqHap]=counthaplotype(hap);
    p=sizHap./sum(sizHap);
    s=0;
    for i=1:numHap-1
    for j=i+1:numHap
        %d=sum(seqHap(i,:)~=seqHap(j,:))./length(seqHap(j,:));
        d=pdist(seqHap([i j],:),'hamming');
        s=s+2*p(i).*p(j).*d;
    end
    end
    s=s./(1-1/numHap);
    %s=(sum(p.^2)-1/numHap)./(1-1/numHap);   % homozygosity
else
    if isempty(winsize)
        winsize=21;
    end
    m=size(hap,2);
    if m<=winsize
        error('Window size too big')
    end
    s=[];
    px=[];
    for k=1:steplen:m-winsize+1
        hap1=hap(:,k:k+winsize-1);
        p1=mean(pos(k:k+winsize-1));
        s1=snp_hapfreqxdist(hap1);        % s1=nucdiv(hap1);
        s=[s s1];
        px=[px p1];
    end
    if nargout<1
        if pos
            plot(px,s);
        else
            plot(1:length(s),s);
        end
        xlabel('Position')
        ylabel('Haplotype Freq*Dist/(i+j)')
    end
end
