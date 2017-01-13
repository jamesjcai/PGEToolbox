function [ihs,majalle,minalle]=snp_ihs(hapldata,hmarkinfo,coreidx)
%SNP_IHS - Integrated haplotype score.
%
%  Syntax: [ihs]=snp_ihs(hapldata,hmarkinfo,coreidx)
% IHS is normalized score

% if nargin<3, regsize=1.0; end
% if nargin<2, popid='CEU'; end
% if nargin<1, rsid='rs7284767'; end
% if strcmp(popid,'HCB'), popid='CHB'; end
% popid=upper(popid);
% if ~(ismember(popid,{'CEU','CHB','JPT','YRI','JPT+CHB'})), 
%   error('Not available population code.');	 
% end

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin<2)
    error('SNP_IHS needs at least two inputs.');
end
if (nargin<3)
      nmarker=size(hapldata,2);
      coreidx=floor(nmarker/2);
end
[maffreq,majalle,minalle]=snp_maf(hapldata,1);     % MAF of SNPs
majalle=majalle(coreidx);
minalle=minalle(coreidx);


cutoff=0.015;
if maffreq(coreidx)<cutoff
    ihs=nan;
    warn('Core SNP has a low MAF (<0.015). iHS cannot be computed.');
    return;
end

idx=find(maffreq>=cutoff);

if (isempty(idx))
      error('No valid SNPs (MAF>=0.015)');
end
if (length(idx)<3)
      error('Too few valid SNPs (MAF>=0.015)');
end


xpos=hmarkinfo.pos(idx);
%xmaf=hmarkinfo.maf(idx);
xrsid=hmarkinfo.rsid(idx);
xk=find(ismember(xrsid,hmarkinfo.rsid{coreidx}));
xhapldata=hapldata(:,idx);
%xk

	ehh=snp_ehh(xhapldata,xk);
	if size(ehh,1)==2
    ihh1=trapz(xpos,ehh(1,:));
    ihh2=trapz(xpos,ehh(2,:));
        if ihh2>0
            ihs=log(ihh1/ihh2);
        else
            ihs=nan;
        end
    end
        


%{
disp('iHS scores were standardized empirically with the CEU HapMap data (mean = 0.7970; s.d.= 0.6040)');
figure;
plot(pos,ihs,'.');
ylabel('|iHS|')
xlabel('Genomic Postion')
title('|iHS| Scatter Plot')
hline(2.0)
%}

%iHS scores can be standardized using estimates of the mean and s.d. 
%obtained via coalescent simulation under a variety of demographic models.
%These simulations were tailored to match the frequency spectrum, SNP 
%density and recombination profile of the observed data. 
%Alternative demographic models included either exponential growth or a 
%bottleneck (which varied in onset, severity, duration and population size 
%recovery after the bottleneck). 

if (nargout<1)
i_dispheader('Integrated Haplotype Score (iHS)')
        fprintf('%s\t(%s/%s %d)\t%f\n', hmarkinfo.rsid{coreidx},...
            int2nt(majalle),int2nt(minalle),hmarkinfo.pos(coreidx),ihs);
i_dispfooter    
end

