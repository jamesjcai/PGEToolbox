function [h,hvar,coresets]=snp_ehhm(hapldata,n1,n2,showit)
%SNP_EHH - EHH for core with multiple markers between N1 and N2
%
%  Syntax: [h,hvar,coresets]=snp_ehh(hapldata,n1,n2,showit)
%
% n1   -   Position of first core marker
% n2   -   Position of last core marker
% H        - EHH value
% HVAR     - variance of EHH
% CORESETS - core haplotype set (used to be determined ancestral state)
%
%Haplotype homozygosity (HH) is an effective measure of linkage 
%disequilibrium (LD) for more than 2 markers. If we want to see, how LD 
%breaks down with increasing distance to a specified core region, we can 
%calculate HH in a stepwise manner as EHH (Sabeti et al. 2002). EHH is 
%calculated between a distance x and the specified core region for a 
%chromosome population carrying a single core haplotype. Distance x 
%increases stepwise to the most outlying marker. The procedure is repeated 
%for each core haplotype. HH is evaluated as 
%
%            HH = [sum(pi^2) - 1/n] / [1 - 1/n] 
%
%with pi being the relative haplotype frequency and n the sample size. 
%It corrects for sampling effects (Sabatti & Risch 2002). The variance 
%of HH is estimated according to Nei (1975).
%
% geno, mark, phas
% genodata, markinfo, phas
%
%EHH estimates the level of haplotype splitting due to recombination and 
%mutation at extended regions on both sides of a specified core region. I
%t may be used to explore haplotype-specific LD patterns, e.g. for disease 
%associated haplotypes. In combination with the core haplotype frequency it 
%may also serve as an indicator of recent positive selection. Frequent core 
%haplotypes with an unusually high long-range LD are supposed to be 
%positively selected. The various core haplotypes can serve as internal 
%controls. HH=(nchoosek(p*n,2)+nchoosek((1-p)*n,2))/nchoosek(n,2);

%if nargin<3, n2=n1; end
%if (n2<n1), error('n2 must be greater than n1'); end

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (nargin<4), showit=false; end
if (nargin<1)
    disp('Using example data...')
    % using example data
    olddir=pwd;
    cdpge;
    %hapldata=load('example_data/ehhexample')
    load('example_data/ehhexamplemat');
    cd(olddir);
end
if (nargin<2)
      nmarker=size(hapldata,2);
      n1=floor(nmarker/2);
end
if (nargin<3), n2=n1; end

if ~(n1>0 && n2>0 && n2>=n1)
    error('Incorrect n1 or/and n2');
end


[maffreq,allemajor,alleminor]=snp_maf(hapldata,1);     % MAF of SNPs


[core]=hapldata(:,n1:n2);
[coresets,temp,coretype]=unique(core,'rows');

% sort coresets and coretype according to frequencies (descending order)
[temp,idx]=sort(grpstats(coretype,coretype,'length'));
idx=idx(end:-1:1);
coresets=coresets(idx);
newcoretype=zeros(size(coretype));
for k=1:length(idx)
    newcoretype(coretype==idx(k))=k;
end
n=size(coresets,1);


data_left=hapldata(:,1:n1-1);
data_right=hapldata(:,n2+1:end);

[temp,m]=size(hapldata);
m2=m-(n2-n1);   % n2-n1 is collasped into core

h=ones(n,m2);     % initialize result H. In the case of K = 1, let EHH = 1.
hvar=zeros(n,m2); % initialize result var(H).

core_p=zeros(n,1);


        
if (showit), wbar = waitbar(0,'Please wait...'); end

for (k=1:n),
      idx=(newcoretype==k);
      core_p(k)=mean(idx);
      dl=data_left(idx,:);
      dr=data_right(idx,:);
      for (i=1:n1-1),
        if (showit), waitbar(i/m2,wbar); end
	    [temp1,y,temp2]=counthaplotype(dl(:,i:end));
        z=sum(y); 
        if z>1
            p=y./z; a=sum(p.^2);          % p=relative haplotype frequency
            h(k,i)=(a-1/z)/(1-1/z);
            if nargin>1
            hvar(k,i)=(2*(z-1)/z^3)*(2*(z-2)*(sum(p.^3)-a^2) + a - a^2 );
            end
        end        
      end
      for (j=n2+1:m2),        
	    [temp1,y,temp2]=counthaplotype(dr(:,1:j-n1));
	    z=sum(y);
        if z>1
            p=y./z; a=sum(p.^2);
            h(k,j)=(a-1/z)/(1-1/z);
            if nargin>1
            hvar(k,j)=(2*(z-1)/z^3)*(2*(z-2)*(sum(p.^3)-a^2) + a - a^2 );
            end
        end
        if (showit), waitbar(j/m2,wbar); end
      end
end
if (showit), close(wbar); end


if (nargout<1 || showit)
    figure;    
    hold on
    for (k=1:n-1)
        plot(h(k,:),'-')
        corehaplotxt{k}=sprintf('Allele: %s; Freq.: %.3f',...
                int2nt(coresets(k,:)),core_p(k));
    end
    plot(h(end,:),'-m')
        corehaplotxt{n}=sprintf('Allele: %s; Freq.: %.3f',...
                int2nt(coresets(n,:)),core_p(n));    
    vline(n1+floor((n2-n1)/2))
    hline(0.05);
    hold off    
    %legend(corehaplotxt,'location','NorthOutside')
    legend(corehaplotxt)
    %legend(num2str([1:3]),-1)
    ylabel('EHH')
    xlabel('No. of extended SNPs')
end

