function [binvaluev,rawbinv]=histsfs(freqvv,itype,cutoffval)
%HISTSFS - Histogram of Site-Frequency Spectrum
% histsfs(freqvv,itype)
%
% itype = 1 folded
% itype = 2 unfolded

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-12-27 00:25:00 -0600 (Fri, 27 Dec 2013) $
% $LastChangedRevision: 755 $
% $LastChangedBy: jcai $

if ~iscell(freqvv)
    freqvv={freqvv};
end

if nargin<3
    cutoffval=0.2;
end

binvaluev=[];
rawbinv=[];

for k=1:length(freqvv)
    freqv=freqvv{k};
    %freqv=freqv(freqv>0 & freqv<1);
    %freqv=freqv(freqv<1);
    %freqv=[freqv(freqv<0.5);1-freqv(freqv>=0.5)];
    if nargin<2, itype=1; end


    switch itype
        case 1
       binvalue=histc(freqv,0:0.05:0.5);
    %   txt={'0-0.05','0.05-0.1','0.1-0.15','0.15-0.2',...
    %    '0.2-0.25','0.25-0.3','0.3-0.35','0.35-0.4','0.4-0.45','0.45-0.5'};
       txt={'0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50'};
       xtxt='Mininum Allele Frequence';
        case 2
       binvalue=histc(freqv,0:0.1:1);
    %   txt={'0-0.1','0.1-0.2','0.2-0.3','0.3-0.4',...
    %    '0.4-0.5','0.5-0.6','0.6-0.7','0.7-0.8','0.8-0.9','0.9-1.0'};
       txt={'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'};
       xtxt='Derived Allele Frequency';
        case 3
    %   binvalue=histc(freqv,[0 eps 0.1 1]);
    %   txt={'0','<0.1','>=0.1'};
       binvalue=histc(freqv,[0 cutoffval 1]);
       txt='bb';
    %   txt={sprintf('<%.2f',cutoffval(1)),...
    %       sprintf('%.2f<= & <%.2f',cutoffval(1),cutoffval(2)),...
    %       sprintf('>=%.2f',cutoffval(2))};

       txt={sprintf('<%g',cutoffval(1)),...
           sprintf('>=%g',cutoffval(1))};

       xtxt='Derived Allele Frequency';
        otherwise
            error('Wrong switch')
    end

    [n,m]=size(binvalue);
    if m>n
        binvalue=binvalue';
    end

    binvalue(end-1)=binvalue(end-1)+binvalue(end);
    binvalue(end)=[];

    rawbinv=[rawbinv,binvalue];
    binvalue=binvalue./sum(binvalue);
    binvaluev=[binvaluev, binvalue];
end


if nargout<1
    bar(binvaluev,1)
    set(gca,'XTickLabel',txt)
    xlim([0 length(binvalue)+1])
    xlabel(xtxt)
    ylabel('Fraction of SNPs')
end

% See also: BINOCI
% Error bars denote 95% confidence intervals for binomial expectations.