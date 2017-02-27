function [freqall,nsampl] = sfs(seq,ancseq)
%SFS - Site-Frequency Spectrum
%Generates the frequencies of segregating sites of various types:
%Example: f1=3 f2=2 f3=1 f8=1 f9=1; where, for example, f1=2 means that
%the number of mutations of size 1 (externalmutations) is 2.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<3
    showit=0;
end
if (isstruct(seq)), seq=seq.seq; end
[n,m] = size(seq);

ancidx=0;
if (nargin<2),
    unfolded=0;
    freqall=sfsfolded(seq);
else
    unfolded=1;
    if (length(ancseq(:))==1)      % a reference index of ancestor seq
        ancidx=ancseq;
        idx=setdiff([1:n],ancseq);
        ancseq=seq(ancseq,:);
        seq=seq(idx,:);
    end
    freqall=sfsunfolded(seq,ancseq);
end

nsampl=size(seq,1);

if (nargout<1),
if (unfolded),
    i_dispheader('Unfolded Site-Frequency Spectrum')
    fprintf('Number of input sequences = %d\n',n);
    if ancidx
    fprintf('Sequence %d was used as ancestral sequence\n',ancidx);
    end
    fprintf('Number of processed sequences = %d\n',nsampl);
    fprintf('------------------------------------------\n');
else
    i_dispheader('Folded Site-Frequency Spectrum')
    fprintf('Number of input sequences = %d\n',n);
    fprintf('No sequence was used as ancestral sequence\n');
    fprintf('Number of processed sequences = %d\n',nsampl);
    fprintf('------------------------------------------\n');
end

	for (k=1:length(freqall)),
	      if (freqall(k)>0)
		    fprintf (['f%d=%d '],k,freqall(k));
	      end
    end
 fprintf('\n');
 fprintf('------------------------------------------\n');
  disp('fN=M means that the number of mutations of')
  disp('size N (externalmutations) is M.')

%	for (k=1:length(freq)),
%	      if (freq(k)>0)
%		    fprintf (['%d %d, '],k,freq(k));
%	      end
%	end
%    fprintf ('\n');
i_dispfooter
end
