function [t]=time2abs(p,N)
% [t]=time2abs(p,N)

twosided=true;
[t]=i_time2prob(p,N,twosided);

function [t]=i_time2prob(p,N,twosided)
%
% REF: Kimura Ohta 1973 (equ 13 and page 10)
% http://www.pubmedcentral.nih.gov/picrender.fcgi?artid=1212997&blobtype=pdf
% Applications e.g., Stephens JC (1998) CCR5-delta32
% Book: Mathematical Population Genetics, WJ Ewens, (equ. 3.4 and 3.11)
% 
% See also: TIME1ST2PROB

if nargin<3, twosided=false; end
if nargin<2, N=10000; end

if twosided
    t=-4.*N.*(p.*log(p)+(1-p).*log(1-p));
else
    %t=-4.*N.*(p.*log(p)./(1-p));             % (page 10)
    t=-4.*N.*(1-p).*log(1-p)./p;              % (3.11), i.e., time to fix
end

%{
x=0:0.01:1;
y1=time2prob(x,10000,true);
y2=time2prob(x,10000,false);
plot(x,y1,'r-'); hold on; plot(x,y2,'b-')
gtext('-4.*N.*(p.*log(p)+(1-p).*log(1-p))')
gtext('-4.*N.*(1-p).*log(1-p)./p')
gtext('t(1/2)~=2.8N')
xlabel('p'); ylabel('Generations')
title('Time to absorption (red) or fixation (blue)')
%}
