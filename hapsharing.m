function [hss]=hapsharing(s,pos)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

%s=h1==h2;
if nargin<2
%   hss=zeros(1,length(h1));
    hss=zeros(size(s));
    st=sprintf('%d',s);    
    [d1]=i_runscount(st,'0');
    [d0]=i_runscount(st,'1');
    leadingzerosn=0;
    if ~s(1)
        leadingzerosn=d0(1);
        d0=d0(2:end);
        hss=hss(leadingzerosn+1:end);
    end
    cd1=cumsum(d1);
    cd0=cumsum(d0);
    
    % need moving leading zeros    
    if length(cd1)>length(cd0)
        idx=[0 cd1(1:end-1)+cd0]+1;
        for k=1:length(idx)
            hss(idx(k):idx(k)+d1(k)-1)=d1(k);
        end
    else       % equal length        
        idx=[0 cd1+cd0]+1;
        for k=1:length(idx)-1
            hss(idx(k):idx(k)+d1(k)-1)=d1(k);
        end
    end
    if leadingzerosn>0
        hss=[zeros(1,leadingzerosn),hss];
    end
else
    hss=0;
    if s(pos)
        s1=s(1:pos-1);
        s1=s1(end:-1:1);
        s2=s(pos+1:end);
        [n1]=i_getn(s1);
        [n2]=i_getn(s2);
        hss=n1+n2+1;
        hss=hss-(hss==1);
    end
end



function [n1]=i_getn(s1)
   n1=0;
   if ~isempty(s1) && s1(1)
       s1=sprintf('%d',s1);
       t1=textscan(s1,'%s','delimiter','0','multipleDelimsAsOne',true);
       t1=t1{:}; n1=length(t1{1});
   end

function [d1]=i_runscount(st,delim)    
    t1=textscan(st,'%s','delimiter',delim,'multipleDelimsAsOne',true);
    t1=t1{:};
    d1=zeros(1,length(t1));
    for k = 1:length(t1)
        d1(k) = length(t1{k});
    end
