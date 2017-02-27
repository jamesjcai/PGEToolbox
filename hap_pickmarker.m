function [haplotype2,markinfo2]=hap_pickmarker(haplotype,markinfo,s)
%[haplotype2,markinfo2]=hap_pickmarker(haplotype,markinfo,idx)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

haplotype2=haplotype;
markinfo2=[];
if ~isempty(markinfo)
    markinfo2=markinfo;
end

if nargin<3    
[s,v] = choosebox('Name','Pick including marker(s)','PromptString',...
    'Markers available:','SelectString','Selected markers:',...
    'ListString',markinfo.rsid'); 
else
    v=1;
end
if (v==1),
    if islogical(s)
        s=find(s);
    end    
    haplotype2 = haplotype(:,s); 
if ~isempty(markinfo)    
    a=fieldnames(markinfo);
    for k=1:length(a)
        if strcmp(a{k},'popid')
            markinfo2.popid=markinfo.popid;
        else
            markinfo2=i_saftmapping(markinfo2,markinfo,a{k},s);
        end
    end
end

end

function markinfo2=i_saftmapping(markinfo2,markinfo1,fieldtxt,s)
    if isfield(markinfo2,fieldtxt)&&isfield(markinfo1,fieldtxt)
        v=getfield(markinfo1,fieldtxt,{s});        
        %markinfo2=setfield(markinfo2,fieldtxt,v);
        markinfo2.(fieldtxt)=v;
    end

