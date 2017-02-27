function [C]=readstructureoutput(filename)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

disp(['Reading STRUCTURE output file ',filename,' ...'])
txt = textread(filename,'%s','delimiter','\n','whitespace','');
idx1=find(cellfun(@isempty,strfind(txt,'Inferred ancestry of individuals:'))==0);
idx2=find(cellfun(@isempty,strfind(txt,'Estimated Allele Frequencies in each cluster'))==0);
idx3=find(cellfun(@isempty,strfind(txt,'Estimated Ln Prob of Data   ='))==0);

LnPtxt=txt{idx3};
LnPtxt=LnPtxt(31:end);
sumout.LnP=str2double(LnPtxt);


txt=txt(idx1+2:idx2-3);
C=[];
for k=1:length(txt);
    linetxt=txt{k};
    linefrc=strread(linetxt(strfind(linetxt,':')+2:end),'%f')';
    [itemp,iname]=strread(linetxt(1:strfind(linetxt,'(')-1),'%d%s');
    sumout.idvnames{k}=iname{1};
    C=[C;linefrc];
end