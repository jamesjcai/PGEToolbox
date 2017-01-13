function [D,txt]=snp_dbsnpinfo(rsid,datasrc)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<2
    datasrc='FLT';
end
%datasrc='DocSet';

rsidx=sprintf('%d,',rsid);
switch datasrc
    case 'FLT'
urlFetch=sprintf('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&retmode=text&report=FLT&id=%s',...
    rsidx);
    case 'DocSet'
urlFetch=sprintf('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&retmode=text&report=DocSet&id=%s,',...
    rsidx);
end

try
    pagecontent=urlread(urlFetch);
catch
    %errordlg(lasterr)
    disp(urlFetch)
    error(lasterror)
end
    
txt = char(strread(pagecontent,'%s','delimiter','\n','whitespace',''));
txt = cellstr(txt);
mt = find(cellfun('isempty',txt));
txt(mt)=[];
idx=find(cellfun('isempty',strfind(txt,'SNP | alleles'))==0);


    
if size(idx,1)~=length(rsid)
    urlFetch
    rsid
    txt(idx)
    error('SNP_DBSNPINFO cannot parse retrieved information.')
end
txt=txt(idx);



for k=1:length(txt)
    allele{k}=txt{k}(:,16:18);
end    
D.allele=allele;