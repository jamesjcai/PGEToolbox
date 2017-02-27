function a=snp_ancallele(rsid)
%SNP_ANCALLELE - ancestral allele of a SNP (info from dbSNP)
% Syntax: a=snp_ancallele(rsid)
%
% The original function is written by Dr. Ence Yang

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-04-28 14:12:08 -0500 (Sun, 28 Apr 2013) $
% $LastChangedRevision: 533 $
% $LastChangedBy: jcai $

if ~strcmpi(rsid(1:2),'rs')
    error(['ERROR: snp_id does not start with "rs":',rsid]);
end
url=['http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=',rsid(3:end)];
cnt=urlread(url);
ix=strfind(cnt,'Ancestral Allele:');
if cnt(strfind(cnt,'Ancestral Allele:')+70)=='<'
    a=cnt(strfind(cnt,'Ancestral Allele:')+69);
elseif isempty(ix)
    if ~strncmp(cnt,'This snp_id was merged into',27)
        error(['ERROR: unknown content of url:',url]);
    end
    ix1=find(cnt=='<');
    ix2=find(cnt=='>');
    new_snp_id=cnt(ix2(1)+1:ix1(2)-1);
    fprintf('Waring: snp id "%s" has been merged into snp id "%s"\n',rsid,new_snp_id);
    a=snp_ancallele(new_snp_id);
else
    a='';   
end

