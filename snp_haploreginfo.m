function snp_haploreginfo(rsid)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin<1),
   prompt={sprintf('HaploReg is a tool for exploring annotations of a SNP\n Please enter the SNP RSID: ')};
   def={'rs2693665'};
   dlgTitle='Input for SNP RSID';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);
 
    if ~(isempty(answer)),
      rsid=answer{1};
    else
         return;
    end
end

rsidx=sprintf('%s',rsid);
urlFetch=sprintf('http://www.broadinstitute.org/mammals/haploreg/detail.php?query=&id=%s',...
    rsidx);
web(urlFetch,'-browser')

%try
%    pagecontent=urlread(urlFetch);
%catch ME
%    error(ME.message);
%end
    
