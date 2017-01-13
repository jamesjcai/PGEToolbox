function snp_weblink(markinfo)
%SNP_WEBLINK - Website links for given SNP markers

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $



if (nargin<1),
   prompt={sprintf('This command will return the chromosomal position of a SNP.\n Please enter a SNP ID: ')};
   def={'rs2693665'};
   dlgTitle='Input for SNP ID';
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,def);

    if ~(isempty(answer)),
      crsid=answer{1};
    else
         return;
    end


filename=tempname;
filename(end-2:end)='htm';
fid=fopen(filename,'w');

fprintf(fid,'<!doctype html public "-//w3c//dtd html 4.0 transitional//en">\n');
fprintf(fid,'<html>\n');
fprintf(fid,' <head>\n');
fprintf(fid,'  <title>PGEToolbox SNP weblink</title>\n');
fprintf(fid,' </head>\n');
fprintf(fid,' <body><pre>\n');


fprintf(fid,[sprintf('%10s',crsid),...
	'\t[<a href="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=',...
         crsid,'" target="_blank">dbSNP</a>|<a href="http://hapmap.org/cgi-perl/gbrowse/hapmap27_B36/?name=SNP%%3A',...
         crsid,'" target="_blank">HapMap</a>|<a href="http://www.ensembl.org/Homo_sapiens/snpview?snp=',...
         crsid,'" target="_blank">Ensembl</a>|<a href="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&report=Brief&id=',...
         crsid(3:end),'" target="_blank">Brief</a>]\n']);



%      fprintf('%10s\t%s\t%s\t%s\t%s\t%d\n',['Mrk_',num2str(k)],rsids{k},alleles{k},chrom{k},strands{k},position(k));
fprintf(fid,'  </pre></body>\n');
fprintf(fid,'</html>\n');
fclose(fid);
%web(filename,'-browser');
web(filename);


else
    alleles=markinfo.allele;
    strands=markinfo.strand;
    rsids=markinfo.rsid;
    position=markinfo.pos;
    chrom=markinfo.chr;

filename=tempname;
filename(end-2:end)='htm';
fid=fopen(filename,'w');

fprintf(fid,'<!doctype html public "-//w3c//dtd html 4.0 transitional//en">\n');
fprintf(fid,'<html>\n');
fprintf(fid,' <head>\n');
fprintf(fid,'  <title>PGEToolbox SNP weblink</title>\n');
fprintf(fid,' </head>\n');
fprintf(fid,' <body><pre>\n');

for (k=1:length(position)),
crsid=rsids{k};
fprintf(fid,[sprintf('%10s',crsid),'\t',alleles{k},'\t',chrom{k},'\t',...
        strands{k},'\t', num2str(position(k)),...
	'\t[<a href="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=',...
         crsid,'" target="_blank">dbSNP</a>|<a href="http://hapmap.org/cgi-perl/gbrowse/hapmap_B36/?name=SNP%%3A',...
         crsid,'" target="_blank">HapMap</a>|<a href="http://www.ensembl.org/Homo_sapiens/snpview?snp=',...
         crsid,'" target="_blank">Ensembl</a>|<a href="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&report=Brief&id=',...
         crsid(2:end),'" target="_blank">Brief</a>]\n']);

%      fprintf('%10s\t%s\t%s\t%s\t%s\t%d\n',['Mrk_',num2str(k)],rsids{k},alleles{k},chrom{k},strands{k},position(k));

end
fprintf(fid,'  </pre></body>\n');
fprintf(fid,'</html>\n');
fclose(fid);
%web(filename,'-browser');
web(filename);
end
