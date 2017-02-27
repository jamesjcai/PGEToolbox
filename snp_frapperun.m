function [res]=snp_frapperun(genodata,k)
%SNP_FRAPPE - runs frappe
%
%snp_frapperun(genodata)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<2
    k=3;
end

oldpath=pwd;
%cdpge; cd('addins/frappe');
[exedir,dlgshown]=pge_getprgmdir(sprintf('%s_prgmdir',mfilename));
if isempty(exedir)||dlgshown, return; end
cd(exedir);

fprintf('Writing input file...');
status1=i_writefrappeinput(genodata,'input.ped');
status2=i_writefrappeparam(genodata,'param.txt',k);

if status1==1 && status2==1    % good
    fprintf('done.\n');
else
    cd(oldpath);
    return;
end

if ispc
    cplotcmd='frappe_1.1.exe param.txt';
elseif ismac
    cplotcmd='./frappe.MAC param.txt';
else
    cplotcmd='./frappe.LINUX param.txt';
end

fprintf('Running FRAPPE...');
try
[~,~]=system(cplotcmd);
catch ME
    warning(ME.message);
    cd(oldpath);
    return;
end
fprintf('done.\n');
%res=i_i_writefrappeoutput('input_result.txt');
x=importdata('input_result.txt');
res=x.data;
cd(oldpath);
bar(res,1,'stacked','edgecolor','none')
axis ij;
ylim([0 1])
xlim([0 size(genodata,1)+1])


function [status]=i_writefrappeparam(g,f,k)
    fid=fopen(f,'w');
    fprintf(fid,'GenotypeFile  ="input.ped" ## Mandatory genotype file name\n');
    fprintf(fid,'MaxIter=  10000    ## maximum number of EM iterations\n');
    fprintf(fid,'K      = %d  ## Number of ancestral individuals\n',k);
    fprintf(fid,'M      = %d  ## Number of markers\n',size(g,2)/2);
    fprintf(fid,'I      = %d  ## Number of individuals;\n',size(g,1));
    fprintf(fid,'step   =  200 ### Define either step or Nout\n');
    fprintf(fid,'Nout   =  0\n');
    %fprintf(fid,'IndividualFile = "simIndfile.txt"  ### Optional file \n');
    %fprintf(fid,'printP	= 1	## optional \n');
    %fprintf(fid,'threshold= 10000  ## optional convergence threshold\n');
    fclose(fid);
    status=1;


function [status]=i_writefrappeinput(geno,filename)
%SNP_WRITELINKAGE - saves as linkage format
%snp_writelinkage(geno,mark,filename)
if (isempty(geno)), status=0; return; end
fid = fopen(filename,'wt');
if fid==-1, status=0; end

[samplen,marklen]=snp_samplen(geno);
indvlen=samplen/2;

ACGT='12340';
for k=1:indvlen
      fprintf(fid,'%d %d 0 0 1 1 ',k,k);
      for j=1:marklen*2-1
	      fprintf(fid,'%s ',ACGT(geno(k,j)));
      end
      fprintf(fid,'%s\n',ACGT(geno(k,marklen*2)));
      %fprintf(fid,'\n');
end
fclose(fid);
status=1;




