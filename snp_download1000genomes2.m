function [genodata,markinfo] = snp_download1000genomes2(query,popid)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<1
    [query]=select1000GenomesRegion;
end
if nargin<2
    popid='ALL';
end
chrstr=i_getchrstr(query);

genodata=[]; markinfo=[];

url=['http://browser.1000genomes.org/Homo_sapiens/UserData/Wizard?region=',...
     query,'&url=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr',...
     chrstr,'.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz&vcffilter=1&defaults=0&_backtrack=SelectSlice&wizard_next=/Homo_sapiens/UserData/SliceFile',...
     '&ind_select=NA18525,NA18526,NA18527,NA18528,NA18530,NA18532,NA18534,NA18535,NA18536,NA18537,NA18538,NA18539,NA18541,NA18542,NA18543,NA18544,NA18545,NA18546,NA18547,NA18548,NA18549,NA18550,NA18552,NA18553,NA18555,NA18557,NA18558,NA18559,NA18560,NA18561,NA18562,NA18563,NA18564,NA18565,NA18566,NA18567,NA18570,NA18571,NA18572,NA18573,NA18574,NA18576,NA18577,NA18579,NA18582,NA18592,NA18593,NA18595,NA18596,NA18597,NA18599,NA18602,NA18603,NA18605,NA18606,NA18608,NA18609,NA18610,NA18611,NA18612,NA18613,NA18614,NA18615,NA18616,NA18617,NA18618,NA18619,NA18620,NA18621,NA18622,NA18623,NA18624,NA18626,NA18627,NA18628,NA18630,NA18631,NA18632,NA18633,NA18634,NA18635,NA18636,NA18637,NA18638,NA18639,NA18640,NA18641,NA18642,NA18643,NA18645,NA18647,NA18740,NA18745,NA18747,NA18748,NA18749,NA18757',...
     '&wizard_submit=Next >'];

%s=urlread(url);
[stat,h]=web(url);
if stat~=0, return; end

%p_url=urlread('http://www.bioinformatics.org/pgetoolbox/snp_downloadhapmap_url');
%{
p_url=sprintf('ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz',chrstr);
p_region=query;
p_panelurl='ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel';

params = {'filter_screen','1','defaults','0','filter','0','Submit','Next >',...
          'url',p_url,'region',p_region,'panelurl',p_panelurl};
s=urlread('http://browser.1000genomes.org/Homo_sapiens/UserData/SelectSlice','GET',params);
%}

%%

%http://browser.1000genomes.org/tmp/slicer/1.1-50000.ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz

pw1=pwd;
cdpge;
cd('addins');
cd('gzip');

ME=1; t=1;
downURL=sprintf('http://browser.1000genomes.org/tmp/slicer/filtered_%s.ALL.chr%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz',...
                strrep(query,':','.'),chrstr);
           
while ~isempty(ME)&&t<5
    try
        fprintf('Trying to download %s\n',downURL);
        urlwrite(downURL,'infile.vcf.gz');
        ME=[];
    catch ME
        if strcmp(ME.identifier,'MATLAB:urlwrite:FileNotFound')
            disp(ME.message)
            fprintf('\nWait for %d seconds and try again...(%d/5)\n',5^t,t+1);
            pause(5^t);
            t=t+1;
        end        
    end
end
if isempty(ME)
    fprintf('Downloaded %s\n',downURL);
else
    fprintf('Unable to download %s\n',downURL);
    return;
end

%pw1=fileparts(which(mfilename));
[status]=system('gzip -f -d infile.vcf.gz');
%gunzip('infile.vcf.gz')   % not use this because of a java bug
if status==0
    fprintf('Unzipped infile.vcf.gz\n');
    [genodata,markinfo]=snp_readvcf('infile.vcf');
    fprintf('Read infile.vcf\n');
end
cd(pw1);
close(h);

%%
%gunzip(downURL)
%snp_readvcf('1.1-50000.ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf')

function chrstr=i_getchrstr(region)
    c=regexp(region,':','split');
    chrstr=c{1};
    