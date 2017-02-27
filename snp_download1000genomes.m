function [genodata,markinfo] = snp_download1000genomes(query,popid,warnwebwindow)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-12-27 00:25:00 -0600 (Fri, 27 Dec 2013) $
% $LastChangedRevision: 755 $
% $LastChangedBy: jcai $

if nargin<1
    [query]=select1000GenomesRegion;
end
if nargin<2, popid='ALL'; end
if nargin<3, warnwebwindow=false; end

chrstr=i_getchrstr(query);
genodata=[]; markinfo=[];


%http://browser.1000genomes.org/tmp/slicer/1.1-50000.ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
olddir=pwd;
cdpge;
%pw1=fileparts(which(mfilename));
%cd(pw1);
cd('addins/gzip');


ME=1; t=1;
downURL=sprintf('http://browser.1000genomes.org/tmp/slicer/%s.ALL.chr%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz',...
                 strrep(query,':','.'),chrstr);
downFile=sprintf('%s.ALL.chr%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz',...
                 strrep(query,':','.'),chrstr);
unzipFile=sprintf('%s.ALL.chr%s.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf',...
                 strrep(query,':','.'),chrstr);
             

if ~exist(downFile,'file')&&~exist(unzipFile,'file')
    
    %{
    'http://browser.1000genomes.org/Homo_sapiens/UserData/SelectSlice?
    filter_screen=1&
    defaults=0&
    url=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz&
    region=1:1-50000&
    filter=0&
    panelurl=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel&
    submit=Next >'
    %}

    %['http://browser.1000genomes.org/Homo_sapiens/UserData/SliceFeedback?baisize=;bai=;newsize=ND;newname=',...
    %  strrep(query,':','.'),'.ALL.chr',...
    %  chrstr,'.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz;region=',
    %  urlencode(query)];

    %'http://browser.1000genomes.org/Homo_sapiens/UserData/SelectSlice?filter_screen=1&defaults=0&url=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz&region=1:1-50000&filter=0&panelurl=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel&submit=Next >

    url=['http://browser.1000genomes.org/Homo_sapiens/UserData/SelectSlice?filter_screen=1&defaults=0&url=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr',...
        chrstr,'.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz&region=',...
        query,'&filter=0&panelurl=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel&submit=Next >'];

    if warnwebwindow
        uiwait(helpdlg('A web browser will open -- wait until the web page is completely loaded. Close the web brower. Genotype data will be automatically downloaded.'));
    end
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
    while ~isempty(ME)&&t<3
        try
            fprintf('\n====== SNP_DOWNLOAD1000GENOMES ======\n');
            fprintf('Downloading zipped genotype data file %s\n',downURL);
            urlwrite(downURL,downFile);
            ME=[];
        catch ME
            if strcmp(ME.identifier,'MATLAB:urlwrite:FileNotFound')
                disp(ME.message)
                fprintf('\nWait for %d seconds and try again...(%d/3)\n',5^t,t+1);
                pause(5^t);
                t=t+1;
            end
        end
    end
    if isempty(ME)
        fprintf('Genotype file downloaded and saved as %s\n located at %s\n',downFile,pwd);
    else
        fprintf('Unable to download gzipped genotype file at %s\n',downURL);
        return;
    end
end


if ~exist(unzipFile,'file')
    %pw1=fileparts(which(mfilename));
    fprintf('Unzip file %s\n',downFile);
    if ispc
        [status]=system(sprintf('gzip -f -d -k %s',downFile));   % -k keep the input file
    else
        [status]=system(sprintf('gzip -f -d %s',downFile));
    end
    %gunzip('infile.vcf.gz')   % not use this because of a java bug
    if status==0
        fprintf('File unzipped.\n');
    end
end

try
    [genodata,markinfo]=snp_readvcf(unzipFile);
catch ME
    fprintf(ME.message);
    return;
end
fprintf('Reading VCF file %s\n',unzipFile);

YesNo = any(strcmp(who,'h'));

if YesNo&&exist(h,'var')&&~isempty(h)
    try
        close(h);
    catch
    end
end

try
    cd(olddir);
catch
end

%%
%gunzip(downURL)
%snp_readvcf('1.1-50000.ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf')

function chrstr=i_getchrstr(region)
    c=regexp(region,':','split');
    chrstr=c{1};
    