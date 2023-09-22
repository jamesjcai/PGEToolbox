function [genodata, markinfo] = snp_readperlegen(filename, noise)
%SNP_READPERLEGEN - Reads Perlegen genotype data
%[genodata,markinfo] = snp_readperlegen(filename,noise)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 1
    [filename, pathname] = uigetfile( ...
        {'*.perlegen;*.per;*.dat;*.txt', 'Perlegen Files (*.perlegen, *.per, *.dat, *.txt)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Pick a Perlegen (SNP genotype data) format file');
    if ~(filename), genodata = [];
        markinfo = [];
        return;
    end
    filename = [pathname, filename];
end

if nargin < 2, noise = 0; end

if ~(exist(filename, 'file') || exist(fullfile(cd, filename), 'file')),
    error('InvalidInput', 'Input must be a valid file')
end


%    FID = fopen(filename);
%    hapmapTline = fgets(FID);         %read header line
%	word = gettext(Tline);
%	if strcmp(word,'rs#')==1
%	  Tline = fgets(FID);
%	end
%    fclose(FID)

if noise,
    disp(['Reading "PERLEGEN" file ', filename, ' ...'])
end

[rsid, alleles, chrom, pos, strand, assemblyid, center, marker] = textread(filename, ...
    '%s%s%s%d%s%s%s%[ATCGN\t ]%*[^\n]', 'delimiter', '\t', ...
    'commentstyle', 'shell', 'headerlines', 5);


%markerinfo= struct('rsid',rsid,'alleles',alleles,'chrom',chrom,'pos',pos,'strand',strand,'assemblyid',assemblyid,'center',center,...
%'protLSID',protLSID,'assayLSID',assayLSID,'panelLSID',panelLSID,'QCcode',QCcode);

n = length(marker);
genodata = [];
for (k = 1:n),
    %      marker{k}
    %      pause
    %      sscanf(marker{k},'%s')
    %      pause
    genodata = [genodata; sscanf(marker{k}, '%s')];
end
genodata = i_perlegen_transpose(genodata);
genodata = i_numeric(genodata);


markinfo.allele = alleles;
markinfo.strand = strand;
markinfo.rsid = rsid;
markinfo.chr = chrom;
markinfo.pos = pos;

%indvinfo.fatherid=zeros(1,size(genodata,1));
%indvinfo.motherid=zeros(1,size(genodata,1));


    function [geno_out] = i_perlegen_transpose(geno_in)
        %geno_in=['GCGGGC';'AAATAT'];
        [n, m] = size(geno_in);
        %geno_out=char(ones(m/2,n*2)*42);
        geno_out = zeros(m/2, n*2);

        for (j = 1:n),
            for (k = 1:2:m),
                x = (k + 1) / 2;
                y = j * 2 - 1;
                geno_out(x, [y, y + 1]) = geno_in(j, [k, k + 1]);
            end
        end


            function [genodata] = i_numeric(genodata)
                genodata = double(genodata);
                genodata(find(genodata == double('A'))) = 1;
                genodata(find(genodata == double('C'))) = 2;
                genodata(find(genodata == double('G'))) = 3;
                genodata(find(genodata == double('T'))) = 4;
                genodata(find(genodata == double('N'))) = 5;
                genodata(genodata < 1 | genodata > 4) = 5;
