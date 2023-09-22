function [ldinfo] = emldrun(genodata)
%EMLDRUN - run EMLD to calculate pairwise LD between SNPs
%
% [ldinfo]=emldrun(genodata)
%
% SNPs are from independent individuals
%
% EMLD is a program to calculate pair-wise linkage disequlibrium from
% genotype data. EM algorithm is used to estimate haplotype frequencies.
%
%SEE ALSO: SNP_LDPLOT

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-04-23 19:44:14 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 529 $
% $LastChangedBy: jcai $

oldpath = pwd;
ldinfo = [];
cdpge;
cd 'addins/EMLD';
%[exedir,dlgshown]=pge_getprgmdir(sprintf('%s_prgmdir',mfilename));
%if isempty(exedir)||dlgshown, return; end
%cd(exedir);

filename = 'emldinput.txt';
try
    disp('Converting GENODATA into {1,2,0} format...')
    geno = snp_12geno(genodata);
    disp('Writing EMLD input file...')
    i_writelinkage(geno, filename);
    disp('Running java EMLD...')
    system(sprintf('java EMLD %s', filename));
catch exception
    ldinfo = [];
    return;
end

disp('Reading EMLD output...')
[~, ~, dv, dpv, r2v] = textread('LD.xt', '%d%d%f%f%f', 'headerlines', 1);
ldinfo = struct;
ldinfo.d = triu(squareform(dv'));
ldinfo.dprime = triu(squareform(dpv'));
ldinfo.r2 = triu(squareform(r2v'));
cd(oldpath);


    function i_writelinkage(geno, filename)

        fid = fopen(filename, 'wt');
        if fid == -1
            warning('Unable to open file.');
            return;
        end
        [samplen, marklen] = snp_samplen(geno);
        indvlen = samplen / 2;
        fprintf(fid, '%d\t%d\n', indvlen, marklen);

        for k = 1:indvlen
            fprintf(fid, '100%d\t', k);
            for j = 1:marklen * 2
                fprintf(fid, '%d\t', geno(k, j));
            end
            fprintf(fid, '\n');
        end
        fclose(fid);