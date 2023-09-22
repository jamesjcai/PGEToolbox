function [het, fst] = fdist2run(ndem, npop, fst0, nsmp, muttype, nrep)
%FDISTRUN - executes FDIST2
%
%[het,fst] = fdist2run(ndem,npop,fst0,nsmp,muttype,nrep)
%[het,fst] = fdist2run(ndem,npop,fst0,nsmp);


% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-04-23 19:44:14 -0500 (Tue, 23 Apr 2013) $
% $LastChangedRevision: 529 $
% $LastChangedBy: jcai $

% ndem - Total number of demes (100 max).
% npop - No of populations sampled (must be less than or equal to total number of demes).
% fst0 - Expected Fst for infinite allele, infinite island model.
% nsmp Sample size (assumed the same in all populations)
% type - Indicator - 1 for stepwise mutations; 0 for infinite alleles
% nrep - No of realizations (loci)
%
%Example:
%
% number of demes is 100
% number of samples  is 15
% expected (infinite allele, infinite island) Fst is 0.170000
% median sample size is 50 - suggest give 50 if >50
% infinite alleles mutation model assumed
% 20000 realizations (loci) requested

oldpath = pwd;
cdpge;
cd('addins/fdist2');
%[exedir,dlgshown]=pge_getprgmdir(sprintf('%s_prgmdir',mfilename));
%if isempty(exedir)||dlgshown, return; end
%cd(exedir);

%if (nargin<6), nrep=20000; end
if (nargin < 6), nrep = 2000; end
if (nargin < 5), muttype = 0; end
if (nargin < 4), nsmp = 50; end
if (nargin < 3), fst0 = 0.17; end
if (nargin < 2), npop = 15; end
if (nargin < 1), ndem = 20; end

if ndem > 100 || npop > ndem
    error('');
end

% NDEME - Total number of demes (100 max).
% NPOP - No of populations sampled (must be less than or equal to total number of demes).
% FST0 - Expected Fst for infinite allele, infinite island model.
% NSAMP - Sample size (assumed the same in all populations)
% muttype - Indicator - 1 for stepwise mutations; 0 for infinite alleles
% NLOCI - No of realizations (loci)
%fdist2(20,15,0.17,50,0,30)
%y=input('Load existing output (Yes=1; No=0): ');


y = 1;
if y
    x = textread('fdist2out.dat');
    het = x(:, 1);
    fst = x(:, 2);
    [ndem, npop, fst0, nsmp, muttype, nrep] = i_readpara;
else
    i_writepara(ndem, npop, fst0, nsmp, muttype, nrep);
    olddir = pwd;
    cdpge;
    cd('addin/fdist2');
    [het, fst] = fdist2(ndem, npop, fst0, nsmp, muttype, nrep);
    cd(olddir);
end

fprintf('\nNumber of demes: %d\n', ndem);
fprintf('Number of samples: %d\n', npop);
fprintf('Expected (infinite allele, infinite island) Fst: %f\n', fst0);
fprintf('Median sample size: %d - suggest give 50 if >50\n', nsmp);
    if muttype == 0
        fprintf('Infinite alleles mutation model assumed\n');
    else
        fprintf('Step-stone mutation model assumed\n');
    end
    fprintf('Realizations (loci) requested: %d\n', nrep)

    %fprintf('average Fst is %f\n',sum(fst.*het)./sum(het));
    %disp('Running FDIST2 CPLOT...')
    %fst(fst<0)=0;
    cd(oldpath);


    function [ndem, npop, fst0, nsmp, muttype, nrep] = i_readpara
        [ndem, npop, fst0, nsmp, muttype, nrep] = textread('fdist_params2.dat', ...
            '%d%d%f%d%d%d', 'delimiter', '\n');

            function i_writepara(ndem, npop, fst0, nsmp, muttype, nrep)
                fid = fopen('fdist_params2.dat', 'w');
                fprintf(fid, '%d\n', ndem);
                fprintf(fid, '%d\n', npop);
                fprintf(fid, '%f\n', fst0);
                fprintf(fid, '%d\n', nsmp);
                fprintf(fid, '%d\n', muttype);
                fprintf(fid, '%d\n', nrep);
                fclose(fid);
