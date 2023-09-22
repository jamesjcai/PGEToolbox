function [h, hraw] = snp_faywu00h(geno, p)
%SNP_FAYWU00H - Fay & Wu's H statistics for SNP
% Syntax: [h]=snp_faywu00h(geno,p)
%
% p   -  derived allele frequency, see SNP_DAF

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

% ref: http://genome.cshlp.org/content/15/11/1553.long

if nargin < 2
    p = []; % derived allele frequency, see SNP_DAF
end


nsam = size(geno, 1) * 2; % num of gametes

[theh] = snp_thetah(geno, [], 0, p);
[thew] = snp_thetaw(geno, [], 0);
[thepi] = snp_thetapi(geno, [], 0);
[Sn] = snp_segsites(geno);
%[h] = faywu00h(smpln, Sn, theh);

hraw = thepi - theh;

nx = 1:(nsam - 1);
an = sum(1./nx);
bn = sum(1./(nx.^2));
bn2 = sum(1./((1:nsam).^2));
t1 = Sn ./ an;
t2 = Sn * (Sn - 1) ./ (an.^2 + bn);

n = nsam;
hvar = t1 * (n - 2) / (6 * (n - 1)) + t2 * ((18 * n^2) * (3 * n + 2) * bn2 ...
    -(88 * n^3 + 9 * n^2 - 13 * n + 6)) / (9 * n * (n - 1)^2);


% if hvar>0
if Sn >= 1
    h = 0.5 * hraw ./ sqrt(hvar);
else
    % hvar==0
    %hraw
    %snp_viewgeno(geno)
    %pause
    h = nan;
end


if (nargout < 1),
    i_dispheader('Fay and Wu''s H Test')
    disp('Mode: SNP');
    fprintf('\n');
    fprintf('No. of Segregating sites (S): %d\n', Sn);
    fprintf('Fay''s theta_H: %f\n', theh);
    fprintf('\n');
    %fprintf ('Diff = %f, s.e. = %f\n', Diff, DiffSE);
    fprintf('Fay and Wu''s H (not normalized): %f\n', hraw);
    fprintf('Fay and Wu''s H (normalized): %f\n', h);
    [H] = faywu00h_simu(smpln, 1000, thew, 0);
    p = sum(hraw > H) ./ 1000;
    fprintf('Statistical significance:\n P = %f%s\n', p, sigtag(p));
    i_dispfooter
end
