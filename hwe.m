function [P] = hwe(nhom1, nhet, nhom2)
%HWE - Exact test of Hardy-Weinberg Equilibrium (Wigginton et al. 2005)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if (nargin < 1)
    prompt = {'Numbers of AA:', 'Numbers of Aa:', 'Numbers of aa:'};
    def = {'1467', '381', '23'};
    dlgTitle = 'Enter the observed numbers of phenotypes';
    lineNo = 1;
    answer = inputdlg(prompt, dlgTitle, lineNo, def);
    if (isempty(answer)), return; else
        nhom1 = str2num(char(answer{1}));
        nhet = str2num(char(answer{2}));
        nhom2 = str2num(char(answer{3}));
    end
elseif (nargin == 1)
    inputd = nhom1;
    [n, m] = size(inputd);
    if (n == 1 && m == 3) || (n == 3 && m == 1)
        nhom1 = inputd(1);
        nhet = inputd(2);
        nhom2 = inputd(3);
    else
        return;
    end
end


P = snphwe_mex(nhom1, nhet, nhom2);
if nargout == 1
    return;
end

% The Hardy-Weinberg law states that the frequencies of the genotypes AA, Aa,
% and aa will be equal to p2, 2pq, and q2 after a single generation, regardless
% of the initial genotype frequencies, and that they will remain in these
% frequencies forever.

nA = 2 * nhom1 + nhet;
na = 2 * nhom2 + nhet;

p = nA ./ (nA + na);
q = 1 - p;

% If we know p and q, then we can calculate the freq. of AA (p^2)
% Aa (2pq) and aa (q^2) that would be expected if the population
% is in HWE as follows:
t = nhom1 + nhet + nhom2;
eAA = (p * p) * t;
eAa = (2 * p * q) * t;
eaa = (q * q) * t;

if (eAA == 0 || eAa == 0 || eaa == 0)
    P1 = 1;
    P2 = 1;
    X2 = nan;
    X2C = nan;
else
    X2 = (nhom1 - eAA)^2 / eAA + (nhet - eAa)^2 / eAa + (nhom2 - eaa)^2 / eaa;
    X2C = (abs(nhom1-eAA) - 0.5)^2 / eAA + (abs(nhet-eAa) - 0.5)^2 / eAa ...
        +(abs(nhom2-eaa) - 0.5)^2 / eaa;
    df = 1;
    P1 = 1 - chi2cdf(X2, df);
    P2 = 1 - chi2cdf(X2C, df);
end

if (nargout < 1)
    disp(' ')
    disp('Chi-square tests of Hardy-Weinberg Equilibrium:')
    fprintf('----------------------------------------------\n');
    fprintf(['Chi-square value : %10.3f\n'], X2);
    fprintf(['      P-value : %10.5f\n'], P1);
    fprintf(['Chi-square with Yates'' correction : %10.3f\n'], X2C);
    fprintf(['      P-value : %10.5f\n'], P2);
    fprintf('----------------------------------------------\n');
    disp('Wigginton Exact Tests of Hardy-Weinberg Equilibrium');
    fprintf('      P-value : %10.5f\n', P);
    fprintf('----------------------------------------------\n');

end

%fprintf('SNP1,%d,%d,%d,%d,%d,%d\n',nhom1,nhet,nhom2,round(eAA),round(eAa),round(eaa))


%{

/*
    // This function implements an exact SNP test of Hardy-Weinberg
        // Equilibrium as described in Wigginton, JE, Cutler, DJ, and
        // Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
        // Equilibrium. American Journal of Human Genetics. 76: 000 - 000
        //
        // Written by Jan Wigginton
        */

        %}
