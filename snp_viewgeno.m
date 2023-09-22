function snp_viewgeno(geno, rowflag)
%SNP_VIEWGENO - Displays SNP genotype information
%snp_viewgeno(geno,'markers')

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin == 1), rowflag = []; end

if ~(isempty(rowflag))
    if ~strcmpi(rowflag, 'markers')
        error('MATLAB:UNIQUE:UnknownRowFlag', 'Unknown rowflag.');
    end
    rowtag = 'Mrk_';
else
    rowtag = 'Idv_';
end


[n, m] = size(geno);
%if (nargin<2),
names = {};
for k = 1:n
    %names{k}=['Idv_',sprintf('%d',k)];
    names{k} = [rowtag, num2str(k)];
end
%end

geno(find(geno == 1)) = 'A';
geno(find(geno == 2)) = 'C';
geno(find(geno == 3)) = 'G';
geno(find(geno == 4)) = 'T';
geno(find(geno == 5)) = 'N';
geno = char(geno);

i_dispheader('Genotype View')
for (i = 1:n),
    ssx = sprintf('%10s\t', names{i});
    for k = 1:m / 2
        ssx = [ssx, sprintf('%2s ', geno(i, [2 * k - 1, 2 * k]))];
    end
    disp(ssx)
end
i_dispfooter
