function PGEGUI
%PGEGUI

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

ButtonName = questdlg('What is your data type?', ...
    'Population Genetics & Evolution Toolbox', ...
    'Genotype', 'Haplotype', 'Sequence', 'Genotype');
switch ButtonName
    case 'Genotype'
        snptoolg
    case 'Haplotype'
        snptoolh
    case 'Sequence'
        seqtool
end
