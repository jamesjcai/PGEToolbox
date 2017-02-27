function snptool
%SNPTOOL - snptool

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

   ButtonName = questdlg('What is your data type?', ...
                         'SNP Tool', ...
                         'Genotype', 'Haplotype', 'Genotype');
   switch ButtonName,
     case 'Genotype',
      snptoolg
     case 'Haplotype',
     snptoolh
   end
