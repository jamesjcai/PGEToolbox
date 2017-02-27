function [marker]=inputMarker()

items(1).name = 'SNP id:   ';
items(1).default = 0;        
items(1).values = {''};
items(1).help = 'Input a SNP identification such as, rs006859.';
prompt='Input a SNP identification such as, rs006859.';

dlg_title = 'Download from WWW.HAPMAP.ORG';
msg = sprintf('The wildcard character * is allowed. \n\nExamples: Chr20, Chr9:660,000..760,000, SNP:rs6870660,\n NM_153254, BRCA2, 5q31, ENm010.');
out = inputdlg(msg,dlg_title,1)

if ~(isempty(out)),
	marker=out{1};
else
	marker='';
end
