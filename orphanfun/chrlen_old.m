function [bp,LMap] = chrlen_old(id)
%CHRLEN - returns human chromosome length
%USAGE: [bp] = chrlen(id)
%
%-Version: from UCSC Genome Browser based on Human Mar. 2006 Assembly.
% http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=91767390&chromInfoPage=
% Human Mar. 2006 (hg18) Browser Sequences

% Population Genetics & Evolution Toolbox, (C) 2009
% Author: James J. Cai
% Email: jamescai@stanford.edu
% Website: http://bioinformatics.org/pgetoolbox/
% Last revision: 10/30/2009

if (ischar(id)), 
	id=upper(id);
	if (strcmp(id,'23')), id='X'; end
	if (strcmp(id,'24')), id='Y'; end
elseif(id==23), 
	id='X';
elseif(id==24), 
	id='Y';
else
	id=int2str(id); 
end

clen=[247199719, 135374737, 134452384, 132289534, 114127980,...
      106360585, 100338915, 88822254, 78654742, 76117153, 63806651,...
      242751149, 62435964, 46944323, 49591432, 199446827, 191263063,...
      180837866, 170896992, 158821424, 146274826, 140273252, 154913754,...
      54733917];
cnam = {'1','10','11','12','13','14','15','16','17','18','19','2','20','21','22','3','4','5','6','7','8','9','X','Y'};



[x,y]=ismember(upper(id),cnam);
if (y>0)
    bp=clen(y);
else
    error(sprintf('%s is not a valid chromosome name.',id));
end


if (nargout>1),
	LMap={	'1', 247199719;
		'2', 242751149;
		'3', 199446827;
		'4', 191263063;
		'5', 180837866;
		'6', 170896992;
		'7', 158821424;
		'8', 146274826;
		'9', 140273252;
		'10', 135374737;
		'11', 134452384;
		'12', 132289534;
		'13', 114127980;
		'14', 106360585;
		'15', 100338915;
		'16', 88822254;
		'17', 78654742;
		'18', 76117153;
		'19', 63806651;
		'20', 62435964;
		'21', 46944323;
		'22', 49591432;
		'X', 154913754;
		'Y', 54733917};
end

