function writefasta(aln,filename)
%WRITEFASTA - Write alignment structure into a FASTA formatted file
%
% Syntax:  writefasta(aln,'filename')
%
% Inputs:
%    aln         - Alignment structure
%    filename    - FASTA formatted file (ASCII text file).
%
% See also: READFASTA

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if nargin < 2
    [filename, pathname,filterindex] = uiputfile( ...
       {'*.fasta;*.fas', 'FASTA Format Files (*.fasta, *.fas)';
        '*.*',  'All Files (*.*)'}, ...
        'Save as');
	if ~(filename), return; end
	filename=[pathname,filename];
	if filterindex==1
		if (isempty(find(filename=='.', 1)))
		filename=[filename,'.fas'];
		end
	end

end

s=aln.seq;
n=size(s,1);
[NT,AA] = seqcode;
switch (aln.seqtype)
    case (3)	% Protein
        s(isnan(s))=i_getcode4gap('PROTEIN');
        SEQ=AA;
    otherwise	% nucleotides
        s(isnan(s))=i_getcode4gap('DNA');
        SEQ=NT;
end

fid = fopen(filename,'wt');
if fid == -1
   error('Unable to open file.');
end

mt = 1:50:size(s,2);
mt = cat(1,mt',size(s,2)+1);

for i=1:n,
    fprintf(fid, '>%s\n',aln.seqnames{i});
	for j=1:length(mt)-1
	   fprintf(fid, '%s\n', SEQ(s(i,mt(j):mt(j+1)-1)));
	end
end
fclose(fid);
