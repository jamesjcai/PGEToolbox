function pge_printaln(aln, filename)
%PGE_PRINTALN - Write alignment structure into a numbered sequential PHYLIP formatted file
%
% Syntax:  pge_printaln(aln,'filename')
%
% Inputs:
%    aln         - Alignment structure
%    filename    - Sequencial PHYLIP file.
%
% See also:

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if nargin < 2
    [filename, pathname, filterindex] = uiputfile( ...
        {'*.phylip;*.phy', 'Phylip Format Files (*.phylip, *.phy)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Save as');
    if ~(filename), return; end
    filename = [pathname, filename];

    if (filterindex == 1)
        if (isempty(find(filename == '.')))
            filename = [filename, '.phy'];
        end
    end
end


[n, m] = size(aln.seq);
if ~isfield(aln, 'pos')
    pos = 1:m;
else
    pos = aln.pos;
end

p = 1:n;
q = 1:m;
[NT, AA] = seqcode;

switch (aln.seqtype)
    case (3) % Protein
        aln.seq(find(isnan(aln.seq))) = i_getcode4gap('PROTEIN');
        Seq(p, q) = AA(aln.seq(p, q));
    otherwise % nucleotides
        aln.seq(find(isnan(aln.seq))) = i_getcode4gap('DNA');
        Seq(p, q) = NT(aln.seq(p, q));
end


fid = fopen(filename, 'wt');
if fid == -1
    disp('Unable to open file.');
    return
end


vlabel = numvlabel(pos);

fprintf(fid, [' %d %d\n'], n, m);
for (k = 1:length(vlabel)),
    fprintf(fid, '%10s  %s\n', ' ', vlabel{k});
end

for i = 1:n,
    name = i_name10(char(aln.seqnames(i)));
    xseq = Seq(i, :);
    if (i == 1)
        fprintf(fid, '%10s  %s\n', name, char(xseq));
    else
        idx = xseq == Seq(1, :);
        xseq(idx) = '.';
        fprintf(fid, '%10s  %s\n', name, char(xseq));
    end
end
fclose(fid);


%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

    function [n] = i_getcode4gap(seqtype)

        % Molecular Biology & Evolution Toolbox, (C) 2005
        % Author: James Cai
        % Email: jamescai@hkusua.hku.hk
        % Website: http://web.hku.hk/~jamescai/
        % Last revision: 5/28/2005


        switch (upper(seqtype))
            case {'PROTEIN'}
                [NT, AA] = seqcode;
                n = find(AA == '-');
                % n=22;
            case {'DNA'}
                [NT] = seqcode;
                n = find(NT == '-');
                % n=5;
            otherwise
                error('Wrong seqtype.')
        end


            function [name2] = i_name10(name)

                % Molecular Biology & Evolution Toolbox, (C) 2005
                % Author: James Cai
                % Email: jamescai@hkusua.hku.hk
                % Website: http://web.hku.hk/~jamescai/
                % Last revision: 5/28/2005

                len = length(name);
                if (len < 10)
                    name(len+1:10) = ' ';
                elseif (len > 10)
                    name = char(name(1:10));
                end
                name2 = name;