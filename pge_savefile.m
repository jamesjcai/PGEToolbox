function pge_savefile(aln, filename, formatid)
%PGE_SAVEFILE - Save sequence file

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin < 1), error('Need alignment.'); end

if (nargin == 1),
    [filename, pathname, filterindex] = uiputfile( ...
        {'*.fasta;*.fas', 'FASTA Files (*.fasta, *.fas)'; ...
        '*.phylip;*.phy', 'PHYLIP Files (*.phylip, *.phy)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Save as');
    if (filterindex == 3),
        formatid = i_ask4formatid;
    else
        formatid = filterindex;
    end
    if (formatid == 3), aln = [];
        return;
    end
    if isequal(filename, 0) || isequal(pathname, 0)
        aln = [];
        return;
    else
        filename = fullfile(pathname, filename);
    end
end

if (nargin == 2), formatid = i_ask4formatid; end

switch (formatid)
    case (1)
        i_writefasta(aln, filename);
    case (2)
        i_writephylip(aln, filename);
    case (3)
        disp('unknonw file format');
end


    function i_writefasta(aln, filename)

        [n, m] = size(aln.seq);
        p = 1:n;
        q = 1:m;
        [NT, AA] = seqcode;

        switch (aln.seqtype)
            case (3) % Protein
                aln.seq(find(isnan(aln.seq))) = i_getcode4gap('PROTEIN');
                % AA = 'ARNDCQEGHILKMFPSTWYV*-';
                % Seq(p,q)=AA(aln.seq(p,q));
                Seq = AA(aln.seq);
            otherwise % nucleotides
                aln.seq(find(isnan(aln.seq))) = i_getcode4gap('DNA');
                % NT = ['ACGT-'];
                Seq(p, q) = NT(aln.seq(p, q));
        end

        fid = fopen(filename, 'wt');
        if fid == -1
            disp('Unable to open file.');
            return;
        end
        mt = 1:60:size(Seq, 2);
        mt = cat(1, mt', size(Seq, 2)+1);
        for i = 1:n,
            name = char(aln.seqnames(i));
            fprintf(fid, ['>%s\n'], name);
            for (j = 1:length(mt) - 1),
                % fprintf(fid, ['%s\n'],char(Seq(i,:)));
                fprintf(fid, ['%s\n'], char(Seq(i, [mt(j):mt(j+1) - 1])));
            end
        end
        fclose(fid);


            function i_writephylip(aln, filename)

                [n, m] = size(aln.seq);
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

                fprintf(fid, [' %d %d\n'], n, m);
                mt = 1:60:size(Seq, 2);
                mt = cat(1, mt', size(Seq, 2)+1);

                for (j = 1:length(mt) - 1),
                    for (i = 1:n),
                        s = char(Seq(i, [mt(j):mt(j+1) - 1]));

                        idx = 1:3:length(s);
                        idx = cat(1, idx', length(s)+1);
                        ss = '';
                        for (k = 1:length(idx) - 1),
                            ss = [ss, ' ', char(s([idx(k):idx(k+1) - 1]))];
                        end

                        if (j == 1)
                            name = char(aln.seqnames(i));
                            fprintf(fid, ['%10s %s\n'], i_name10(name), ss);
                        else
                            fprintf(fid, ['%10s %s\n'], ' ', ss);
                        end
                    end
                    fprintf(fid, '\n');
                end

                fclose(fid);