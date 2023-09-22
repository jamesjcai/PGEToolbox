function [aln] = pge_openfile(filename, formatid)
%PGE_OPENFILE - Open sequence file

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 1
    [filename, pathname, filterindex] = uigetfile( ...
        {'*.fasta;*.fas', 'FASTA Files (*.fasta, *.fas)'; ...
        '*.phylip;*.phy', 'PHYLIP Files (*.phylip, *.phy)'; ...
        '*.mat', 'MAT-file (*.mat)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Open a file');

    if (filterindex == 4),
        formatid = i_ask4formatid;
    else
        formatid = filterindex;
    end

    if (formatid == 4), aln = [];
        return;
    end

    if isequal(filename, 0) || isequal(pathname, 0)
        aln = [];
        return;
    else
        filename = fullfile(pathname, filename);
    end

end


disp(['Reading ', filename]);
switch (formatid)
    case (1)
        [aln] = i_readfasta(filename);
    case (2)
        [aln] = i_readphylip(filename);
    case (3)
        [alnstrc] = load(filename, 'aln');
        aln = alnstrc.aln;
        return;
    case (4)
        disp('unknonw');
end


choice = questdlg('Select sequence type', ...
    'SEQUENCE TYPE:', ...
    'Non-coding', 'Coding', 'Coding');
% Handle response
switch choice
    case 'Non-coding'
        seqtype = 1;
    case 'Coding'
        seqtype = 2;
end

aln.seqtype = seqtype;
aln.geneticcode = 1;
m = size(aln.seq, 2);
aln.pos = 1:m;


    function [aln] = i_readphylip(filename)
        file = fopen(filename, 'r');
        [nm, x] = fscanf(file, '%d', [1, 2]);
        if (x ~= 2), error('NOT PHYLIP FORMAT'); end
        n = nm(1);
        m = nm(2);
        fclose(file);

        txt = textread(filename, '%s', 'delimiter', '\n', 'whitespace', '');
        % remove header line
        txt(1) = [];

        % remove empty lines
        while isempty(txt{1}), txt(1) = []; end

        % find first empty string in cell array, which occurs after the first
        % consensus line
        mt = find(cellfun('isempty', txt));
        % eliminate empty lines
        txt(mt) = [];
        num_seq = n;
        names = {};
        for s = 1:num_seq,
            string = txt{s};
            token = [];
            remainder = [];
            token = string(1:10);
            remainder = removeblanks(string(11:length(string)));
            names{s} = removeblanks(token);
            disp(token)

            for r = s + num_seq:num_seq:size(txt, 1),
                % make sure that there aren't sequence numbers at the end
                remainder = [remainder, removeblanks(txt{r})];
            end

            if ~(length(remainder) == m), error('INCORRECT PHYLIP FORMAT');
                return;
            end
            S(s, :) = upper(remainder);
        end


        NT = 'ACGT-';
        [n, m] = size(S);
        seq = ones(n, m) .* 5;
        for k = 1:5, seq(S == NT(k)) = k; end

        aln.seqnames = names;
        aln.seq = seq;

        aln.locus = ones(n, 1);
        aln.population = ones(n, 1);
        aln.count = ones(n, 1);


            function [aln] = i_readfasta(filename)
                MAXNAME = 200;
                file = fopen(filename, 'r');

                % Now we are looking for the maximum length of the sequence
                n = 0; % the number of sequences
                m = 0; % the maximum length
                cm = 0; % current sequence length

                while 1
                    [x, nr] = fscanf(file, '%c', 1);
                    if nr == 0 break;
                    end;
                    if x == '>' % new sequence started
                        if cm > m m = cm;
                        end;
                        cm = 0;
                        fgets(file);
                        n = n + 1;
                    else
                        if isletter(x) || x == '-'
                            cm = cm + 1;
                        end;
                    end;
                end

                if cm > m m = cm;
                end;

                % go throught the file
                Ss = char(m);
                S = [];
                str = zeros(1, MAXNAME);
                sizes = zeros(1, n);
                frewind(file);
                % names=[];
                names = {};
                i = 0;
                j = 1;
                while 1
                    [x, nr] = fscanf(file, '%c', 1);
                    if nr == 0, break;
                    end;
                    if x == '>' % new sequence started
                        if i ~= 0 % save the sequence
                            [x, sizes(i)] = size(Ss);
                            S = strvcat(S, Ss);
                            Ss = [];
                            Ss = char(m);
                        end;
                        str = fgetl(file); % read the name, we remove the '>' symbol
                        % names=strvcat(names,str);
                        pos = find(str == ' ');
                        if ~(isempty(pos))
                            str = str(1:pos(1, 1));
                        end
                        i = i + 1;
                        names{i} = str;
                        disp(str);
                        j = 1;
                    else
                        if isletter(x) || x == '-'
                            % processing the sequence symbol
                            Ss(j) = upper(x);
                            j = j + 1;
                        end;
                    end;
                end

                S = strvcat(S, Ss);
                [x, sizes(i)] = size(Ss);
                fclose(file);

                NT = 'ACGT-';
                [n, m] = size(S);
                seq = ones(n, m) .* 5;
                for k = 1:5, seq(S == NT(k)) = k; end

                aln.seqnames = names;
                aln.seq = seq;
                %
                aln.locus = ones(n, 1);
                aln.population = ones(n, 1);
                aln.count = ones(n, 1);
