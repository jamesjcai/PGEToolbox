function [chrid, pos] = snp_locator(rsidstr, showit)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-04-26 23:08:19 -0500 (Fri, 26 Apr 2013) $
% $LastChangedRevision: 532 $
% $LastChangedBy: jcai $

%if strcmp(rsidstr,'rs1637337')
%   chrid=7;
%   pos=74834162;
%   return;
%end

chrid = 0;
pos = 0;

if (nargin < 2)
    showit = 0;
end
if (nargin < 1),
    prompt = {sprintf('This command will return the chromosomal position of a SNP.\n Please enter a SNP ID: ')};
    def = {'rs2693665'};
    dlgTitle = 'Input for SNP ID';
        lineNo = 1;
        answer = inputdlg(prompt, dlgTitle, lineNo, def);

        if ~(isempty(answer)),
            rsidstr = answer{1};
        else
            return;
        end
    end

    if isstr(rsidstr)
        rsid = str2num(rsidstr(3:end));
    end

    urlFetch = sprintf('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&report=DocSet&id=%d', ...
        rsid);
    urlFetch
    try
        pagecontent = urlread(urlFetch);
    catch ME
        %errordlg(lasterr)
        disp(urlFetch)
        error(ME.message);
    end


    fetchResults = char(strread(pagecontent, '%s', 'delimiter', '\n', 'whitespace', ''));
    fetchResults = cellstr(fetchResults);

    numLines = strmatch('CHROMOSOME BASE POSITION=', fetchResults);
    if ~(isempty(numLines))
        theline = fetchResults(numLines, :);
        theline = theline{1};
        %[mat, idx] = regexp(theline,'\d','match','start');    % matlab 7 only
        %[mat2, idx2] = regexp(theline,':','match','start');   % matlab 7 only

        theline = strrep(theline, 'X', '23');
        theline = strrep(theline, 'Y', '24');

        %CHROMOSOME BASE POSITION=1:16818324|1:17153164
        %CHROMOSOME BASE POSITION=1:16818324

        % [mat, idx] = regexp(theline,'\d');
        % [mat2, idx2] = regexp(theline,':');
        % chridstr=theline(idx(find(idx<idx2)));
        [mat, idx] = regexp(theline, '=\d+:');
        chridstr = theline(mat+1:idx-1);
        %    if upper(chridstr)=='X'
        %        chridstr='23';
        %    elseif upper(chridstr)=='Y'
        %        chridstr='24';
        %    end
        chrid = str2num(chridstr);
        [mat2, idx2] = regexp(theline, ':\d+');
        %pos=str2num(theline(idx(find(idx>idx2))))+1;
        %    if length(mat2)>1,
        %        warning('The SNP appears at multiple locations');
        %    end
        %pos=str2num(theline(mat2(1)+1:idx2(1)))+1;
        pos = str2num(theline(mat2(1)+1:idx2(1)));
    end

    if (nargout < 1 || showit)
        i_dispheader('SNP Locator Result')
        %fprintf('Genome Build: 36.3\n');
        %fprintf('Genome Build: 37.1/GRCh37 (as of 6/21/2010)\n');
        fprintf('Genome Build: GRCh38.p2 (as of Aug 1 2015\n');
        fprintf('SNP: %s\n', rsidstr);
        if chrid > 0
            fprintf('Chromosome: %d\n', chrid);
            fprintf('Position: %d\n', pos);
        else
            fprintf('cannot be found in dbSNP.\n');
        end
        i_dispfooter
    end
