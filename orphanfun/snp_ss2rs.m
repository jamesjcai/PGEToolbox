function [rsid]=snp_ss2rs(ssid)

% http://www.ncbi.nlm.nih.gov/projects/SNP/snp_retrieve.cgi?subsnp_id=ss48407963

%Reference SNP Id(rs#): <a href="snp_ref.cgi?rs=rs750625">rs750625</a>

urlFetch=sprintf('http://www.ncbi.nlm.nih.gov/projects/SNP/snp_retrieve.cgi?subsnp_id=%s',...
    ssid);
try
    pagecontent=urlread(urlFetch);
catch
    disp(urlFetch)
    rethrow(lasterror);
end
    
fetchResults = char(strread(pagecontent,'%s','delimiter','\n','whitespace',''));
fetchResults = cellstr(fetchResults);

numLines = strmatch('  Reference SNP Id(rs#): <a href="snp_ref.cgi?rs=',fetchResults);
if ~(isempty(numLines))
    theline=fetchResults(numLines,:);
    theline=theline{1};
    %[mat, idx] = regexp(theline,'\d','match','start');    % matlab 7 only
    %[mat2, idx2] = regexp(theline,':','match','start');   % matlab 7 only
    [mat]=regexp(theline,'rs[0-9]+','match');
    rsid=mat{1};
else
    rsid='';
end
