function [ldinfo] = snp_ldbyhaploview(genodata, gmarkinfo)
%[ldinfo]=snp_ldbyhaploview(genodata,gmarkinfo)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 2
    gmarkinfo = [];
end
haploviewrun(genodata, gmarkinfo, method);

%{
oldpath=pwd;
cdpge;
cd 'addins/Haploview';
try
    disp('Writing input files (input.ped and input.map)...')
    i_writelinkage(genodata,gmarkinfo,'input');
    a='java -jar Haploview.jar -n -pedfile input.ped -info input.map -skipcheck -blockoutput GAB -dprime';
    disp('Running java -jar Haploview.jar...')
    system(a);
catch
    ldinfo=[];
    return;
end
%}

disp('Reading Haploview output...')
[i, j, dprime, lod, rsqure, cilow, cihi, dist, tint] = textread('input.ped.LD', '%s%s%f%f%f%f%f%d%f', 'headerlines', 1);
ldinfo = struct;
%ldinfo.d = triu(squareform(dv'));
%ldinfo.dprime = triu(squareform(dpv'));
%ldinfo.r2= triu(squareform(r2v'));

% cd(oldpath);


    function [status] = i_writelinkage(geno, mark, filename)
        %saves as linkage format
        %snp_writelinkage(geno,mark,filename)

        if (isempty(geno)), status = 0;
            return;
        end
        if nargin < 2
            mark = [];
        end
        %if (isempty(mark)), status=0; return; end
        if (nargin < 3),
            error('need file name');
        else
            filenamemarker = [filename, '.map'];
            filename = [filename, '.ped'];
        end

        fid = fopen(filename, 'wt');
        if (fid == -1),
            status = 0;
            warning('Unable to open file.');
            return;
        end
        [samplen, marklen] = snp_samplen(geno);
        indvlen = samplen / 2;

        ACGT = '12340';
        for (k = 1:indvlen),
            fprintf(fid, ['Fami%d\tIndv%d\t0\t0\t1\t0\t'], ...
                k, k);

            for (j = 1:marklen * 2),
                fprintf(fid, '%s\t', ACGT(geno(k, j)));
            end
            fprintf(fid, '\n');
        end
        fclose(fid);

        if (~isempty(mark))
            fid = fopen(filenamemarker, 'wt');
            %fprintf(fid,'MARKER_ID\tSNP_rs#\tbp_POSITION\n');
            for k = 1:marklen
                %fprintf(fid,'%d\t%s\t',k,mark.rsid{k});
                %fprintf(fid,'%d\t',k);
                fprintf(fid, '%s\t', mark.rsid{k});
                fprintf(fid, '%s\n', num2str(mark.pos(k)));
                %fprintf(fid,'%s\n',num2str(k));
            end
            fclose(fid);
        end
        status = 1;
