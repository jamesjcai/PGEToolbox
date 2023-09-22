function [ihs, ihsori, ehh, markinfo, coreset, plotdata] = snp_ihssingle(rsid, radius, popcode, showit, chrid, pos, inputfile)
%SNP_IHSSINGLE - Integrated haplotype scores around a single SNP.
%
%  Syntax: [ihs,ihsori]=snp_ihssingle(rsid,radius,popcode,showit)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-05-12 11:07:32 -0500 (Sun, 12 May 2013) $
% $LastChangedRevision: 538 $
% $LastChangedBy: jcai $

ihs = 0;
ihsori = 0;
if nargin < 7, inputfile = ''; end
if nargin < 5, chrid = []; end
if nargin < 6, pos = []; end
if nargin < 4, showit = 0; end
if nargin < 3
    popcode = 'CEU';
    %popcode='JPT%2BCHB';
else
    popcode = upper(popcode);
    if ~(ismember(popcode, {'CEU', 'CHB', 'JPT', 'YRI', 'JPT+CHB'})),
        error('Not available population code.');
    end
    if strcmp(popcode, 'JPT+CHB')
        popcode = 'JPT%2BCHB';
    end
end

if nargin < 2, radius = 1000000; end
if (nargin < 1),
    prompt = {'RSID'};
    def = {'rs996521'}; %def={'rs6060371'}; %def={'rs7662684'};
    dlgTitle = 'Input RSID';
    lineNo = 1;
    answer = inputdlg(prompt, dlgTitle, lineNo, def);
    if ~(isempty(answer)),
        rsid = answer{1};
    else
        return;
    end
end


if isempty(chrid) || isempty(pos)
    disp(sprintf('Obtaining chromosomal position of SNP %s...', rsid))
    [chrid, pos] = snp_locator(rsid, 1);
end

if ~(chrid > 0 && pos > 0)
    return;
end

startn = max(1, pos-radius);
endn = min(chrlen(chrid), pos+radius-1);


if chrid > 22
    if chrid == 23
        chrid = 'X';
    elseif chrid == 24
        chrid = 'Y';
    else
        error('Wrong chromosome number.')
    end
    query = sprintf('Chr%s:%d..%d', chrid, startn, endn);
else
    query = sprintf('Chr%d:%d..%d', chrid, startn, endn);
end

fprintf('Dowloading phased haplotype data from %s...\n', query);
fprintf('Length of the region = %d bp\n', endn-startn+1);

if ~isempty(inputfile)
    [hapldata, markinfo] = snp_readhaplotype(inputfile);
else
    pause(1)
    [hapldata, markinfo] = snp_downloadhaplotype(query, popcode, true);
    pause(1)
end

[hapldata, markinfo] = hap_pickmarker(hapldata, markinfo, markinfo.maf > 0.015);

chrpos = markinfo.pos;

[yes, idx] = ismember(rsid, markinfo.rsid);
if ~yes
    error('Supplied rsid is not found in the region.')
end

fprintf('Calculating EHH for %s...\n', rsid);
    [ehh, ~, coreset] = snp_ehh(hapldata, idx, idx, 0);

    [b] = i_ehh_leftrightboundries(ehh, idx);
    %ehh=(ehh-0.05+1/72);
    if (size(ehh, 1) == 2),
        ihh1 = trapz(chrpos(b(1, 1):b(1, 2)), ehh(1, b(1, 1):b(1, 2)));
        ihh2 = trapz(chrpos(b(2, 1):b(2, 2)), ehh(2, b(2, 1):b(2, 2)));
        %ihh1=trapz(ehh(1,:));
        %ihh2=trapz(ehh(2,:));
        if ihh2 > 0
            ihsori = log(ihh1/ihh2);
            % ihs=(ihsori-0.7970)/0.6040;
        end
    end

    %disp('iHS scores were standardized empirically with the CEU HapMap data (mean = 0.7970; s.d.= 0.6040)');
    %http://www.nature.com/nature/journal/v449/n7164/extref/nature06250-s1.pdf
    %page 2: " In particular, their (Voight et al) iHHA is actually calculated
    %by integrating the quantity (EHH-0.05+1/N), with N being the number of
    %chromosomes carrying A (personal communication), and similarly for iHHD.
    ihs = ihsori;

    if (nargout < 1 || showit)
        i_dispheader('Integrated Haplotype Score (iHS), unnormalized')
        fprintf('SNP: %s\n', rsid);
        fprintf('Population: %s\n', upper(popcode));
        fprintf('Area under haplotype 1: %f\n', ihh1)
        fprintf('Area under haplotype 2: %f\n', ihh2)
        fprintf('Ratio of two areas: %f\n', abs(ihh1/ihh2))
        %fprintf('|iHS| without noramlization: %f\n', ihsori);
        fprintf('iHS: %f\n', ihs);
        %else
        %	fprintf ('N/A.\n');

        figure;
        plotdata.x1 = chrpos(b(1, 1):b(1, 2));
        plotdata.x2 = chrpos(b(2, 1):b(2, 2));
        plotdata.y1 = ehh(1, b(1, 1):b(1, 2));
        plotdata.y2 = ehh(2, b(2, 1):b(2, 2));
        plotdata.vlinepos = chrpos(idx);
        plotdata.hlinepos = 0.05;

        hold on
        plot(chrpos(b(1, 1):b(1, 2)), ehh(1, b(1, 1):b(1, 2)), '-m+');
        plot(chrpos(b(2, 1):b(2, 2)), ehh(2, b(2, 1):b(2, 2)), '-bx');
        vline(chrpos(idx))
        hline(0.05);
        hold off
        box on
        opengl software
        ylabel('EHH')
        xlabel('SNPs')
        x = 'ACGTN';
        legend({x(coreset(1)), x(coreset(2))})
        plotdata.legend = {x(coreset(1)), x(coreset(2))};
    end
    i_dispfooter;
