function snp_plotmarkpos(markinfo, showchr)
%SNP_PLOTMARKPOS - plots position of SNP markers
%snp_plotmarkpos(markinfo,showchr)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

% Consider the following simplification:
%h=stem(markinfo.pos,markinfo.maf)
%set(h,'Marker','none')

if (nargin < 2), showchr = false; end
if (showchr)
    if (isfield(markinfo, 'chr'))
        chrid = markinfo.chr{1}; % chrid='chr13'
        chridx = chrid(4:end);
        if strcmpi('X', chridx)
            chrid = 23;
        elseif strcmpi('Y', chridx)
            chrid = 24;
        else
            chrid = str2num(chridx);
        end
    else
        error('Unkonwn chromosome.')
    end
    subplot(2, 1, 1)
    ideogram(chrid);
    title(sprintf('Chromosome %s', chridx))
    subplot(2, 1, 2)
end

pos = markinfo.pos;
nb = length(pos); %  number of boxplots
if (isfield(markinfo, 'maf'))
    mafp = markinfo.maf;
else
    mafp = 0.5 * ones(1, nb);
end
if (isfield(markinfo, 'daf'))
    dafp = markinfo.daf;
else
    dafp = [];
end
hold on
for ii = 1:nb
    % xx = [pos(ii) pos(ii)]./1000000;
    xx = [pos(ii), pos(ii)];
    %yy = [1 1+mafp(ii)];
    if isempty(dafp)
        yy = [0, mafp(ii)];
    else
        yy = [0, dafp(ii)];
    end
    plot(xx, yy, '-')
end
xlabel('Position')
ylabel('Allele Freq.')
%plot([min(min(pos)) max(max(pos))],[1 1],'-r')


if (showchr)
    %  make some extra space
    axlim = axis;
    %axlim(1) = axlim(1)-1;
    %axlim(2) = axlim(2)+1;

    axlim(1) = 0;
    axlim(2) = chrlen(21);

    %axlim(3) = axlim(3)-2;
    %axlim(4) = axlim(4)+2;
    axis(axlim)
end

%axis off

%i_dispheader('Plot SNP position')
%disp('LEGEND: Green lines indicate relative positions of SNPs')
%i_dispfooter