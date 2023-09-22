function [centromere] = centromerecoord(hg)
%CENTROMERECOORD - coordinates of centromeres of human chromosomes

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2014-03-26 09:30:26 -0500 (Wed, 26 Mar 2014) $
% $LastChangedRevision: 758 $
% $LastChangedBy: jcai $

if nargin < 1
    hg = 19;
end

if hg == 19
    centromere = [121535434, 124535434; ...
        92326171, 95326171; ...
        90504854, 93504854; ...
        49660117, 52660117; ...
        46405641, 49405641; ...
        58830166, 61830166; ...
        58054331, 61054331; ...
        43838887, 46838887; ...
        47367679, 50367679; ...
        39254935, 42254935; ...
        51644205, 54644205; ...
        34856694, 37856694; ...
        16000000, 19000000; ...
        16000000, 19000000; ...
        17000000, 20000000; ...
        35335801, 38335801; ...
        22263006, 25263006; ...
        15460898, 18460898; ...
        24681782, 27681782; ...
        26369569, 29369569; ...
        11288129, 14288129; ...
        13000000, 16000000; ...
        58632012, 61632012; ...
        10104553, 13104553];

elseif hg == 18

    centromere = [121236957, 123476957; ...
        91689898, 94689898; ...
        90587544, 93487544; ...
        49354874, 52354874; ...
        46441398, 49441398; ...
        58938125, 61938125; ...
        58058273, 61058273; ...
        43958052, 46958052; ...
        47107499, 50107499; ...
        39244941, 41624941; ...
        51450781, 54450781; ...
        34747961, 36142961; ...
        16000000, 17868000; ...
        15070000, 18070000; ...
        15260000, 18260000; ...
        35143302, 36943302; ...
        22187133, 22287133; ...
        15400898, 16764896; ...
        26923622, 29923622; ...
        26267569, 28033230; ...
        10260000, 13260000; ...
        11330000, 14330000; ...
        58598737, 61598737; ...
        11253954, 11653954];

end

%M10=10000000;
%cenA=centromere(chrid,1)-M10;
%cenB=centromere(chrid,2)+M10;

% SOURCE:
%http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz
%http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
