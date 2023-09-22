function [zz] = coolloci

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

ColNames = {'Locus', 'Position', 'Focal Point', 'PMID'};
ColWidth = [100, 200, 90, 90];

zz = cell(2, 3);
zz(1, 1) = {'KITLG'};
zz(1, 2) = {'Chr12:86000000..89000000'};
zz(1, 3) = {'87450000'};
zz(1, 4) = {'17542651'};

zz(2, 1) = {'MLPH'};
zz(2, 2) = {'Chr2:237900000..238200000'};
zz(2, 3) = {'238066523'};
zz(2, 4) = {'20416085'};

zz(3, 1) = {'ATRNL1'};
zz(3, 2) = {'Chr2:116600000..118000000'};
zz(3, 3) = {'117667124'};
zz(3, 4) = {'20416085'};

zz(4, 1) = {'VDR'};
zz(4, 2) = {'Chr12:46450000..46650000'};
zz(4, 3) = {'46529644'};
zz(4, 4) = {'20416085'};

zz(5, 1) = {'RXRA'};
zz(5, 2) = {'Chr9:136300000..136500000'};
zz(5, 3) = {'136377116'};
zz(5, 4) = {'20416085'};

zz(6, 1) = {'GIP'};
zz(6, 2) = {'Chr17:44390917..44400954'};
zz(6, 3) = {'44390917'};
zz(6, 4) = {'999'};


%if nargout<1
zz = tableGUI('FigName', 'Cool Loci', 'array', zz, 'ColNames', ColNames, 'ColWidth', ColWidth, 'NumRows', 2, ...
    'RowNumbers', 'y', 'checks', {0, 1, 0});
%end
