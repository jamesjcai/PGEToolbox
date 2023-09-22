%%PGEToolbox DEMO
%
% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


clc;
disp(' ')
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%                            %%')
disp('%% Welcome to PGEToolbox DEMO %%')
disp('%%                            %%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
disp('(1) DEMO 1 - Neutrality Tests')
disp('(2) DEMO 2 - SNP analysis')
disp('(3) EXIT')
disp(' ')
pickI = input('Please select: ');

switch (pickI)
    case (1)
        playshow pge_demo1;
    case (2)
        playshow pge_demo2;
end