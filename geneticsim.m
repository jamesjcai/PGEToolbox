function [s] = geneticsim(fst)
%GENETICSIM - a measure of genetic similarity (Rousset 1997) based on Fst
%
% FST/(1-FST) as measure of genetic similarity (Rousset 1997).
%

s = fst ./ (1 - fst);


%ezplot('fst/(1-fst)',0,1)
%ylabel('Genetic Similarity')

%REF:
%Rousset, F 1997. Genetic differentiation and estimation of gene flow from F-statistics under
%isolation by distance. Genetics 145: 1219-1228.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $
