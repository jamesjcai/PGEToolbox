function [driftmatrix] = kimuradrift(pmatrix, nloci, npops, numalleles)
%
%pseudo-sampling device proposed by Kimura (1980) to simulate the diffusion
%process

%Code copyright 2002 by M. B. Hamilton, J. R. Miller and Georgetown University.
%
%populations drift, uses Kimura's psuedo random variable method.
%pmatirx is the original matrix of p values for each population (columns)
%and for each locus (rows).  nloci is the number of loci, npops is the number
%of populations, numalleles is the number of alleles sampled (2*N for nuclear alleles
%and N for organelle alleles).

%Here, we used the pseudo-sampling device
%proposed by Kimura (1980) to simulate the diffusion process
%directly (see Griffiths [2003]; Coop and Griffiths
%[2004]; and Spencer and Coop [2004] for another way
%of simulating the conditional diffusion process).

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

driftmatrix = pmatrix;
for k = 1:nloci
    for j = 1:npops %allele frequencies drift

        delta = (2 * rand - 1) * sqrt((3 * pmatrix(k, j) * (1 - pmatrix(k, j)))/(numalleles));
        if (pmatrix(k, j) + delta) > 1.0
            driftmatrix(k, j) = 1.0;
        elseif (pmatrix(k, j) + delta) < 0.0
            driftmatrix(k, j) = 0.0;
        else
            driftmatrix(k, j) = pmatrix(k, j) + delta;
        end

    end %end npops loop
end %end loci loop
