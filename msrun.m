function [G, pos] = msrun(nsam, nrep, theta, segs, rho, nsites)
%MSRUN - executes Hudson's MS program to do coalescent simulations
%
%[G] = msrun(nsam,nrep,theta,segs,rho,nsites)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin < 6), nsites = 2; end
if (nargin < 5), rho = 0; end
if (nargin < 4), segs = 0; end
if (nargin < 3),
    error('MSRUN:MissingInput', ...
        'MSRUN requires at least 3 arguments, e.g., [G,pos]=msrun(6,10,12.5);');
    end

    seedms_mex;
    if (nargout == 1)
        G = cell(nrep, 1);
        for k = 1:nrep,
            try
                [io] = ms_mex(nsam, theta, segs, rho, nsites);
            catch exception
                io = [];
                error(exception.message);
            end
            G{k} = io;
        end
    elseif (nargout == 2)
        G = cell(nrep, 1);
        pos = cell(nrep, 1);
        for k = 1:nrep,
            try
                [io, posi] = ms_mex(nsam, theta, segs, rho, nsites);
            catch exception
                error(exception.message);
                io = [];
                posi = [];
            end
            G{k} = io;
            pos{k} = posi;
        end
    end


    %http://mbe.oxfordjournals.org/cgi/content/full/18/6/1134

    %Standard coalescent simulations first produce random genealogies, then place
    %mutations at constant rate {theta}/2 ({theta} = 4Nµ is the population mutation
    %parameter, where N is the effective population size and µ is the per-locus
    %mutation rate per generation) along each of the branches (Kingman 1982a,
    %1982b; Hudson 1990).

    % Builds a random tree every time next_tree is called or up to -maxcount
    % times with branch lengths and provides the ability to randomly add
    % mutations onto the tree with a probabilty proportional to the branch
    % lengths.
    %
    % This algorithm is based on the make_tree algorithm from Richard Hudson 1990.
    %
    % Hudson, R. R. 1990. Gene genealogies and the coalescent
    %        process. Pp. 1-44 in D. Futuyma and J.  Antonovics, eds. Oxford
    %        surveys in evolutionary biology. Vol. 7. Oxford University
    %        Press, New York.
    %
    % This module was previously named Bio::Tree::RandomTree


    % The long term evolution of DNA can be modelled using the coalescent. The Markov
    % chain simulation technique developed by Griffiths and Tavaré can be used to
    % produce likelihood estimates for DNA sequence data under the
    % infinitely-many-sites model. The program GENETREE will calculate not only
    % likelihoods but also simulate likelihood surfaces for parameter estimationin
    % subdivided populations with or without population growth. A whole host of
    % statistics can be calculated including times to the most recent common ancestor,
    % probabilities of the location of the most recent common ancestors, times to the
    % different mutations and the probabilities of these having occurred in the
    % various subpopulations.