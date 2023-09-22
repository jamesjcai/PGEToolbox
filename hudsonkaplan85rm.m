function [rmin, D] = hudsonkaplan85rm(seq, showit)
%HUDSONKAPLAN85RM - Minimum number of recombination events (Rm)
%
% aka: Four-gamete test
%
%see also: myersgriffiths03r

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin < 2
    showit = 0;
end

[haplodata] = extractsegregatingsites(seq, 1);
[m] = size(haplodata, 2);
Dcount = zeros(m);

for i = 1:m - 1
    for j = i + 1:m
        %[numHap]=counthaplotype([haplodata(:,i),haplodata(:,j)]);
        numHap = size(unique([haplodata(:, i), haplodata(:, j)], 'rows'), 1);
        Dcount(i, j) = numHap;
    end
end
%Dcount

D = (Dcount == 4);


[x, y] = find(triu(D));

if isempty(x)
    rmin = 0;
else
    idx = zeros(1, length(x)+1);
    x = [x; 0];
    y = [y; 0];
    idx(end) = 1;
    while any(idx)
        id = find(idx);
        x(id) = [];
        y(id) = [];
        idx(id) = [];
        if ~isempty(x)
            for k = 2:length(x)
                if x(k) >= x(k-1) && y(k) <= y(k-1)
                    idx(k-1) = 1;
                end
            end
        end
    end

    if length(x) > 1
        x = [x(end:-1:1); 0];
        y = [y(end:-1:1); 0];
        idx(end) = 1;
        while any(idx)
            id = find(idx);
            x(id) = [];
            y(id) = [];
            idx(id) = [];
            if ~isempty(x)
                for k = 2:length(x)
                    if x(k) >= x(k-1) && y(k) <= y(k-1)
                        idx(k-1) = 1;
                    end
                end
            end
        end
    end
    rmin = length(x);
end

if (nargout < 1 || showit)
    i_dispheader('Four-gamete test (Hudson and Kaplan 1985)')
    fprintf('Minimum number of recombination events\n');
    fprintf('   Number of biallelic segregating sites: %d\n', m);
    fprintf('   Number of pairwise comparisons analyzed: %d\n', nchoosek(m, 2));
    fprintf('   Number of pairs of sites with four gametic types: %d\n', sum(D(:)));
    fprintf('   Minimum number of recombination events, Rm: %d\n', rmin);
    i_dispfooter
end


%d=Array of pairwise 4GT test results
%d=0,   if  there are less then 4 gametes
%d=1,   if  there are 4 gametes

%What is the minimal number of recombinations that could
%explain observed data ?
%Statistics FR    (Hudson and Kaplan, 1985)


%Wang et al., 2002 - Study
%R. Hudson’s program for simulating genealogies with mutation, drift and
% recombination under various demographic scenarios
%Study of dependence of average lengths of blocks on different factors
%Comparison of simulation results to data from Patil et al., 2002


%R_M, denotes the minimum number of recombination events implied by
%the data using the four-gamete test
%R, denotes the total number of recombination events in the history of
%the sample.

% We describe here an algorithm for determining RM. Recall the matrix D =
% (d(i, j)), where d(i, j ) = 1 if all four gametes involving sites i and j
% are present in the sample; otherwise d(i, j) = 0. To each nonzero element
% d(i, j) of D above the diagonal we associate the open interval (i, j) and
% form a list of these intervals ordered so that the starting points of the
% intervals are not decreasing. The method for finding RM deletes certain
% members of this list. The first type of intervals deleted are those that
% completely contain other intervals. For example, if (i, j ) and (m, n)
% are on the list, and m <= i < j <= n then (m, n) is deleted. For the
% remaining intervals on the list, let (i2, j2) be the first interval with
% i2>=j1 such that (i2,j2) is not disjoint from all of the remaining
% intervals and delete all intervals whose first component is <j2 and >i2.
% This process is continued until it is not possible to find an interval
% that has a nonempty intersection with some other interval on the list. At
% this point all of the intervals are disjoint, and to each interval we
% must assign at least one recombination event. Hence, RM equals the final
% number of intervals on the list.
