function [d] = tajima89d_test(aln, fromS)
%TAJIMA89D_TEST - Tajima's Test of Neutrality
%This function conducts Tajima's test of neutrality (Tajima 1989). Tajima's test
%compares the number of segregating sites per site with the nucleotide diversity.
%(A site is considered a segregating site if there are two or more nucleotides
%at that site in a comparison of m sequences, and nucleotide diversity is
%defined as the average number of nucleotide differences per site between two
%sequences.)  If all the alleles are selectively neutral, then the product 4Nv
%(where N is the effective population size and v is the mutation rate per site)
%can be estimated in two ways, and the difference in the estimate obtained in
%these two ways provides indication of non-neutral evolution.
%
%Minimum number of sequences in data files: The data file must contain at least
%four sequences. Sites containing alignment gaps (or sites with missing data)
%are not used (these sites are completely excluded).
%
% Syntax: [d]=tajima89d_test(aln)
%
% Inputs:
%    aln   - Alignment structure
%
% Outputs:
%    d     - Tajima's D statistic
%
% See also: tajima89d

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


if nargin < 2, fromS = 0; end
% estimate theta-W from S instead of m_mut.

if (isstruct(aln)), seq = aln.seq;
else seq = aln;
end

[n, m] = size(seq);
if n < 4, error('Four or more sequences are need to compute Tajima''s statistics.'); end

[pi_value] = thetapi(aln, 0, 0); % [pi_value] = estimatetheta(aln,1); or pi_value=theta_pi(aln);
[S, V, m_num, sn_num, sm_num] = countsegregatingsites(aln);

%[sfv] = sfs(aln);
%x=[1:length(sfv)];
%m_num=sum(x.*sfv);
%S=sum(sfv); V=m;

if (fromS), Sk = S;
else Sk = m_num;
end
d = [];

if ~(Sk == 0 && pi_value == 0)
    d = tajima89d(n, Sk, pi_value);
end


if (nargout < 1),
    i_dispheader('Tajima''s Neutrality Test')
    disp('Mode: Nucleotides');
    disp('Gaps/Missing data: Complete Deletion');

    fprintf('No. of Sequences (n): %d\n', n);
    fprintf('No. of Sites (m): %d\n', m);
    fprintf('\n');
    fprintf('No. of Segregating sites (S): %d\n', S);
    fprintf('Total number of mutations, Eta: %d\n', m_num);
    fprintf('\n');
    % fprintf(['Segregating sites per site (pS=S/m) : %f\n'],S/M);
    % Segregating sites per site (pS=S/n)
    fprintf('Average number of nucleotide differences: %f\n', pi_value); % scaled mutation rate theta (pi_value) pi_value=pi_value
    fprintf('Nucleotide diversity (per site), Pi: %f\n', nucdiv(aln));
    fprintf('\n');

    nx = 1:(n - 1);
    a1 = sum(1./nx);
    if (fromS),
        thetaw = S / a1; % Theta (per gene) from S: 5.45455 = wtheta
        thetaw2 = thetaw / V; % Theta (per site) from S: 0.11858
        fprintf('Theta (per gene) from S: %f\n', thetaw);
        fprintf('Theta (per site) from S: %f\n', thetaw2);
    else
        thetaw = m_num / a1; % Theta (per gene) from S: 5.45455 = wtheta
        thetaw2 = thetaw / V; % Theta (per site) from S: 0.11858
        fprintf('Theta (per gene) from Eta: %f\n', thetaw);
        fprintf('Theta (per site) from Eta: %f\n', thetaw2);
    end
    fprintf('\n');

    if isempty(d)
        fprintf('Tajima''s D: -NULL-\n');
    else
        %fprintf (['Diff = %f, s.e. = %f\n'], Diff, DiffSE);
        fprintf('Tajima''s D: %f\n', d);
        %[D] = tajima89d_simu(n,1000,thetaw,0);
        if pi_value > 0
            [bstat] = tajima89d_simu(n, 10000, pi_value, 0);
            p = 2 .* sum(bstat > abs(d)) ./ 10000;
            fprintf('Statistical significance:\n P = %f%s\n', p, sigtag(p));
        end
    end

    i_dispfooter
end
