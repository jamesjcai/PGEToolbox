function [Ps, Pn] = pspn(aln, icode)
%PSPN - Estimates Ps and Pn

%Ps - Number of synonymous polymorphisms
%Pn - Number of nonsynonymous polymorphisms
%see also: DnDs, LnLs

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin < 2), icode = 1; end
if (isstruct(aln)),
    seq = aln.seq;
else
    seq = aln;
end

m = size(seq, 2);
[ns, na] = getsynnonsyndiff(icode);
[cseq] = codonise64(seq);
%nst=zeros(n); nat=zeros(n);

cm = size(cseq, 2);

Ps = 0;
Pn = 0;
for k = 1:cm
    csite = cseq(:, k);
    [numHap, ~, seqHap] = counthaplotype(csite);
    if (numHap > 1),
        [cps, cpn] = i_codonsynnonsyn(seqHap, ns, na);
        %[s_site,v_site,m_num,sn_num,sm_num,sm_numv]=countsegregatingsites(decodonise64(seqHap));
        [m_num] = i_mutnum(decodonise64(seqHap));
        %fprintf ('cps %d *(m_num %d )/(cps %d +cpn %d)\n', cps, m_num, cps, cpn);

        if ((cps + cpn) == 0)
            xps = 0;
            xpn = 0;
        else
            xps = cps * m_num / (cps + cpn);
            xpn = cpn * m_num / (cps + cpn);
        end
        Ps = Ps + xps;
        Pn = Pn + xpn;
    end
end


%%%%%%%%%%%%%
%%% SUBS  %%%
%%%%%%%%%%%%%

    function [cps, cpn] = i_codonsynnonsyn(csite, ns, na)
        [n, m] = size(csite);
        cps = 0;
        cpn = 0;

        for i = 1:n - 1
            for j = i + 1:n
                x = csite(i);
                y = csite(j);
                xx = ns(x, y);
                yy = na(x, y);
                if (xx > 0), cps = cps + xx; end
                if (yy > 0), cpn = cpn + yy; end
            end
        end

            function [m_num] = i_mutnum(seq)
                %e.g., seq=[1 3 4; 1 4 4; 1 3 3]
                m_num = 0;
                for k = 1:3
                    thissite = seq(:, k);
                    x = length(unique(thissite)) - 1;
                    m_num = m_num + x;
                end
