function [p,dprs] = mktestcmd(aln)
%MKTESTCMD - McDonald-Kreitman test
%USAGE: [p,dprs] = mktestcmd(aln)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[Ds,Ps,Dn,Pn] = dspsdnpn(aln);
[p]=fisherextest(round(Ds),round(Ps),round(Dn),round(Pn));
if (nargout>1),
	dprs=[Ds,Ps;Dn,Pn];
end

if (nargout<1),
i_dispheader('McDonald-Kreitman''s Test')
%fprintf ([' Syn. Subs.: \n']);
%fprintf (['     Fixed differences between species:    %2.1f   Polymorphic sites:    %2.1f\n'],Ds,Ps);
%fprintf ([' NSyn. Subs.:  \n']);
%fprintf (['     Fixed differences between species:    %2.1f   Polymorphic sites:   %2.1f\n'],Dn,Pn);

fprintf ('Ds: %5.1f\t\tPs:%5.1f\n',Ds,Ps);
fprintf ('Dn: %5.1f\t\tPn:%5.1f\n\n',Dn,Pn);
fprintf ('Ds - Syn. Divergences; Ps - Syn. Polymorphisms\n');
fprintf ('Dn - Nonsyn. Divergences; Pn - Nonsyn. Polymorphisms\n');
fprintf ('\n');
fprintf ('Fisher''s exact test\n P-value (two tailed): %f%s\n',p, sigtag(p));
[P,G]=gtest(Ds,Ps,Dn,Pn);
fprintf(['G Test, G value: %10.3f\n'], G);
fprintf(['   P-value: %10.5f %s\n\n'], P, sigtag(P));

	if (all([Ps, Dn, Ds])),
		%NI = (Ds/Ps)/(Dn/Pn);
        NI = (Pn/Dn)/(Ps/Ds);
        alphax = 1-(Pn*Ds)/(Ps*Dn);
	else
	    NI=0;
        alphax = NaN;
	end
	fprintf ('Neutrality Index (NI), [(Pn/Dn)/(Ps/Ds)]: %f\n',NI);
	fprintf ('Alpha [1-(Pn*Ds)/(Ps*Dn)]: %f\n',alphax);
i_dispfooter
end

% NI = (Pn/Dn)/(Ps/Ds); in original paper (NI; Rand and Kann 1996, 1998)
% NI<1 indicates positive selection
% the ratio of nonsynonymous polymorphisms to fixations divided by the
% ratio of synonymous polymorphisms to fixations

% Neutrality Index, NI: = (a/b)/(c/d);
% A simple test of the null hypothesis is a Chi-square test, which is
% X2 =n(ad-bc)2 /[(a+b)(a+c)(b+d)(c+d)] where n=a+b+c+d is the total
% number of polymorphic sites. When n is not small, X2 follows
% approximately Chi-square distribution with one degree of freedom. So
% if the value of X2 is larger that 3.841, the null hypothesis can be
% rejected at 5%




%Input Data File: C:\...\seq_examples\mktest.fas
% Intraspecific Data: aaa
%   Number of sequences: 4
% Interspecific Data: bbb
%   Number of sequences: 4
%
% Selected region: 1-27     Number of sites: 27
% Total number of codons: 9     Total number of codons analyzed: 9
%
% Genetic Code: Nuclear Universal
% Protein Coding, and Noncoding Regions:
%   Number of protein coding regions (exons): 1
%   Number of noncoding regions (intronic and flanking regions): 0
%     Protein coding region, from site: 1  to  27
%
% Some codons with multiple evolutionary paths were detected.
%    The number of synonymous and nonsynonymous mutations were estimated
%    using the DnaSP conservative criteria (see the help file).
%      Sites affected:  9 15 18 21 22 24 25 27
%
% ======= Polymorphic Changes Within Species 1 =======
% Segregating sites: 16   Total number of mutations: 19
% Total number of Synonymous changes: 7
%           3 3 9 15 18 24 27
% Total number of Nonsynonymous changes: 12
%           3 7 11 12 14 14 16 21 22 23 25 26
%
% ======= Polymorphic Changes Within Species 2 =======
% Segregating sites: 3   Total number of mutations: 3
% Total number of Synonymous changes: 2
%           21 27
% Total number of Nonsynonymous changes: 1
%           11
%
% ======== Polymorphic Changes Within Species ========
% Segregating sites: 16   Total number of mutations: 19
% Total number of Synonymous changes: 7
%           3 3 9 15 18 24 27
% Total number of Nonsynonymous changes: 12
%           3 7 11 12 14 14 16 21 22 23 25 26
%
% ======= Fixed Differences Between Species ========
% Total number of Synonymous changes: 1
%           15
% Total number of Nonsynonymous changes: 4
%           4 13 18 19
%
% ============== McDonald and Kreitman Table ==============
% Synonymous Substitutions:
%     Fixed differences between species:    1   Polymorphic sites:    7
% Nonsynonymous Substitutions:
%     Fixed differences between species:    4   Polymorphic sites:   12
%
% Neutrality Index, NI: 0.429
%
% Fisher's exact test. P-value (two tailed): 0.631094  (not significant)
%
% G test. G value: 0.540      P-value: 0.46224  (not significant)
% G value with Williams' correction: 0.481
%      P-value: 0.48786  (not significant)
% G value with Yates' correction: 0.032
%      P-value: 0.85804  (not significant)
%
% * 0.01<P<0.05;  ** 0.001<P<0.01;  *** P<0.001