function [res, data] = line4mkrpf(aln)
%LINE4MKRPF - writes input line for MKRPF test

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[Ds, Ps, Dn, Pn, Ls, Ln] = dspsdnpn(aln);
disp('1. class name')
disp('2. gene name')
disp('3. FS : # fixed silent silents, Ds')
disp('4. SS : # segregating silent sites, Ps')
disp('5. FR : # fixed replacement sites, Dn')
disp('6. SR : # segregating replacement sites, Pn')
disp('7. n1 : number of sequences sampled in species 1.')
disp('8. n2 : number of sequences sampled in species 2.')
disp('9. TS : total number of silent sites in the alignment.')
disp('10. TR : total number of replacement sites in the alignment.')
disp('11. Ratio : haploid ratio (1.0 for autosomal, 0.75 for X-linked, 0.25 for Y-linked).')

            res = sprintf('ClassA\tLocusA\t%f\t%f\t%f\t%f\t%d\t%d\t%f\t%f\t1', ...
                Ds, Ps, Dn, Pn, 90, 1, Ls, Ln);
            if (nargout < 1)
                disp(res)
            end
            if (nargout > 1)
                data = [Ds, Ps, Dn, Pn, 90, 1, Ls, Ln];
            end
