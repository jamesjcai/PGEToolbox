

%
%http://mbe.oxfordjournals.org/cgi/reprint/24/7/1562.pdf
%
% The symmetric finite island model (Wright 1931)
%with 2 or 3 subpopulations (demes). The level of differentiation
%is measured by the FST statistic, which can be calculated
%as F_{ST}
%
%
% Ns = the number of breeding individuals per deme, 
% m = the probability that each gene is an emigrant, 
% d = the number of demes (d = 2 or 3) (e.g., Slatkin 1991).

d=d.^2/(d-1).^2;
Fst = 1./(1+4*Ns.*m.*d);



% Nm = Ns.*m = number of migratants per generation
Fst = 1./(1+4*Ns.*m);     % Wright (1931)
Nm = (1-Fst)./4*Fst;

%When Fst=0.2, Nm=1
%
%Slatkin M. 1991. Inbreeding coefficients and coalescence times.
%Genet Res. 58:167?175.



%%
%
% Sampling theory of neutral mutations by Ewens (1972)
%
% The expected number of alleles per locus (a) in a sample of size n would
% be  


b=0:(2*n);
a = sum(theta./(theta+b-1));



a/sum(1./1:(n-1));




