function snp_fstat(geno1,geno2)


r=2;                                    % two popuplations
%hok1=snp_obshet(geno1);
%hok2=snp_obshet(geno2);
[hok1,hv]=snp_heterozygosity(geno1)
[hok2,hv]=snp_heterozygosity(geno2)

n1=snp_samplen(geno1);
n2=snp_samplen(geno2);
n_bar = mean([n1 n2]);                  % average sample size
h_bar = (n1*hok1+n2*hok2)/(r*n_bar)     % average heterozygote freq

maf1=snp_maf(geno1)
maf2=snp_maf(geno2)


%[x,y,z]=counthaplotype([maf1,maf2]');
%p_bar=sum((y.*z')./sum(y))
p_bar = (n1*maf1+n2*maf2)/(r*n_bar);     % average allele freq


s2 = (n1*(maf1-p_bar)^2+n2*(maf2-p_bar)^2)/((r-1)*n_bar);   % sampl var of allel freq over pop
nc = (r*n_bar - (n1^2/(r*n_bar)+n2^2/(r*n_bar)))/(r-1);   % var of sample sizes
nc2 = n_bar*(1-(var([n1 n2])^2)/r);    %
if (nc~=nc2) error('xxx'); end




%sigw = (n1*hok1 + n2*hok2)/(n1 + n2)
siga = (n_bar/nc)*(s2-(1/(n_bar-1))*(p_bar*(1-p_bar)-((r-1)/r)*s2-0.25*h_bar ));
sigb = (n_bar/(n_bar-1))*(p_bar*(1-p_bar)-((r-1)/r)*s2-((2*n_bar-1)/4*n_bar)*h_bar);
sigw = 0.5*h_bar



sigb = (n1+n2)*p_bar*(1-p_bar)-(n1*(maf1-p_bar)^2+n2*(maf2-p_bar)^2)-0.25*(n1*hok1+n2*hok2);
sigb = sigb/(n1+n2-r)-0.25*(n1*hok1+n2*hok2)/(n1+n2)

siga = (n1*(maf1-p_bar)^2+n2*(maf2-p_bar)^2)/(r-1);
siga = siga - ((n1+n2)*p_bar*(1-p_bar)-0.25*(n1*hok1+n2*hok2)-(n1*(maf1-p_bar)^2+n2*(maf2-p_bar)^2)/(r-1))/(n1+n2-r);
siga = siga/(n1+n2 - (n1^2+n2^2)/(n1+n2))


siga = 0.316;      % the among sample variance component (between populations)
sigb = -0.002;     % the between individual within sample component (between indiv withing pop)
sigw = 0.335;      % the within individual component (between gametes within indiv.)



Capf = (siga+sigb)/(siga+sigb+sigw);     % Fit  correlation of genes within indiv ('inbreeding')
theta = siga/(siga+sigb+sigw)        % i.e., Fst  correlation of genes of different idv. in same population ('coancestry')
smallf = sigb/(sigb+sigw);           % f   Fis  correlation of genes within individuals within populations

% (Capf-theta)/(1-theta)==smallf

