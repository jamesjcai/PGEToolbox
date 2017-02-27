function [p]=brunnermunzeltest(x,y,alpha)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if nargin<3
    alpha=0.05;
end
n1=length(x);
n2=length(y);
r1=tiedrank(x);
r2=tiedrank(y);
r=tiedrank([x;y]);

m1=mean(r(1:n1));
m2=mean(r((n1+1):end));

% estimate P(X<Y)+.5*P(X=Y) from data
pst=(m2-(n2+1)/2)/n1;


% Calculate the Empirical Variance of the first and second sample ###
v1=sum((r(1:n1)-r1-m1+(n1+1)./2).^2)./(n1-1);
v2=sum((r(n1+1:end)-r2-m2+(n2+1)./2).^2)./(n2-1);

% Calculate the T statistic and degree of freedom of the whole data set ###
statistic=n1*n2*(m2-m1)/(n1+n2)/sqrt(n1*v1+n2*v2);
dfbm=((n1*v1+n2*v2)^2)/(((n1*v1)^2)/(n1-1)+((n2*v2)^2)/(n2-1));


%### Users can choose alternative hypothesis ###
%### Must be one of "two.sided"(default), "greater" or "less" ###
%### Users can specify just the initial letter ###

alternative='two-sided';
switch alternative
    case{'g','greater'}
        p=1-tcdf(abs(statistic),dfbm);        
    case{'l','less'}
        p=tcdf(abs(statistic),dfbm);        
    case{'t','two-sided'}
        statistic
        dfbm
        p=tcdf(statistic,dfbm)
        p=2*min([1-p,p])
end

%pL=tpdf(abs(statistic),dfbm);
%pR=1-pL;
%p=2*min([pR,pL]);

%## Calculate Confidence Interval
%conf.int=c(pst-qt(1-alpha/2,dfbm)*sqrt(v1/(n1*n2^2)+v2/(n2*n1^2)),...
%    pst+qt(1-alpha/2,dfbm)*sqrt(v1/(n1*n2^2)+v2/(n2*n1^2)))

%1-alpha/2
%pst+tinv(alpha/2,dfbm)



%### Display Output ###
%estimate=pst
%ESTIMATE=pst
%names(ESTIMATE) = "P(X<Y)+.5*P(X=Y)"
%%STATISTIC=statistic
%names(STATISTIC) = "Brunner-Munzel Test Statistic"
%PARAMETER = dfbm
%names(PARAMETER) = "df"
%CONF.INT = conf.int
%names(CONF.INT)=c("lower","upper")
%attr(CONF.INT,"conf.level")=(1-alpha)
%METHOD = "Brunner-Munzel Test"
i_dispheader('Brunner-Munzel Test');
fprintf('Statistic = %f\n', statistic);
fprintf('df = %f\n',dfbm);
%fprintf('p-value = %f (%s)\n',pL,'Left tail');
%fprintf('p-value = %f (%s)\n',pR,'Right tail');
fprintf('p-value = %d (%s)\n',p,alternative);
fprintf('sample estimates:\n');
fprintf('P(X<Y)+.5*P(X=Y)\n');
fprintf('        %f\n',pst); 
i_dispfooter;

