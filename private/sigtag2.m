function [tag] = sigtag2(d,bstat)
%SIGTAG2 - 

alpha=0.001;
ci = bootper(bstat,alpha);
%ci = bootcper(bstat,alpha,d)
if (d<=ci(1)|d>=ci(2))
    tag=' (***significant, P<0.001)';
    return
end

alpha=0.01;
ci = bootper(bstat,alpha);
%ci = bootcper(bstat,alpha,d)
if (d<=ci(1)|d>=ci(2))
    tag=' (**significant, 0.001<P<0.01)';
    return
end

alpha=0.05;
ci = bootper(bstat,alpha);
%ci = bootcper(bstat,alpha,d)
if (d<=ci(1)|d>=ci(2))
    tag=' (*significant, 0.01<P<0.05)';
    return
end

alpha=0.1;
ci = bootper(bstat,alpha);
%ci = bootcper(bstat,alpha,d)
if (d<=ci(1)|d>=ci(2))
    tag=' (not significant, P>0.05)';
    return
end
tag=' (not significant, P>0.1)';





function ci = bootper(bstat,alpha)
% percentile bootstrap CI 
pct1 = 100*alpha/2;
pct2 = 100-pct1;
lower = prctile(bstat,pct1); 
upper = prctile(bstat,pct2);
ci =[lower;upper];


%-------------------------------------------------------------------------    
function ci = bootnorm(bstat,alpha,stat)
% normal approximation interval
% A.C. Davison and D.V. Hinkley (1996), p198-200
 
se = std(bstat);   % standard deviation estimate
bias =mean(bstat-stat); % bias estimate
za = norminv(alpha/2);   % normal confidence point
lower = stat - bias + se*za; % lower bound
upper = stat - bias - se*za;  % upper bound
ci = [lower;upper];       


%-------------------------------------------------------------------------
function ci = bootcper(bstat,alpha,stat)
% corrected percentile bootstrap CI
% B. Efron (1982), "The jackknife, the bootstrap and other resampling
% plans", SIAM.
 
% stat is transformed to a normal random variable z0.
% z0 = invnormCDF[ECDF(stat)]
z_0 = norminv(sum(bstat<stat)/length(bstat));
z_alpha = norminv(alpha/2); % normal confidence point
 
% transform z0 back using the invECDF[normCDF(2z0-za)] and
% invECDF[normCDF(2z0+za)] 
pct1 = 100*normcdf(2*z_0-z_alpha); 
pct2 = 100*normcdf(2*z_0+z_alpha);
lower = prctile(bstat,pct2);  % inverse ECDF
upper = prctile(bstat,pct1);
ci = [lower;upper];
 