kimstephan02(haplodata,hmarkinfo)

% Kim, Y. and Stephan, W. 2002, Detecting a local signature of genetic
% hitchhiking along a recombining chromosome. Genetics 160:765-777.



% s, x
% Rn - recombination rate
% n - sample size

afactor=0.5;
resol=0.001;

f_al=sum(1./[1:(n-1)]);

for k=1:seqL
  r=abs(k-x*seqL)*Rn;
  cst=1-power(afactor./s./f_al,0.5*r./s./f_al);
  ncee=cst./resol+0.5;
  
  
end