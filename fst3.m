function [f]=fst3(nv,pv)
%FST3 - Weir's F statistic (Fst) among more than 2 populations
%
% [f]=fst3(nv,pv)
% fst3([5 5 6 4],[0 .3 .167 .250]);

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if any(size(nv)~=size(pv))
    error('Not correct input')
end

  s=length(nv);    % num of subpoulations
  n=sum(nv);       % n1+n2;   Total individuals   
    

  nc = (1/(s-1))*(n-sum(nv.^2)./n);
  
% Weighted average frequency across subpopulations  
  p_hat=sum((nv./n).*pv);


%{
fprintf('\t');
fprintf(' pop%d\t',[1:s]);    
fprintf('All_W\tAll_UW\t');
fprintf('\n');
fprintf('N\t');
fprintf('%5d\t',nv);
fprintf('\n');
fprintf('p:\t');
fprintf('%.3f\t',pv);
fprintf('%.3f\t',p_hat);
fprintf('%.3f\t',mean(pv));
fprintf('\n');
%}

    MSP=(1/(s-1))*sum(nv.*(pv-p_hat).^2);
    MSG=(1./sum(nv-1)).*sum(nv.*pv.*(1-pv));
    if MSG==0
        f=0;
    else
        f = (MSP-MSG)./(MSP+(nc-1).*MSG);
    end
   
% Wright
% http://en.wikipedia.org/wiki/Fixation_index
% http://intl.genome.org/cgi/reprint/15/11/1496

    x=(nv.*(nv-1)/2);
    a=sum(x.*2.*(nv./(nv-1)).*pv.*(1-pv))./sum(x);    
    b=sum(2.*(n./(n-1)).*p_hat.*(1-p_hat));
    f2=1-a./b;
   
%Hua Chen
%Joint allele-frequency spectrum in closely related species
%

%    p_hat=nanmean(pv);
%    a=sum((pv-p_hat).^2)./2;    
%    b=p_hat.*(1-p_hat);
%    f3=a./b;
   
%    fprintf('Fst:\t%.3f (Weir), %.3f (Wright), %.3f (Chen)\n',f,f2,f3);


    return;    
    


%{
	    Locus: loc_2
	     N	     5	     5	     6	     4
	p:   2	 1.000	 0.700	 0.833	 0.750	 0.825	 0.821
	p:   4	 0.000	 0.300	 0.167	 0.250	 0.175	 0.179
	 For locus : loc_2
	Allele	  Capf	 Theta	Smallf	 Relat	Relatc	 Sig_a	 Sig_b	 Sig_w
	    2	 0.495	-0.045	 0.517	-0.060
	    4	 0.495	-0.045	 0.517	-0.060
	  All	 0.495	-0.045	 0.517	-0.060	-2.139	-0.013	 0.160	 0.150


    fst3([5 5 6 4],[0 .3 .167 .250]);

      

for k=1:nmark
    pv=pv1(:,k)';    

    X=(nv.*(nv-1)/2);
    p_hat=sum((nv./n).*pv);

    a=sum(X.*sum(2.*(nv./(nv-1)).*pv.*(1-pv))./X);
    b=sum(2.*(s./(s-1)).*pv.*(1-pv));
    f(k)=1-a./b;
end

    %}



