function [sew,fww,mhz,ni_mh] = sewfww(data)
%SEWFWW - Calculate ratio of adaptive nonsyn. sub., ie., Smith and Eyre-Walker's and Fay, Wyckoff and Wu's alpha_bar and Mantel-Haenszel's z
% Calculate Smith and Eyre-Walker's and Fay, Wyckoff and Wu's alpha_bar and Mantel-Haenszel's z */
%  Syntax: [sew,fww,mhz] = sewfww(data)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

dDn=data(:,1);
dPn=data(:,3);
dDs=data(:,5);
dPs=data(:,7);

%dn pn ds ps
LdDn=data(:,2);
LdPn=data(:,4);
LdDi=data(:,6);
LdPi=data(:,8);

smeasured=sum(data(:,[1,3,5,7]),1);
ssites=sum(data(:,[2,4,6,8]),1);

Dn_bar=smeasured(1);
Pn_bar=smeasured(2);
Ds_bar=smeasured(3);
Ps_bar=smeasured(4);

LDn_bar=ssites(1);
LPn_bar=ssites(2);
LDs_bar=ssites(3);
LPs_bar=ssites(4);
z_bar=sum((dPn.*LdPi)./(dPs+1));

% Adam's method 2
sew=1.0-Ds_bar*LDn_bar/Dn_bar/LDs_bar*z_bar/LPn_bar;
%sew = 1.0 - ((Ds_bar/LDs_bar) / (Dn_bar/LDn_bar)) * z_bar / LPn_bar;
%sew = 1.0 - Ds_bar / Dn_bar / * sum(dPn./(dPs+1));

if (nargout>1),
% Fay, Wyckoff and Wu
fww=1.0-Ds_bar*LDn_bar/Dn_bar/LDs_bar*Pn_bar*LPs_bar/Ps_bar/LPn_bar;
%fww = 1.0 - ((Ds_bar/LDs_bar) / (Dn_bar/LDn_bar)) * (( Pn_bar/LPn_bar) * (LPs_bar/Ps_bar));
% fww = 1.0 - Ds_bar / Dn_bar * Pn_bar / Ps_bar;
end

if (nargout>2),
	% For Mantel-Haenszel z statistic */
	r2=dDn+dPn;
	c2=dDn+dDs;
	N=c2+dPn+dPs;

	r2(N==1)=[]; c2(N==1)=[]; N(N==1)=[];
	dDn(N==1)=[]; dDs(N==1)=[]; dPs(N==1)=[]; dPn(N==1)=[];

	sum_d_m_delta=sum(dDn - (r2.*c2)./N);
	sum_V=sum(c2.*r2.*(dPs + dPn).*(dPs + dDs)./N./N./(N - 1.0));

	% For Mantel-Haenszel z statistic
	mhz = (abs( sum_d_m_delta ) - 0.5) / sqrt( sum_V );
	%fprintf('\nMHz: %0.6g\n', mhz );

end


if (nargout>3),
    dn=dDn./LdDn;
    ds=dDs./LdDi;
    pn=dPn./LdPn;
    ps=dPs./LdPi;

    %[ni_mh]=nimh(dn,ds,pn,ps);
    [ni_mh]=nimh(dDn,dDs,dPn,dPs);
end