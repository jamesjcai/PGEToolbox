function [ci] = binoci(n, p, alpha, type)
%Binomial Confidence Intervals
%N=Total number of trials
%P=prop, population proportion of successes

if size(n, 2) < size(n, 1)
    error('Needs a row vector')
end

if nargin < 4, type = 2; end
if nargin < 3, alpha = 0.05; end
% * standard (Wald) interval *
za = norminv(alpha/2);


switch (type)
    case 1
        if n * min(p, 1-p) < 5,
            warning(sprintf(['The actual coverage probability of the standard\n', ...
                'interval is poor for p near 0 or 1']));
            end

            se = sqrt(p*(1 - p)/n);
            %Normal Approximation using z-distribution
            % normal confidence point  1.96,
            %100(1-alpha/2)th percentile of the standard normal distribution

            %Normal Approximation using t-distribution
            %za = tinv(alpha/2,n-1);
            lower = p + se * za; % lower bound
            upper = p - se * za; % upper bound
            ci = [lower; upper];

        case 2
            % Recommendthe Wilson or the equal-tailed Jeffreys prior interval for small nn ? 40).
            % For larger n, the Wilson, the Jeffreys andthe Agresti–Coull intervals are all comparable,
            % SEE: http://www-stat.wharton.upenn.edu/~tcai/paper/Binomial-StatSci.pdf

            % Wilson CI
            %X=n*p; a=(X+za*za/2)/(n+za*za); b=za*sqrt(n)*sqrt(p*(1-p)+za*za/(4*n))/(n+za*za);
            %ci=[a+b;a-b];
            %or
            temp = za .* sqrt(p.*(1 - p)./n+.25.*(za ./ n).^2);
            tu = (p + .5 .* za .* za ./ n - temp) ./ (1 + za .* za ./ n);
            tl = (p + .5 .* za .* za ./ n + temp) ./ (1 + za .* za ./ n);
            ci = [tl; tu];
        case 3
            % Agresti–Coull intervals
            X = n * p + za * za / 2;
            n = n + za * za;
            p = X / n;
            tu = p - za * sqrt(p*(1 - p)/n);
            tl = p + za * sqrt(p*(1 - p)/n);
            ci = [tl; tu];
        end


        function [ss] = i_BinP(N, p, x1, x2)
            q = p / (1 - p);
            s = 0;
            k = 0;
            tot = 0;
            v = 1;
            while (k <= N)
                tot = tot + v;
                if (k >= x1 & k <= x2), s = s + v; end
                if (tot > 1e30), s = s / 1e30;
                    tot = tot / 1e30;
                    v = v / 1e30;
                end
                k = k + 1;
                v = v * q * (N + 1 - k) / k;
            end
            ss = s / tot;


            % Error bars denote 95% confidence intervals for binomial expectations.


            %var vTL=2.5; var vTU=2.5; var vCL=95
            %function CalcCL(form) {
            %    vTL = eval(form.TL.value)
            %    vTU = eval(form.TU.value)
            %    vCL = 100-(vTL+vTU)
            %    form.CL.value = ''+vCL
            %    }

            %function CalcTails(form) {
            %    vCL = eval(form.CL.value)
            %    vTU = (100-vCL)/2
            %    vTL = vTU
            %    form.TL.value = ''+vTL
            %    form.TU.value = ''+vTU
            %    }

            %function CalcDummy(form) { }

            %function CalcBin(form) {
            %    var vx = eval(form.x.value)
            %    var vN = eval(form.N.value)
            %    var vP = vx/vN
            %    form.P.value = Fmt(vP)
            %    if(vx==0)
            %	{ form.DL.value = "0.0000" } else
            %        { var v=vP/2; vsL=0; vsH=vP; var p=vTL/100
            %        while((vsH-vsL)>1e-5) { if(BinP(vN,v,vx,vN)>p) { vsH=v; v=(vsL+v)/2 } else { vsL=v; v=(v+vsH)/2 } }
            %        form.DL.value = Fmt(v) }
            %    if(vx==vN)
            %        { form.DU.value = "1.0000" } else
            %        { var v=(1+vP)/2; vsL=vP; vsH=1; var p=vTU/100
            %        while((vsH-vsL)>1e-5) { if(BinP(vN,v,0,vx)<p) { vsH=v; v=(vsL+v)/2 } else { vsL=v; v=(v+vsH)/2 } }
            %        form.DU.value = Fmt(v) }
            %    }

            %function BinP(N,p,x1,x2) {
            %    var q=p/(1-p); var k=0; var v = 1; var s=0; var tot=0
            %    while(k<=N) {
            %        tot=tot+v
            %        if(k>=x1 & k<=x2) { s=s+v }
            %        if(tot>1e30){s=s/1e30; tot=tot/1e30; v=v/1e30}
            %        k=k+1; v=v*q*(N+1-k)/k
            %        }
            %    return s/tot
            %    }

            %function CalcPois(form) {
            %    var vZ = eval(form.Z.value)
            %    if(vZ==0)
            %        { form.QL.value = "0.0000" } else
            %        { var v=0.5; var dv=0.5; var p=vTL/100
            %        while(dv>1e-7) { dv=dv/2; if(PoisP((1+vZ)*v/(1-v),vZ,1e10)>p) { v=v-dv } else { v=v+dv } }
            %        form.QL.value = Fmt((1+vZ)*v/(1-v)) }
            %    if(vTU==0)
            %        { form.QU.value = "Infinity" } else
            %        { var v=0.5; var dv=0.5; var p=vTU/100
            %        while(dv>1e-7) { dv=dv/2; if(PoisP((1+vZ)*v/(1-v),0,vZ)<p) { v=v-dv } else { v=v+dv } }
            %        form.QU.value = Fmt((1+vZ)*v/(1-v)) }
            %    }

            %function PoisP(Z,x1,x2) {
            %    var q=1; var tot=0; var s=0; var k=0
            %    while(k<Z || q>(tot*1e-10)) {
            %        tot=tot+q
            %        if(k>=x1 & k<=x2) { s=s+q }
            %        if(tot>1e30){s=s/1e30; tot=tot/1e30; q=q/1e30}
            %        k=k+1; q=q*Z/k
            %        }
            %    return s/tot
            %    }

            %function Fmt(x) {
            %var v
            %if(x>=0) { v=''+(x+0.00005) } else { v=''+(x-0.00005) }
            %return v.substring(0,v.indexOf('.')+5)
            %}
