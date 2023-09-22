function [rlen] = varwinradii(d, howmanyd)
%VARWINRADII - creates sliding widowns with varible sizes.
% It returns borders V for sliding windows with same number of D.

if nargin < 2
    howmanyd = 5;
end

%if islogical(d)
%    d=find(d);
%else
%dori=d;

d = d - min(d) + 1;
L = max(d(:));
V = false(1, L);
V(d) = true;
%howmanyd=round(length(d)/10)


options = optimset('fminbnd');


%    rlen=zeros(1,L);
%    for pos=1:L
%        [x,fval] = fminbnd(@b20,5,L,options,pos,howmanyd,V,L);
%        x=round(x);
%        rlen(pos)=x;
%    end

n = length(d);
rlen = zeros(1, n);
pause
for k = 1:n
    pos = d(k);

    [x, fval] = fminbnd(@b20, 5, L, options, pos, howmanyd, V, L);
    x = round(x);
    rlen(k) = x;

end


    function [res] = b20(x, pos, scalen, V, L)
        a1 = round(max(1, pos-x));
        a2 = round(min(L, pos+x));
        res = abs(sum(V(a1:a2))-scalen-1);