function frdevo(p0,N,s1,s2)
%
%
%AA - 1
%Aa - s*h
%aa - s
%
% s=s2; h=s1/s2;
%
% s2=s; s1=s*h;



N=10000;
p0=0.5;
% s1=0.9; s2=0.6;
s1=-0.1; s2=-0.6;

n0=2*N*p0;

p1v=zeros(1,10);
for t=1:10
    p1=i_p1(p0,s1,s2);
    p1v(t)=p1; 
    p0=p1;
end

plot(p1v,'-o')

n=10; N=50;
M=zeros(10);
for i=1:n-1
for j=i+1:n
    pa=i/2*N;   
    pb=i_p1(pa,s1,s2);
    M(i,j)=binopdf(j,2*N,pb);
end
end
    

M
expm(M*100)

function p1=i_p1(p0,s1,s2)

a=1+s2.*p0+s1.*(1-p0);
b=1+s2.*p0.^2+s1.*p0.*(1-p0);
p1=p0.*a./b;

