%function y=chencohen

Ne=10^4;
u=2.5*10^-8;
p=1./(2*Ne);
t=25*10^4;
k=1.16*10^-3;
M=k/(4*Ne*u*t+2*Ne);


%s=-10^-5;

gam=linspace(-10,0,100);   %-2.72;
%gam=[-2.72 -2.1];

%s=gam./Ne;


h=.5;

y=zeros(1,100);
for k=1:100
  %y(k)=abs(i_fixprob(p,gam(k),h)-M);
  y(k)=abs(fixprobh(p,h,gam(k))-M);
end

%semilogy(gam,y,'-o');
semilogy(gam,y,'-')
hold on
hline(0)

