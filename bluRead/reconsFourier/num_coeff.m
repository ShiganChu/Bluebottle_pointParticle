function [nq_0 sin_coeff cos_coeff]=num_coeff(v,y,M,time,path,name)
%Get the fourier coefficient of v(y)

H=(max(y)-min(y))/2;
%V=X*Y*Z*8;
V=2*H;

N=length(v);%the number of particles

nq=zeros(M,1);
nq_minus=nq;
K_L=(1:M)'*pi/H;

for l=1:M
  k_l=l*pi/H;
  nq(l)       =sum(1/V*exp(i*k_l*y).*v);
  nq_minus(l) =sum(1/V*exp(-i*k_l*y).*v);
end
sin_coeff=nq-nq_minus;
cos_coeff=nq+nq_minus;
nq_0=mean(v)*N/V;


% the y must be fixed!!
%nq=nq_0+sum((nq(l)+nq_minus(l))*cos(k_l*y))+i*sum((nq(l)-nq_minus(l))*sin(k_l*y))
NUM=100;%divide Y into NUM parts, to see the vel distribution
yy=linspace(min(y),max(y),NUM);


cesaro_coeff=1-(1:M)'/(M+1);

for j=1:NUM
%  NQ(j)=nq_0+sum(cos_coeff.*cos(K_L*yy(j)))-i*sum(sin_coeff.*sin(K_L*yy(j)));
  NQ(j)=nq_0+sum(cesaro_coeff.*cos_coeff.*cos(K_L*yy(j)))-i*sum(cesaro_coeff.*sin_coeff.*sin(K_L*yy(j)));
end


