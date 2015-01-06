clear all
close all
clc
%addpath(genpath('/home/asiera/bluebottle/src/tools/matlab'));
%rmpath(genpath('/home/asiera/bluebottle/src/tools/matlab'));
addpath(genpath('/home/shichu/Documents/MATLAB/blueRead'));

folder='/home/asiera/bluebottle/cases/concentration_gradient';


time=0.01;
bluePicPath='/home/shichu/Documents/MATLAB/bluePic';
X=10;Y=10;Z=10; %half length of 3 dir
H=X;
%V=X*Y*Z*8;
V=2*H;
M=8;%the maximum number for fourier expansion
reali=8;
nq0=cell(reali,1);
sin_coeff=nq0;
cos_coeff=nq0;
N=500;

for kk=1:reali

z=2*X*rand(N,1)-X;
x=z;

z=asin(z/H)*2/pi*H;


density=ones(N,1);
[nq0{kk} sin_coeff{kk} cos_coeff{kk}]=num_coeff(density,z,M,time,bluePicPath,'numDen');
end


%fourier_coeff(density,x,M,time,bluePicPath,'numDen');
%fourier_coeff(density,y,M,time,bluePicPath,'numDen');
fourier_coeff(density,z,M,time,bluePicPath,'numDen');


for kk=1:reali
clear buf1 buf2
c_0(kk)=nq0{kk};

buf2=sin_coeff{kk};
buf1=cos_coeff{kk};

c_1(kk)=buf1(1);
c_2(kk)=buf1(2);
c_3(kk)=buf1(3);

s_1(kk)=buf2(1);
s_2(kk)=buf2(2);
s_3(kk)=buf2(3);
end


for kk=1:reali
C_0(kk)=sum(c_0(1:kk))/kk;

C_1(kk)=sum(c_1(1:kk))/kk;
C_2(kk)=sum(c_2(1:kk))/kk;
C_3(kk)=sum(c_3(1:kk))/kk;

S_1(kk)=sum(s_1(1:kk))/kk;
S_2(kk)=sum(s_2(1:kk))/kk;
S_3(kk)=sum(s_3(1:kk))/kk;
end


COS_0=0;
COS_COEFF=zeros(M,1);
SIN_COEFF=COS_COEFF;


for jj=1:M
for kk=1:reali
clear buf1 buf2
buf2=sin_coeff{kk};
buf1=cos_coeff{kk};
COS_0=COS_0+nq0{kk};

COS_COEFF(jj)=COS_COEFF(jj)+buf1(jj);
SIN_COEFF(jj)=SIN_COEFF(jj)+buf2(jj);
end
COS_0=COS_0/reali;
COS_COEFF(jj)=COS_COEFF(jj)/reali;
SIN_COEFF(jj)=SIN_COEFF(jj)/reali;
end


K_L=(1:M)'*pi/H;
cesaro_coeff=1-(1:M)'/(M+1);
NUM=100;%divide Y into NUM parts, to see the vel distribution
yy=linspace(min(x),max(x),NUM);

for j=1:NUM
%  NQ(j)=nq_0+sum(cos_coeff.*cos(K_L*yy(j)))-i*sum(sin_coeff.*sin(K_L*yy(j)));
  NQ(j)=COS_0+sum(cesaro_coeff.*COS_COEFF.*cos(K_L*yy(j)))-i*sum(cesaro_coeff.*SIN_COEFF.*sin(K_L*yy(j)));
end
handle=figure;
plot(yy,NQ*V/N,'-');
xlabel('x');
ylabel('num density');
title(['scalar distribution,t= ',num2str(time)]);
filename=[bluePicPath,'/','average_distri','.eps'];
saveas(handle,filename);



handle=figure;
plot((C_0),'r-*');
xlabel('num of realizations');
ylabel('averaged coeff');
legend('0-th order cos coeff');

filename=[bluePicPath,'/','0_cos','.eps'];
saveas(handle,filename);

handle=figure;
plot((C_1),'r-*');
xlabel('num of realizations');
ylabel('averaged coeff');
legend('1-st order cos coeff');
filename=[bluePicPath,'/','1_cos','.eps'];
saveas(handle,filename);

handle=figure;
plot((C_2),'r-*');
xlabel('num of realizations');
ylabel('averaged coeff');
legend('2-nd order cos coeff');
filename=[bluePicPath,'/','2_cos','.eps'];
saveas(handle,filename);



handle=figure;
plot((S_1)*i,'r-*');
xlabel('num of realizations');
ylabel('averaged coeff');
legend('1-st order sin coeff');
filename=[bluePicPath,'/','1_sin','.eps'];
saveas(handle,filename);

handle=figure;
plot((S_2)*i,'r-*');
xlabel('num of realizations');
ylabel('averaged coeff');
legend('2-nd order sin coeff');
filename=[bluePicPath,'/','2_sin','.eps'];
saveas(handle,filename);







