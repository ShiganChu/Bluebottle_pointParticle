clear all
close all
clc
%addpath(genpath('/home/asiera/bluebottle/src/tools/matlab'));
rmpath(genpath('/home/asiera/bluebottle/src/tools/matlab'));
addpath(genpath('/home/shichu/Documents/MATLAB/blueRead'));
casename='/home/asiera/bluebottle/cases/channel_1000';
time=4.32;
bluePicPath='/home/shichu/Documents/MATLAB/bluePic';
for tt=1:4
time=0.01+(tt-1)*0.5;
%time=0.01;
[u, v, w] = cgns_read_part_vel(casename, time);
[x, y, z] = cgns_read_part_position(casename, time);
%[Fx, Fy, Fz] = cgns_read_part_force_hydro(casename, time);
[Fx, Fy, Fz] = cgns_read_part_force_total(casename, time);
[Ox, Oy, Oz] = cgns_read_part_omega(casename, time);

X=30;Y=15;Z=7.5; %half length of 3 dir
H=Y;
%V=X*Y*Z*8;
V=2*H;
M=4;%the maximum number for fourier expansion



N=length(u);%the number of particles
density=ones(N,1);



fourier_coeff(density,y,M,time,bluePicPath,'numDen');
%fourier_coeff(u,y,M,time,bluePicPath,'x-vel');%u is the scalar,y is the position in y, M is the expand order
%fourier_coeff(v,y,M,time,bluePicPath,'y-vel');
%fourier_coeff(Oz,y,M,time,bluePicPath,'z-angVel');
end
