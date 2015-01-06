%Write point particle file
clear all
close all
clc

global caseFolder Home NGA src examples matlab

Home='/data/shichu/'; NGA=[Home,'NGA/'];src=[NGA,'src/']; examples=[NGA,'examples/'];
matlab=[Home,'Documents/MATLAB/'];
addpath(genpath(matlab));

caseName='diffusion';
%caseFolder=[Home,'testPoints/',caseName];

dirc=0;
switch(dirc)
case 0
folder='diffusion';
case 1
folder='xconvDiff';
case 2
folder='yconvDiff';
case 3
folder='zconvDiff';
end

caseFolder=['/data/shichu/testPoints/',folder];

run scalarWriteBlue.m;
run turbWriteBlue.m;
run recordWriteBlue.m;

