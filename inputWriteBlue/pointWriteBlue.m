%Write point particle file

clear all
close all
clc

global caseFolder Home NGA src examples matlab

Home='/data/shichu/'; NGA=[Home,'NGA/'];src=[NGA,'src/']; examples=[NGA,'examples/'];
matlab=[Home,'Documents/MATLAB/'];
addpath(genpath(matlab));

%%CaseFolder path
%caseFolder=[Home,'Points'];

caseName='taylor';
caseFolder=[caseFolder,'/',caseName];

caseFolder='/data/shichu/testPoints/oilSourceScalar';
%caseFolder='/data/shichu/testPoints/testOilSourceScalar';

%caseFolder='/data/shichu/testPoints/sourceScalar';
pointConfig=[caseFolder,'/input/','point.config'];
flowConfig=[caseFolder,'/input/','flow.config'];
turbConfig=[caseFolder,'/input/','turb.config'];

%%particle info, contains how many type of point particles
numOfClass=[100];
radiusClass=[0.15];
densityClass=[0.87];

%radiusClass=[0.5];
%densityClass=[1];

classNum=length(numOfClass);

%%Read turb config
[fid,msg] = fopen(turbConfig, 'r'); % read only 
buf=fgets(fid);
buf1=fscanf(fid,'(Xs, Xe, Xn) %f %f %d');
fscanf(fid,'\n');
buf2=fscanf(fid,'(Ys, Ye, Yn) %f %f %d');
fscanf(fid,'\n');
buf3=fscanf(fid,'(Zs, Ze, Zn) %f %f %d');
X=buf1(1:2);Nx=buf1(3);
Y=buf2(1:2);Ny=buf2(3);
Z=buf3(1:2);Nz=buf3(3);
fclose(fid);

%write the particle config file in bluebottle
[fid,errmsg] = fopen(pointConfig, 'w') % write only 
fprintf(fid,'n %d\n',sum(numOfClass));
for ptype=1:length(numOfClass)
   for k=1:numOfClass(ptype)
   fprintf(fid,'r %f\n',radiusClass(ptype));

xbuf=(X(2)-X(1))*rand()+X(1);
ybuf=(Y(2)-Y(1))*rand()+Y(1);
zbuf=(Z(2)-Z(1))*rand()+Z(1);

zbuf=0;
   fprintf(fid,'(x, y, z) %f %f %f\n',xbuf,ybuf,zbuf);
   fprintf(fid,'rho %f\n',densityClass(ptype));
   fprintf(fid,'rotating %d\n',1);
   fprintf(fid,'\n');
   end
end
fclose(fid);
 
caseFolder
