%Write point particle file

%clear all
%close all
clc

global caseFolder Home NGA src examples matlab

Home='/data/shichu/'; NGA=[Home,'NGA/'];src=[NGA,'src/']; examples=[NGA,'examples/'];
matlab=[Home,'Documents/MATLAB/'];
addpath(genpath(matlab));

%%CaseFolder path
%caseName='diffusion';
%caseFolder=[Home,'testPoints/',caseName];


scalarConfig=[caseFolder,'/input/','scalar.config'];

%%Flow field info


%Device info
devRange=[1 1];

%sc_initial CONDITIONS
if(strcmp(caseName,'taylor')==1)
init=7;
elseif(strcmp(caseName,'oscillat')==1)
init=4;
elseif(strcmp(caseName,'diffusion')==1)
init=8;
else
init=1;
end

sc_initial{1}='QUIESCENT';
sc_initial{2}='SHEAR';
sc_initial{3}='POISEUILLE';
sc_initial{4}='OSCILLATORY';
sc_initial{5}='BULK';
sc_initial{6}='scalar';
sc_initial{7}='TAYLOR';
sc_initial{8}='DIFFUSION';

%Grid info
if(init==8)
X=[-pi pi];Nx=128;
Y=X;Ny=Nx;
Z=X;Nz=Nx;
end

%PHYSICAL PARAMETERS
if(init==8)
DIFF=1;
end




%BOUNDARY CONDITIONS
bc{1}='NEUMANN';
bc{2}='DIRICHLET 0.0 0.0';
bc{3}='PERIODIC';
sc_bc=3;
%sc_bc=3;
sc_bc1=sc_bc;sc_bc2=sc_bc;sc_bc3=sc_bc;


%%Write scalar config
%%==================================================%%
[fid,msg] = fopen(scalarConfig, 'w'); % read only 
fprintf(fid,'n 1\n');
fprintf(fid,'\n');
fprintf(fid,'PHYSICAL PARAMETERS\n');
fprintf(fid,'DIFF %f\n',DIFF);

fprintf(fid,'\n');
fprintf(fid,'BOUNDARY CONDITIONS\n');
fprintf(fid,'SCALAR\n');
fprintf(fid,'sc_bc.scW %s\n',bc{sc_bc1});
fprintf(fid,'sc_bc.scE %s\n',bc{sc_bc1});
fprintf(fid,'sc_bc.scS %s\n',bc{sc_bc2});
fprintf(fid,'sc_bc.scN %s\n',bc{sc_bc2});
fprintf(fid,'sc_bc.scB %s\n',bc{sc_bc3});
fprintf(fid,'sc_bc.scT %s\n',bc{sc_bc3});

fprintf(fid,'\n');
fprintf(fid,'INITIAL CONDITION\n');
fprintf(fid,'sc_init_cond %s\n',sc_initial{init});

fclose(fid);
type(scalarConfig)


