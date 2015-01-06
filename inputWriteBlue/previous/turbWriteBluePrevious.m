%Write point particle file

clear all
close all
clc

global caseFolder Home NGA src examples matlab

Home='/data/shichu/'; NGA=[Home,'NGA/'];src=[NGA,'src/']; examples=[NGA,'examples/'];
matlab=[Home,'Documents/MATLAB/'];
addpath(genpath(matlab));

%%CaseFolder path
caseFolder=[Home,'Points'];
%caseFolder=[Home,'Points/64Turb'];
turbConfig=[caseFolder,'/input/','turb.config'];

%%Flow field info

%Grid info
X=[-pi pi];Nx=64;
Y=X;Ny=Nx;
Z=X;Nz=Nx;

%Device info
devRange=[1 1];

%PHYSICAL PARAMETERS
rho_f=1;
%nu=0.0571;

nu=0.2; %nu should be greater than 0.2; for Nx=64, bluebottle can do at most Re_lambda=17;

%SIMULATION PARAMETERS
duration=10;
CFL=0.5;
turbA=1;

%BOUNDARY CONDITIONS
bc{1}='NEUMANN';
bc{2}='DIRICHLET';
bc{3}='PERIODIC';
p_bc=3;
u_bc=p_bc;v_bc=u_bc;w_bc=u_bc;
p_bc1=p_bc;p_bc2=p_bc;p_bc3=p_bc;
u_bc1=u_bc;u_bc2=u_bc;u_bc3=u_bc;
v_bc1=v_bc;v_bc2=v_bc;v_bc3=v_bc;
w_bc1=w_bc;w_bc2=w_bc;w_bc3=w_bc;


%%Write turb config
[fid,msg] = fopen(turbConfig, 'w'); % read only 
fprintf(fid,'DOMAIN\n');
fprintf(fid,'(Xs, Xe, Xn) %f %f %d\n',X(1),X(2),Nx);
fprintf(fid,'(Ys, Ye, Yn) %f %f %d\n',Y(1),Y(2),Ny);
fprintf(fid,'(Zs, Ze, Zn) %f %f %d\n',Z(1),Z(2),Nz);

fprintf(fid,'\n');
fprintf(fid,'GPU DOMAIN DECOMPOSITION\n');
fprintf(fid,'DEV RANGE %d %d\n',devRange(1),devRange(2));
fprintf(fid,'n 1\n');
fprintf(fid,'(Xs, Xe, Xn) %f %f %d\n',X(1),X(2),Nx);
fprintf(fid,'(Ys, Ye, Yn) %f %f %d\n',Y(1),Y(2),Ny);
fprintf(fid,'(Zs, Ze, Zn) %f %f %d\n',Z(1),Z(2),Nz);
fprintf(fid,'E -1 W -1 N -1 S -1 T -1 B -1\n');

fprintf(fid,'\n');
fprintf(fid,'PHYSICAL PARAMETERS\n');
fprintf(fid,'rho_f %f\n',rho_f);
fprintf(fid,'nu %f\n',nu);

fprintf(fid,'\n');
fprintf(fid,'SIMULATION PARAMETERS\n');
fprintf(fid,'duration %f\n',duration);
fprintf(fid,'CFL %f\n',CFL);
fprintf(fid,'pp_max_iter 1000\n');
fprintf(fid,'pp_residual 1e-6\n');

fprintf(fid,'\n');
fprintf(fid,'BOUNDARY CONDITIONS\n');
fprintf(fid,'PRESSURE\n');
fprintf(fid,'bc.pW %s\n',bc{p_bc1});
fprintf(fid,'bc.pE %s\n',bc{p_bc1});
fprintf(fid,'bc.pS %s\n',bc{p_bc2});
fprintf(fid,'bc.pN %s\n',bc{p_bc2});
fprintf(fid,'bc.pB %s\n',bc{p_bc3});
fprintf(fid,'bc.pT %s\n',bc{p_bc3});
fprintf(fid,'X-VELOCITY\n');
fprintf(fid,'bc.uW %s\n',bc{u_bc1});
fprintf(fid,'bc.uE %s\n',bc{u_bc1});
fprintf(fid,'bc.uS %s\n',bc{u_bc2});
fprintf(fid,'bc.uN %s\n',bc{u_bc2});
fprintf(fid,'bc.uB %s\n',bc{u_bc3});
fprintf(fid,'bc.uT %s\n',bc{u_bc3});
fprintf(fid,'Y-VELOCITY\n');
fprintf(fid,'bc.vW %s\n',bc{v_bc1});
fprintf(fid,'bc.vE %s\n',bc{v_bc1});
fprintf(fid,'bc.vS %s\n',bc{v_bc2});
fprintf(fid,'bc.vN %s\n',bc{v_bc2});
fprintf(fid,'bc.vB %s\n',bc{v_bc3});
fprintf(fid,'bc.vT %s\n',bc{v_bc3});
fprintf(fid,'Z-VELOCITY\n');
fprintf(fid,'bc.wW %s\n',bc{w_bc1});
fprintf(fid,'bc.wE %s\n',bc{w_bc1});
fprintf(fid,'bc.wS %s\n',bc{w_bc2});
fprintf(fid,'bc.wN %s\n',bc{w_bc2});
fprintf(fid,'bc.wB %s\n',bc{w_bc3});
fprintf(fid,'bc.wT %s\n',bc{w_bc3});

fprintf(fid,'\n');
fprintf(fid,'SOLVABILITY ENFORCEMENT PLANE\n');
fprintf(fid,'out_plane TOP\n');

fprintf(fid,'\n');
fprintf(fid,'SIMULATION DRIVING CONDITIONS\n');
fprintf(fid,'turbA %f\n',turbA);

fclose(fid);
 

