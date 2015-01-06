%Write point particle file

%clear all
%close all
clc

global caseFolder Home NGA src examples matlab

Home='/data/shichu/'; NGA=[Home,'NGA/'];src=[NGA,'src/']; examples=[NGA,'examples/'];
matlab=[Home,'Documents/MATLAB/'];
addpath(genpath(matlab));

%%CaseFolder path
%caseFolder=[Home,'Points'];
%caseFolder=[Home,'Points/addMass_stress_drag'];
%caseFolder=[Home,'Points/64Turb'];
%caseName='taylor';caseFolder=[Home,'Points/',caseName];
%caseFolder=[Home,'Points/drag'];
%caseFolder=[Home,'Points/oscillat'];

%caseName='oscillat';dirc=3;caseFolder=[Home,'Points/',caseName,'/',num2str(dirc)];
%caseFolder='/data/shichu/testPoints/scalarDiff/2D/perDiffusion';
%caseName='diffusion';
%dirc=0;

caseName='taylor';
caseFolder='/data/shichu/testPoints/Taylor';


turbConfig=[caseFolder,'/input/','turb.config'];

%%Flow field info

%Device info
devRange=[1 1];

%INITIAL CONDITIONS
if(strcmp(caseName,'taylor')==1)
init=7;
elseif(strcmp(caseName,'oscillat')==1)
init=4;
elseif(strcmp(caseName,'diffusion')==1)
switch(dirc)
case 0
init=8;
case 1
init=5;
case 2
init=5;
case 3
init=5;
end

else
init=1;
end

initial{1}='QUIESCENT';
initial{2}='SHEAR';
initial{3}='POISEUILLE';
initial{4}='OSCILLATORY';
initial{5}='BULK';
initial{6}='TURB';
initial{7}='TAYLOR';
initial{8}='QUIESCENT';

%Grid info
if(init==6||init==7)
X=[-pi pi];Nx=64;
Y=X;Ny=Nx;
Z=X;Nz=Nx;
elseif(init==8)
X=[-pi pi];Nx=128;
Y=X;Ny=Nx;
Z=X;Nz=Nx;

%X=[-2*pi 2*pi];Nx=128;
%Y=X;Ny=Nx;
%Z=X;Nz=Nx;
elseif(init==5)
X=[-pi pi];Nx=128;
Y=X;Ny=Nx;
Z=X;Nz=Nx;
else
X=[-2 2];Nx=64;
Y=X;Ny=Nx;
Z=X;Nz=Nx;
end

%PHYSICAL PARAMETERS
rho_f=1;
%nu=0.0571;
if(init==6)
nu=0.2; %nu should be greater than 0.2; for Nx=64, bluebottle can do at most Re_lambda=17;
else
nu=1;
end


%SIMULATION PARAMETERS
if(init==7)
duration=4;
elseif(init==8)
duration=1;
elseif(strcmp(caseName,'diffusion')==1)
duration=1;
else
duration=10;
end
CFL=0.5;


%SIMULATION DRIVING CONDITIONS
if(init==6)
turbA=1;
else
turbA=0;
end

if(strcmp(caseName,'oscillat')==1)
switch dirc
case 1
gradP=[-1 0 0];
case 2
gradP=[0 -1 0];
case 3
gradP=[0 0 -1];
end

else
gradP=[0 0 0];
end

g=[0 0 0];
%C_add=0;
%C_stress=0;

C_add=0.5;
C_stress=1;
C_drag=1;
%C_drag=0;

%BOUNDARY CONDITIONS
bc{1}='NEUMANN';
bc{2}='DIRICHLET 1.0 0.0';
bc{3}='PERIODIC';


p_bc=3;
u_bc=p_bc;v_bc=u_bc;w_bc=u_bc;
p_bc1=p_bc;p_bc2=p_bc;p_bc3=p_bc;
u_bc1=u_bc;u_bc2=u_bc;u_bc3=u_bc;
v_bc1=u_bc1;v_bc2=u_bc2;v_bc3=u_bc3;
w_bc1=u_bc1;w_bc2=u_bc2;w_bc3=u_bc3;

if(init==4)
u_bc=3;
switch dirc
case 1
u_bc1=u_bc;u_bc2=u_bc;u_bc3=2;
case 2
u_bc1=u_bc;u_bc2=u_bc;u_bc3=2;
case 3
u_bc1=2;u_bc2=u_bc;u_bc3=u_bc;
end

elseif(init==6||init==7)
u_bc=3;
u_bc1=u_bc;u_bc2=u_bc;u_bc3=u_bc;
elseif(init==8)
u_bc=1;
u_bc1=u_bc;u_bc2=u_bc;u_bc3=u_bc;
%else
%error('wrong initial condition');
end



if(strcmp(caseName,'diffusion')==1)
bc{2}='DIRICHLET 10 0.0';

p_bc=1;
p_bc1=p_bc;p_bc2=p_bc;p_bc3=p_bc;
switch(dirc)
case 0
u_bc=p_bc;v_bc=u_bc;w_bc=u_bc;
u_bc1=u_bc;u_bc2=u_bc;u_bc3=u_bc;
v_bc1=u_bc1;v_bc2=u_bc2;v_bc3=u_bc3;
w_bc1=u_bc1;w_bc2=u_bc2;w_bc3=u_bc3;
case 1
u_bc=2; v_bc=1;w_bc=1;
u_bc1=2;u_bc2=1;u_bc3=1;
v_bc1=1;v_bc2=u_bc2;v_bc3=u_bc3;
w_bc1=1;w_bc2=u_bc2;w_bc3=u_bc3;
case 2
u_bc=1; v_bc=2;w_bc=1;
u_bc1=1;u_bc2=1;u_bc3=1;
v_bc1=1;v_bc2=2;v_bc3=u_bc3;
w_bc1=1;w_bc2=u_bc2;w_bc3=u_bc3;
case 3
u_bc=1; v_bc=1;w_bc=1;
u_bc1=1;u_bc2=1;u_bc3=1;
v_bc1=1;v_bc2=u_bc2;v_bc3=u_bc3;
w_bc1=1;w_bc2=u_bc2;w_bc3=2;
end
end


%SOLVABILITY ENFORCEMENT PLANE
%For inflow-outflow problems, set out_plane to be the outflow plane.
%    For periodic problems with a mean flow direction, set out_plane to be one of the periodic directions not in the mean flow direction.
%    For periodic problems with no mean flow direction, set out_plane as HOMOGENEOUS, which will spread the modification equally over all sides of the domain.
solvability{1}='EAST';
solvability{2}='NORTH';
solvability{3}='TOP';
solvability{4}='HOMOGENEOUS';
if(init==4)
switch dirc
case 1
out_plane=3;
case 2
out_plane=3;
case 3
out_plane=1;
end

else
out_plane=4;
end

%%Write turb config
%%==================================================%%
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
fprintf(fid,'INITIAL CONDITION\n');
fprintf(fid,'init_cond %s\n',initial{init});

fprintf(fid,'\n');
fprintf(fid,'SOLVABILITY ENFORCEMENT PLANE\n');
%fprintf(fid,'out_plane TOP\n');
%fprintf(fid,'out_plane EAST\n');
fprintf(fid,'out_plane %s\n',solvability{out_plane});


fprintf(fid,'\n');
fprintf(fid,'SIMULATION DRIVING CONDITIONS\n');
fprintf(fid,'turbA %f\n',turbA);
fprintf(fid,'gradP.x %f 0.\n',gradP(1));
fprintf(fid,'gradP.y %f 0.\n',gradP(2));
fprintf(fid,'gradP.z %f 0.\n',gradP(3));
fprintf(fid,'g.x %f 0.\n',g(1));
fprintf(fid,'g.y %f 0.\n',g(2));
fprintf(fid,'g.z %f 0.\n',g(3));

fprintf(fid,'\n');
fprintf(fid,'Add mass %f\n',C_add);
fprintf(fid,'Fluid stress %f\n',C_stress);
fprintf(fid,'Drag force %f\n',C_drag);

fclose(fid);
type(turbConfig)


