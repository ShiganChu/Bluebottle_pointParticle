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
%pointConfig=[caseFolder,'/input/','point.config'];
%flowConfig=[caseFolder,'/input/','flow.config'];
%turbConfig=[caseFolder,'/input/','turb.config'];
recordConfig=[caseFolder,'/input/','record.config'];

%%particle info, contains how many type of point particles
if(strcmp(caseName,'diffusion')==1)
flow_cgns_dt=0.005;
point_cgns_dt=-1;
scalar_dt=flow_cgns_dt;
else
flow_cgns_dt=0.01;
point_cgns_dt=0.01;
scalar_dt=0.01;
end
paraview_dt=-0.01;
precursor_dt=-0.01;
%%TODO Figure out the exact value for those
%- RESTART_STOP is a new boolean record.config parameter that specifies whether to stop the simulation after writing the restart file. The choices are: 
%restart_stop_dt= 0 = do not stop after writing output file  
%restart_stop_dt= 1 = stop after writing output file
restart_dt=14;
restart_stop_dt=1;

%scalar_dt=1e-3;


%write the particle config file in bluebottle
[fid,errmsg] = fopen(recordConfig, 'w'); % write only 
fprintf(fid,'FLOW_FIELD %f\n',flow_cgns_dt);
fprintf(fid,'\n');
fprintf(fid,'POINT_PARTICLE %f\n',point_cgns_dt);
fprintf(fid,'\n');
fprintf(fid,'PARAVIEW %f\n',paraview_dt);
fprintf(fid,'\n');
fprintf(fid,'PRECURSOR %f\n',precursor_dt);
fprintf(fid,'\n');
fprintf(fid,'RESTART %f\n',restart_dt);
fprintf(fid,'\n');
fprintf(fid,'RESTART_STOP %d\n',restart_stop_dt);
fprintf(fid,'\n');
fprintf(fid,'SCALAR %f\n',scalar_dt);
fprintf(fid,'\n');
fclose(fid);
 
type(recordConfig);
