function [t_flow_str t_flow] = cgns_read_flow_time(casename);

%% Reads the flow times of the cgns files given in the output
% file, as specified by part-**.****.cgns etc (* is not
% indicative of sigfigs in this case, there can be more or less as
% specified in the record.config file

% t_flow_str - time written in cgns filename, i.e. flow-3.14159.cgns
% t_flow - time written in cgns file, using h5read

path = [casename '/output'];
od = cd(path);
contents = dir;
% Delete the . and .. directories from the contents
contents(1) = [];
contents(1) = [];
% Sort the contents by date added and remove unnecessary info
S = [contents(:).datenum];
[S,S] = sort(S);
contents = {contents(S).name};

% initialize time struct
t_flow_str = struct([]);

% Initialize k which increments to the next spot in the struct
j = 1;

% Look through directory contents and place flow times in the
% correct structure
for i = 1:length(contents)
    if strncmp('flow', contents{i},4) == 1
        time = sscanf(contents{i}, 'flow-%s');
        time = strrep(time, '.cgns', '');
        t_flow_str(j).time = time;
        tsol = '/Base/Zone0/Etc/Time/ data';
        t_path = [path '/' contents{i}];
        t_flow(j) = h5read(t_path, tsol);
        j = j+1;
    end
end

cd(od);
