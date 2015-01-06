function [t_part_str t_part] = cgns_read_part_time(casename)

%% Reads the part times of the cgns files given in the output
% file, as specified by part-**.****.cgns etc (* is not
% indicative of sigfigs in this case, there can be more or less as
% specified in the record.config file

% t_part_str is the time written in the cgns file name, i.e. part-12.31415.cgns
% t_part is the actual time written in the cgns file, found using h5read

path = [casename '/output'];
od = cd(path);
contents = dir;
% Delete the . and .. directories from contents
contents(1) = [];
contents(1) = [];
% Sort the contents by date added and remove unnecessary info
S = [contents(:).datenum];
[S,S] = sort(S);
contents = {contents(S).name};

% Initialize time struct matrix
t_part_str = struct([]);

% Initialize j which increments to the next spot in the struct
j = 1;

% Look through directory contents and place flow times in the
% correct structure
for i = 1:length(contents)
    if strncmp('part', contents{i},4) == 1
        time = sscanf(contents{i}, 'part-%s');
        time = strrep(time, '.cgns', '');
        t_part_str(j).time = time;
        tsol = '/Base/Zone0/Etc/Time/ data';
        t_path = [path '/' contents{i}];
        t_part(j) = h5read(t_path, tsol);
        j = j+1;
    end
end

cd(od);
