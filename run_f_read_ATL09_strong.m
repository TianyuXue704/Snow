load pointing.mat
clc;clear;close all;
path = 'C:\Users\TianyuXue\Desktop\Np_correction\ATL09_read\ATL09_20230102055504_01941801_007_01.h5';
a = dir(path);
num_h5 = length(a);
for i = 1:num_h5
    h5file = [a(i).folder, '\', a(i).name];
    % Define the corresponding .mat file name by extracting the base name
    h5file_basename = erase(a(i).name, '.h5');
    folder = [a(i).folder,'\'];
    % Check if the .mat file already exists
    if exist(h5file_basename, 'file') 
        fprintf('The .mat file for %s already exists. Skipping...\n', h5file_basename);
    else
        % Process the HDF5 file using f_read_ATL03_strong function
        f_read_ATL09_strong(h5file);
        fprintf('succeed read %s\n',h5file_basename);
        close all;
    end
end


