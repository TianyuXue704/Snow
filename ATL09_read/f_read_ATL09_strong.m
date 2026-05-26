function [] = f_read_ATL09_strong(h5file)
a = dir(h5file);
folder = [a.folder,'\'];
if ~exist(folder,'dir')
    mkdir(folder)
end

[initial_lat, initial_lon, column_od_asr,time] = read_ATL09_strong1(h5file);
file_strong1 = [folder,a.name(1:end-3),'.mat'];

save(file_strong1,'initial_lat', 'initial_lon', 'column_od_asr','time');
clearvars   initial_lat  initial_lon  column_od_asr time



