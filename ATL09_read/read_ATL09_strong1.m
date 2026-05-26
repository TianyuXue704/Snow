function [initial_lat, initial_lon, column_od_asr,time] = read_ATL09_strong1(file)
initial_lat = h5read(file,'/profile_1/high_rate/latitude');
initial_lon = h5read(file,'/profile_1/high_rate/longitude');
column_od_asr = h5read(file,'/profile_1/high_rate/column_od_asr');
column_od_asr_qf = h5read(file,'/profile_1/high_rate/column_od_asr_qf');
time = 0:0.04:0.04*(length(initial_lat)-1);
ind_noise = find(column_od_asr_qf ~= 4);
if ~isempty(ind_noise)
    initial_lat(ind_noise) = [];
    initial_lon(ind_noise) = [];
    column_od_asr(ind_noise) = [];
    time(ind_noise) = [];
end
end







