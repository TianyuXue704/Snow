function [] = f_read_ATL03_strong(h5file,SavePath)
 a = dir(h5file);
folder = SavePath;
if ~exist(folder,'dir')
    mkdir(folder)
end


[initial_h, initial_t, initial_lat, initial_lon, ref_elev,solar_elevation,initial_distance,initial_index,Et,back_ground_save, signal_conf_ph] = read_ATL03_strong3(h5file);
solar_elevation(solar_elevation>90)=[];
solar_elev = mean(solar_elevation);
if(solar_elev > 0)
    flag = 'Day';  
else
    flag = 'Ngt'; 
end
file_strong1 = [folder,flag,'_','h5data_strg3_',a.name(1:end-3),'.mat'];
save(file_strong1,'initial_h','initial_t','initial_lat','initial_lon','ref_elev','solar_elevation','initial_distance','initial_index','Et','back_ground_save','signal_conf_ph');
clearvars initial_h  initial_t  initial_lat  initial_lon  ref_elev solar_elevation initial_distance initial_index Et back_ground_save signal_conf_ph
end

