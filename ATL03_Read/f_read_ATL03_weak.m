function [] = f_read_ATL03_weak(h5file)
a = dir(h5file);
folder = [a.folder,'\'];
if ~exist(folder,'dir')
    mkdir(folder)
end
%% 读取光子数据
[initial_h, initial_t, initial_lat, initial_lon, ref_elev, signal_conf_ph,solar_elevation] = read_ATL03_weak1(h5file);
% 判断昼夜（先剔除无效数据，再取平均）
solar_elevation(solar_elevation>90)=[];
solar_elev = mean(solar_elevation);
if(solar_elev > 0)
    flag = 'Day';  
else
    flag = 'Ngt'; 
end
% 保存数据
file_weak1 = [folder,flag,'_','h5data_weak1_',a.name(1:end-3),'.mat'];  %a.name(11:end-3) 是为了去掉processed
save(file_weak1,'initial_h','initial_t','initial_lat','initial_lon','ref_elev','signal_conf_ph','solar_elevation');
clearvars initial_h initial_t  initial_lat initial_lon ref_elev signal_conf_ph solar_elevation

[initial_h, initial_t, initial_lat, initial_lon, ref_elev, signal_conf_ph,solar_elevation] = read_ATL03_weak2(h5file);
solar_elevation(solar_elevation>90)=[];
solar_elev = mean(solar_elevation);
if(solar_elev > 0)
    flag = 'Day';  
else
    flag = 'Ngt'; 
end
file_weak2 = [folder,flag,'_','h5data_weak2_',a.name(1:end-3),'.mat'];
save(file_weak2,'initial_h','initial_t','initial_lat','initial_lon','ref_elev','signal_conf_ph','solar_elevation');
clearvars initial_h initial_t  initial_lat initial_lon ref_elev signal_conf_ph solar_elevation

[initial_h, initial_t, initial_lat, initial_lon, ref_elev, signal_conf_ph,solar_elevation] = read_ATL03_weak3(h5file);
solar_elevation(solar_elevation>90)=[];
solar_elev = mean(solar_elevation);
if(solar_elev > 0)
    flag = 'Day';  
else
    flag = 'Ngt'; 
end
file_weak3 = [folder,flag,'_','h5data_weak3_',a.name(1:end-3),'.mat'];
save(file_weak3,'initial_h','initial_t','initial_lat','initial_lon','ref_elev','signal_conf_ph','solar_elevation');
clearvars initial_h initial_t  initial_lat initial_lon ref_elev signal_conf_ph solar_elevation

end

