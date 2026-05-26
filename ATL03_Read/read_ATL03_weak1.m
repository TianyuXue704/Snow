% 读取.h5数据集中的强光束的时间、高度、经纬度信息
% sc_orient = 1, gtxR是强光束
% sc_orient = 0, gtxL是强光束
% 3个强光束内，读取优先级：3
function [initial_h, initial_t, initial_lat, initial_lon, ref_elev, signal_conf_ph,solar_elevation] = read_ATL03_weak1(file)
fileInfo = h5info(file);
groups = fileInfo.Groups;
for i = 1:length(groups)
    groupName{i} =groups(i).Name;
end
s_orient = h5read(file,'/orbit_info/sc_orient');
% 读取弱光束
if(s_orient == 1 & ismember('/gt1l', groupName))
    initial_h = h5read(file,'/gt1l/heights/h_ph');
    initial_t = h5read(file,'/gt1l/heights/delta_time');
    initial_lat = h5read(file,'/gt1l/heights/lat_ph'); %纬度
    initial_lon = h5read(file,'/gt1l/heights/lon_ph'); %经度
    ref_elev = h5read(file,'/gt1l/geolocation/ref_elev'); %地底角
    signal_conf_ph = h5read(file,'/gt1l/heights/signal_conf_ph'); % 光子置信度
    solar_elevation = h5read(file,'/gt1l/geolocation/solar_elevation'); % 太阳高度角
elseif(s_orient == 0 & ismember('/gt1r', groupName))
    initial_h = h5read(file,'/gt1r/heights/h_ph');
    initial_t = h5read(file,'/gt1r/heights/delta_time');
    initial_lat = h5read(file,'/gt1r/heights/lat_ph'); %纬度
    initial_lon = h5read(file,'/gt1r/heights/lon_ph'); %经度
    ref_elev = h5read(file,'/gt1r/geolocation/ref_elev'); %地底角
    signal_conf_ph = h5read(file,'/gt1r/heights/signal_conf_ph'); % 光子置信度
    solar_elevation = h5read(file,'/gt1r/geolocation/solar_elevation'); % 太阳高度角
else
    initial_h = [];
    initial_t = [];
    initial_lat = [];
    initial_lon = [];
    ref_elev= [];
    signal_conf_ph = [];
    solar_elevation = [];
    
end

end


