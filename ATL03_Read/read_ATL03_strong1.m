
function [initial_h, initial_t, initial_lat, initial_lon, ref_elev,solar_elevation,initial_distance,initial_index,Et,back_ground_save] = read_ATL03_strong1(file)
    fileInfo = h5info(file);
    groups = fileInfo.Groups;
    for i = 1:length(groups)
        groupName{i} =groups(i).Name;
    end
    s_orient = h5read(file,'/orbit_info/sc_orient');
    if(s_orient == 1 && ismember('/gt1r', groupName))
        initial_h = h5read(file,'/gt1r/heights/h_ph');
        initial_t = h5read(file,'/gt1r/heights/delta_time');
        initial_lat = h5read(file,'/gt1r/heights/lat_ph'); 
        initial_lon = h5read(file,'/gt1r/heights/lon_ph'); %
        initial_distance = lldistkm(initial_lat(1),initial_lon(1),initial_lat,initial_lon);
        ref_elev = h5read(file,'/gt1r/geolocation/ref_elev'); %?????
        signal_conf_ph = h5read(file,'/gt1r/heights/signal_conf_ph'); %?????
        solar_elevation = h5read(file,'/gt1r/geolocation/solar_elevation'); % ?0?9???????
        alongtrack_distance = h5read(file,'/gt1r/heights/dist_ph_along');%??????????segment
        segment_length = h5read(file,'/gt1r/geolocation/segment_length');%???segment?????
        Et = h5read(file,'/ancillary_data/atlas_engineering/transmit/tx_pulse_energy');%???????????????????????
        qualityph = h5read(file,'/gt1r/heights/quality_ph');
        index_beg = h5read(file,'/gt1r/geolocation/ph_index_beg');
        back_ground_rate = h5read(file,'/gt1r/bckgrd_atlas/bckgrd_counts_reduced');
        back_height = h5read(file,'/gt1r/bckgrd_atlas/bckgrd_int_height_reduced');
        back_ground_rate = double(back_ground_rate);
        back_height = double(back_height);
        mean_rate  = back_ground_rate./back_height;%average photon rate
        back_time = h5read(file,'/gt1r/bckgrd_atlas/delta_time');
        t_start = initial_t(1);
        initial_t = initial_t - t_start;
        back_time = back_time - back_time(1);
        back_ground_save = interp1(back_time, mean_rate, initial_t, 'nearest');%插值实现噪声
        index_null = find(isnan(back_ground_save), 1);
        back_ground_save(index_null:end) = back_ground_save(index_null-1);
        Et = max(Et(1,1),Et(1,2))*10^6;
        index_beg = single(index_beg);
        ref_elev = ref_elev(index_beg~=0);
        index_beg = index_beg(index_beg~=0);
        initial_beg = 1:length(initial_t)';
        [a,initial_index] = histc(initial_beg,index_beg);
        initial_index = initial_index';
        index1 = find(initial_index == 0,1);
        initial_index(index1:end) = initial_index(index1-1);
        ind_quality = find(qualityph == 0);
        initial_h = initial_h(ind_quality);
        initial_lat = initial_lat(ind_quality);
        initial_lon = initial_lon(ind_quality);
        initial_index = initial_index(ind_quality);
        initial_distance = initial_distance(ind_quality);
        initial_t = initial_t(ind_quality);
       % ind_confi = find(signal_conf_ph(2,:)<=0)';
   
    elseif(s_orient == 0 & ismember('/gt1l', groupName))
        initial_h = h5read(file,'/gt1l/heights/h_ph');
        initial_t = h5read(file,'/gt1l/heights/delta_time');
        initial_lat = h5read(file,'/gt1l/heights/lat_ph'); 
        initial_lon = h5read(file,'/gt1l/heights/lon_ph'); %
        initial_distance = lldistkm(initial_lat(1),initial_lon(1),initial_lat,initial_lon);
        ref_elev = h5read(file,'/gt1l/geolocation/ref_elev'); %?????
        signal_conf_ph = h5read(file,'/gt1l/heights/signal_conf_ph'); %?????
        solar_elevation = h5read(file,'/gt1l/geolocation/solar_elevation'); % ?0?9???????
        alongtrack_distance = h5read(file,'/gt1l/heights/dist_ph_along');%??????????segment
        segment_length = h5read(file,'/gt1l/geolocation/segment_length');%???segment?????
        Et = h5read(file,'/ancillary_data/atlas_engineering/transmit/tx_pulse_energy');%???????????????????????
        qualityph = h5read(file,'/gt1l/heights/quality_ph');
        index_beg = h5read(file,'/gt1l/geolocation/ph_index_beg');
        back_ground_rate = h5read(file,'/gt1l/bckgrd_atlas/bckgrd_counts_reduced');
        back_height = h5read(file,'/gt1l/bckgrd_atlas/bckgrd_int_height_reduced');
        back_ground_rate = double(back_ground_rate);
        back_height = double(back_height);
        mean_rate  = back_ground_rate./back_height;%average photon rate
        back_time = h5read(file,'/gt1l/bckgrd_atlas/delta_time');
        t_start = initial_t(1);
        initial_t = initial_t - t_start;
        back_time = back_time - back_time(1);
        back_ground_save = interp1(back_time, mean_rate, initial_t, 'nearest');%插值实现噪声
        index_null = find(isnan(back_ground_save), 1);
        back_ground_save(index_null:end) = back_ground_save(index_null-1);
        Et = max(Et(1,1),Et(1,2))*10^6;
        index_beg = single(index_beg);
        ref_elev = ref_elev(index_beg~=0);
        index_beg = index_beg(index_beg~=0);
        initial_beg = 1:length(initial_t)';
        [a,initial_index] = histc(initial_beg,index_beg);
        initial_index = initial_index';
        index1 = find(initial_index == 0,1);
        initial_index(index1:end) = initial_index(index1-1);
        ind_quality = find(qualityph == 0);
        initial_h = initial_h(ind_quality);
        initial_lat = initial_lat(ind_quality);
        initial_lon = initial_lon(ind_quality);
        initial_index = initial_index(ind_quality);
        initial_distance = initial_distance(ind_quality);
        initial_t = initial_t(ind_quality);
    else
        initial_h = [];
        initial_t = [];
        initial_lat = [];
        initial_lon = [];
        ref_elev= [];
        solar_elevation = [];
        initial_distance = [];
        initial_index = [];
        Et = [];
        back_ground_save = [];

    end



    end






