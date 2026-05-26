
clear; close all;

load 'C:\Users\user\Desktop\ICESat-2 MC验证数据\验证区域1\ATL09_20190331003116_00270301_006_02.mat'

load 'C:\Users\user\Desktop\ICESat-2 MC验证数据\验证区域1\6902827_bbp_193.mat'

dis = (initial_lon - lon_choose).^2 + (initial_lat - lat_choose).^2;
[a,b] = min(dis);


% clear; close all;
% 
% load 'C:\Users\user\Desktop\ICESat-2 MC验证数据\验证区域2\ATL09_20190327144812_13620201_006_02.mat'
% 
% load 'C:\Users\user\Desktop\ICESat-2 MC验证数据\验证区域2\2902120_bbp_184.mat'
% 
% dis = (initial_lon - lon_choose).^2 + (initial_lat - lat_choose).^2;
% [a,b] = min(dis);