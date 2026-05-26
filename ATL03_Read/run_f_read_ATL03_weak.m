% 读取一个路径下的所有.h5文件
clc;clear;close all;
path = 'F:\论文画图数据\与Lu同数据包之南极冰缘浮游植物\Ngt\*.h5';
a = dir(path);
num_h5 = length(a);
for i = 1:num_h5
    h5file = [a(i).folder,'\',a(i).name];
    f_read_ATL03_weak(h5file);
    close all;
end