clc; clear; close all;

input_root = 'E:\2026_project\Data\OIB\2019\0406\ATL03_mat_data';
save_dir = 'E:\2026_project\Data\OIB\2019\0406\suggest_latlon';

lat_bounds = [82.416107, 82.518077];
lon_bounds = [-8.869804, -8.180369];

summary = suggest_latlon_range_batch(input_root, ...
    'save_dir', save_dir, ...
    'lat_bounds', lat_bounds, ...
    'lon_bounds', lon_bounds, ...
    'save_plot', false);

disp(summary.best_csv);
