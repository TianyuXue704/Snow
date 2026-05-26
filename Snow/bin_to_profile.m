function [S, Depth, data_S] = bin_to_profile(bin, data_bin, ind_peak, flag)
[m, n] = size(bin);

% --- Parameters ---
delta_h = 0.15; 
delta_t = 0.0002;  
dt = 5;             % Integration duration (e.g., 2 seconds)
step_t = 0.5;       % Sliding step (e.g., 0.1 seconds)

N1 = 50; 
N2 = min(300, m - ind_peak(flag)); 
Depth = -(N1 * delta_h) : delta_h : N2 * delta_h;
range_h = ind_peak(flag) - N1 : ind_peak(flag) + N2;

% --- Initialization ---
S = [];
data_S = [];
ThStep = 0.2 * (dt / delta_t); % Minimum points required to validate a window

% Start sliding
current_start_time = data_bin(3, 1);
last_total_time = data_bin(3, end);

while (current_start_time + dt) <= last_total_time
    % 1. Find indices within the current [T, T + dt] window
    idx_in_window = find(data_bin(3, :) >= current_start_time & ...
                         data_bin(3, :) < (current_start_time + dt));
    
    % 2. Check if window has enough data points
    if length(idx_in_window) > ThStep
        tempS = mean(bin(range_h, idx_in_window), 2);
        tempdataS = mean(data_bin(:, idx_in_window), 2);
        
        S = [S tempS];
        data_S = [data_S tempdataS];
    end
    
    % 3. Slide the window forward by step_t (0.1s)
    current_start_time = current_start_time + step_t;
    
    % Safety break if we exceed the data range
    if current_start_time > last_total_time
        break;
    end
end

if isempty(S)
    return
end
end