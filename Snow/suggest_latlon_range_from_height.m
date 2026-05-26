function result = suggest_latlon_range_from_height(file_path, varargin)
% Suggest latitude/longitude filter ranges from height-based cloud screening.
%
% Example:
  % result = suggest_latlon_range_from_height( ...
  %     'E:\2026_project\Data\OIB\2019\0422\ATL03_mat_data\Day_h5data_strg3_...mat', ...
  %     'lat_bounds', [84.5 86.5], ...
  %     'lon_bounds', [-180 180], ...
  %     'save_dir', 'E:\2026_project\Data\OIB\2019\0422\result');

    p = inputParser;
    addRequired(p, 'file_path', @(x) ischar(x) || isstring(x));
    addParameter(p, 'save_dir', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'lat_bounds', [85.5 85.76], @(x) isnumeric(x) && numel(x) == 2);
    addParameter(p, 'lon_bounds', [-180 180], @(x) isnumeric(x) && numel(x) == 2);
    addParameter(p, 'time_bin_s', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'cloud_height_abs_m', 100, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'cloud_height_offset_m', 120, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'cloud_ratio_threshold', 0.02, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'quantile_low_high', [1 99], @(x) isnumeric(x) && numel(x) == 2);
    addParameter(p, 'min_segment_seconds', 5, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'save_csv', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'save_plot', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'verbose', true, @(x) islogical(x) && isscalar(x));
    parse(p, file_path, varargin{:});
    opt = p.Results;

    data = load(file_path, 'initial_t', 'initial_h', 'initial_lat', 'initial_lon');
    required_fields = {'initial_t', 'initial_h', 'initial_lat', 'initial_lon'};
    for k = 1:numel(required_fields)
        if ~isfield(data, required_fields{k})
            error('Missing variable "%s" in %s', required_fields{k}, file_path);
        end
    end

    t = double(data.initial_t(:));
    h = double(data.initial_h(:));
    lat = double(data.initial_lat(:));
    lon = double(data.initial_lon(:));

    valid = isfinite(t) & isfinite(h) & isfinite(lat) & isfinite(lon);
    t = t(valid);
    h = h(valid);
    lat = lat(valid);
    lon = lon(valid);

    lat_min = min(opt.lat_bounds);
    lat_max = max(opt.lat_bounds);
    lon_min = min(opt.lon_bounds);
    lon_max = max(opt.lon_bounds);
    geo_mask = lat >= lat_min & lat < lat_max & lon >= lon_min & lon < lon_max;

    t = t(geo_mask);
    h = h(geo_mask);
    lat = lat(geo_mask);
    lon = lon(geo_mask);

    if isempty(t)
        error('No points remain after the initial geo bounds filter.');
    end

    t_rel = t - t(1);
    h_median = median(h, 'omitnan');
    cloud_h_th = max(opt.cloud_height_abs_m, h_median + opt.cloud_height_offset_m);
    is_cloud_pt = h > cloud_h_th;

    t_min = min(t_rel);
    t_max = max(t_rel);
    t_edges = t_min:opt.time_bin_s:(t_max + opt.time_bin_s);
    if numel(t_edges) < 2
        t_edges = [t_min, t_min + opt.time_bin_s];
    end

    bin_idx = discretize(t_rel, t_edges);
    in_bin = ~isnan(bin_idx);
    n_bins = numel(t_edges) - 1;

    bin_count = accumarray(bin_idx(in_bin), 1, [n_bins, 1], @sum, 0);
    bin_cloud_count = accumarray(bin_idx(in_bin), double(is_cloud_pt(in_bin)), [n_bins, 1], @sum, 0);
    bin_h_max = accumarray(bin_idx(in_bin), h(in_bin), [n_bins, 1], @max, -Inf);

    bin_cloud_ratio = zeros(n_bins, 1);
    nz = bin_count > 0;
    bin_cloud_ratio(nz) = bin_cloud_count(nz) ./ bin_count(nz);

    is_cloud_bin = (bin_cloud_ratio >= opt.cloud_ratio_threshold) | (bin_h_max >= cloud_h_th);
    is_clean_bin = (~is_cloud_bin) & (bin_count > 0);

    segments = local_find_clean_segments(is_clean_bin, bin_count, t_edges);
    if isempty(segments)
        error('No clean segment found. Try relaxing cloud thresholds.');
    end

    duration_all = [segments.duration_s];
    keep_seg = duration_all >= opt.min_segment_seconds;
    if any(keep_seg)
        segments = segments(keep_seg);
    end

    [~, best_i] = max([segments.n_points]);
    best_seg = segments(best_i);

    candidate_name = "largest_clean";
    seg_list = best_seg;
    if any(is_cloud_bin)
        post_seg = local_post_cloud_segment(is_clean_bin, is_cloud_bin, bin_count, t_edges);
        if ~isempty(post_seg)
            candidate_name(end+1, 1) = "post_cloud";
            seg_list(end+1) = post_seg;
        end
    end

    q = sort(opt.quantile_low_high(:));
    n_candidate = numel(seg_list);
    t_start_s = zeros(n_candidate, 1);
    t_end_s = zeros(n_candidate, 1);
    duration_s = zeros(n_candidate, 1);
    n_points = zeros(n_candidate, 1);
    lat_min_raw = zeros(n_candidate, 1);
    lat_max_raw = zeros(n_candidate, 1);
    lon_min_raw = zeros(n_candidate, 1);
    lon_max_raw = zeros(n_candidate, 1);
    lat_q_low = zeros(n_candidate, 1);
    lat_q_high = zeros(n_candidate, 1);
    lon_q_low = zeros(n_candidate, 1);
    lon_q_high = zeros(n_candidate, 1);

    for i = 1:n_candidate
        seg = seg_list(i);
        in_seg = t_rel >= seg.t_start_s & t_rel < seg.t_end_s;
        if i == n_candidate
            in_seg = t_rel >= seg.t_start_s & t_rel <= seg.t_end_s;
        end

        lat_seg = lat(in_seg);
        lon_seg = lon(in_seg);

        t_start_s(i) = seg.t_start_s;
        t_end_s(i) = seg.t_end_s;
        duration_s(i) = seg.duration_s;
        n_points(i) = seg.n_points;

        lat_min_raw(i) = min(lat_seg);
        lat_max_raw(i) = max(lat_seg);
        lon_min_raw(i) = min(lon_seg);
        lon_max_raw(i) = max(lon_seg);

        lat_q = prctile(lat_seg, q);
        lon_q = prctile(lon_seg, q);
        lat_q_low(i) = lat_q(1);
        lat_q_high(i) = lat_q(2);
        lon_q_low(i) = lon_q(1);
        lon_q_high(i) = lon_q(2);
    end

    T = table(candidate_name, t_start_s, t_end_s, duration_s, n_points, ...
        lat_min_raw, lat_max_raw, lon_min_raw, lon_max_raw, ...
        lat_q_low, lat_q_high, lon_q_low, lon_q_high, ...
        'VariableNames', {'candidate', 't_start_s', 't_end_s', 'duration_s', 'n_points', ...
        'lat_min_raw', 'lat_max_raw', 'lon_min_raw', 'lon_max_raw', ...
        'lat_q_low', 'lat_q_high', 'lon_q_low', 'lon_q_high'});

    [folder, stem, ~] = fileparts(char(file_path));
    if isempty(opt.save_dir)
        save_dir = folder;
    else
        save_dir = char(opt.save_dir);
    end
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    csv_path = '';
    if opt.save_csv
        csv_path = fullfile(save_dir, [stem, '_latlon_suggestion.csv']);
        writetable(T, csv_path);
    end

    plot_path = '';
    if opt.save_plot
        fig = figure('Visible', 'off', 'Color', 'w');
        scatter(t_rel, h, 3, [0.75 0.75 0.75], 'filled');
        hold on;
        scatter(t_rel(is_cloud_pt), h(is_cloud_pt), 4, [0.9 0.1 0.1], 'filled');
        yl = ylim;
        for i = 1:n_candidate
            xline(t_start_s(i), '--', char(candidate_name(i)), 'LineWidth', 1.0);
            xline(t_end_s(i), '--', 'LineWidth', 1.0);
        end
        yline(cloud_h_th, '-.', 'Cloud Height Threshold');
        ylim(yl);
        grid on;
        xlabel('Relative Time (s)');
        ylabel('Height (m)');
        title('Height-Based Segment Suggestion');
        legend({'All points', 'Cloud-like high points'}, 'Location', 'best');
        plot_path = fullfile(save_dir, [stem, '_latlon_suggestion.png']);
        saveas(fig, plot_path);
        close(fig);
    end

    best_row = T(1, :);
    if opt.verbose
        fprintf('\nSuggested range from candidate: %s\n', best_row.candidate);
        fprintf('lat_min = %.6f; lat_max = %.6f;\n', best_row.lat_q_low, best_row.lat_q_high);
        fprintf('lon_min = %.6f; lon_max = %.6f;\n', best_row.lon_q_low, best_row.lon_q_high);
        fprintf('Cloud threshold height used: %.3f m\n', cloud_h_th);
    end

    result = struct();
    result.table = T;
    result.best_row = best_row;
    result.cloud_height_threshold_m = cloud_h_th;
    result.csv_path = csv_path;
    result.plot_path = plot_path;
    result.input_file = char(file_path);
    result.save_dir = save_dir;
end

function segments = local_find_clean_segments(is_clean_bin, bin_count, t_edges)
    edges_change = diff([false; is_clean_bin(:); false]);
    starts = find(edges_change == 1);
    ends = find(edges_change == -1) - 1;

    n = numel(starts);
    segments = repmat(struct('start_bin', 0, 'end_bin', 0, ...
        't_start_s', 0, 't_end_s', 0, 'duration_s', 0, 'n_points', 0), n, 1);

    for i = 1:n
        sb = starts(i);
        eb = ends(i);
        t0 = t_edges(sb);
        t1 = t_edges(eb + 1);
        segments(i).start_bin = sb;
        segments(i).end_bin = eb;
        segments(i).t_start_s = t0;
        segments(i).t_end_s = t1;
        segments(i).duration_s = t1 - t0;
        segments(i).n_points = sum(bin_count(sb:eb));
    end
end

function post_seg = local_post_cloud_segment(is_clean_bin, is_cloud_bin, bin_count, t_edges)
    post_seg = [];
    last_cloud = find(is_cloud_bin, 1, 'last');
    if isempty(last_cloud) || last_cloud >= numel(is_clean_bin)
        return;
    end

    right_clean = is_clean_bin((last_cloud + 1):end);
    rel_start = find(right_clean, 1, 'first');
    if isempty(rel_start)
        return;
    end

    sb = last_cloud + rel_start;
    eb = sb;
    while eb < numel(is_clean_bin) && is_clean_bin(eb + 1)
        eb = eb + 1;
    end

    t0 = t_edges(sb);
    t1 = t_edges(eb + 1);
    post_seg = struct( ...
        'start_bin', sb, ...
        'end_bin', eb, ...
        't_start_s', t0, ...
        't_end_s', t1, ...
        'duration_s', t1 - t0, ...
        'n_points', sum(bin_count(sb:eb)));
end
