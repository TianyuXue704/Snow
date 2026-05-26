function summary = suggest_latlon_range_batch(input_root, varargin)
% Batch-run latitude/longitude suggestions for all MAT files in a folder.
%
% Example:
%   summary = suggest_latlon_range_batch( ...
%       'E:\2026_project\Data\OIB\2019\0412\ATL03_mat_data', ...
%       'save_dir', 'E:\2026_project\Data\OIB\2019\0412\latlon', ...
%       'lat_bounds', [85.5 85.76], ...
%       'lon_bounds', [-119.51 -113.12]);

    p = inputParser;
    p.KeepUnmatched = true;
    addRequired(p, 'input_root', @(x) (ischar(x) || isstring(x)) && isfolder(char(x)));
    addParameter(p, 'save_dir', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'recursive', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'file_pattern', '*.mat', @(x) ischar(x) || isstring(x));
    addParameter(p, 'skip_prefixes', {'dev'}, @(x) ischar(x) || isstring(x) || iscell(x));
    addParameter(p, 'summary_prefix', 'latlon_suggestion', @(x) ischar(x) || isstring(x));
    addParameter(p, 'continue_on_error', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'show_progress', true, @(x) islogical(x) && isscalar(x));
    parse(p, input_root, varargin{:});
    opt = p.Results;

    input_root = char(input_root);
    if isempty(opt.save_dir)
        save_dir = fullfile(input_root, 'latlon');
    else
        save_dir = char(opt.save_dir);
    end
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    file_pattern = char(opt.file_pattern);
    if opt.recursive
        file_list = dir(fullfile(input_root, '**', file_pattern));
    else
        file_list = dir(fullfile(input_root, file_pattern));
    end
    file_list = file_list(~[file_list.isdir]);

    skip_prefixes = local_to_cellstr(opt.skip_prefixes);
    if ~isempty(skip_prefixes)
        keep_mask = true(numel(file_list), 1);
        names = {file_list.name};
        for i = 1:numel(skip_prefixes)
            prefix = skip_prefixes{i};
            if isempty(prefix)
                continue;
            end
            keep_mask = keep_mask & ~strncmpi(names, prefix, numel(prefix))';
        end
        file_list = file_list(keep_mask);
    end

    if isempty(file_list)
        error('No MAT files found in %s', input_root);
    end

    forward_args = local_unmatched_to_pv(p.Unmatched, {'save_dir', 'verbose'});
    n_file = numel(file_list);
    candidate_tables = cell(n_file, 1);
    best_tables = cell(n_file, 1);
    error_file = strings(0, 1);
    error_path = strings(0, 1);
    error_message = strings(0, 1);

    for i = 1:n_file
        file_path = fullfile(file_list(i).folder, file_list(i).name);
        [~, file_stem, ~] = fileparts(file_list(i).name);

        if opt.show_progress
            fprintf('Lat/lon suggestion %d/%d: %s\n', i, n_file, file_list(i).name);
        end

        try
            result = suggest_latlon_range_from_height(file_path, ...
                'save_dir', save_dir, ...
                'verbose', false, ...
                forward_args{:});

            T = result.table;
            T = addvars(T, ...
                repmat(string(file_list(i).name), height(T), 1), ...
                repmat(string(file_stem), height(T), 1), ...
                repmat(string(file_path), height(T), 1), ...
                'Before', 1, ...
                'NewVariableNames', {'source_file', 'source_stem', 'source_path'});
            candidate_tables{i} = T;

            T_best = result.best_row;
            T_best = addvars(T_best, ...
                repmat(string(file_list(i).name), height(T_best), 1), ...
                repmat(string(file_stem), height(T_best), 1), ...
                repmat(string(file_path), height(T_best), 1), ...
                'Before', 1, ...
                'NewVariableNames', {'source_file', 'source_stem', 'source_path'});
            best_tables{i} = T_best;
        catch ME
            error_file(end + 1, 1) = string(file_list(i).name);
            error_path(end + 1, 1) = string(file_path);
            error_message(end + 1, 1) = string(ME.message);

            if opt.show_progress
                fprintf(2, '  Failed: %s\n', ME.message);
            end

            if ~opt.continue_on_error
                rethrow(ME);
            end
        end
    end

    valid_candidate = ~cellfun(@isempty, candidate_tables);
    valid_best = ~cellfun(@isempty, best_tables);

    if any(valid_candidate)
        T_candidates = vertcat(candidate_tables{valid_candidate});
        all_candidates_csv = fullfile(save_dir, [char(opt.summary_prefix), '_all_candidates.csv']);
        writetable(T_candidates, all_candidates_csv);
    else
        T_candidates = table();
        all_candidates_csv = '';
    end

    if any(valid_best)
        T_best = vertcat(best_tables{valid_best});
        best_csv = fullfile(save_dir, [char(opt.summary_prefix), '_best.csv']);
        writetable(T_best, best_csv);
    else
        T_best = table();
        best_csv = '';
    end

    if ~isempty(error_file)
        T_errors = table(error_file, error_path, error_message, ...
            'VariableNames', {'source_file', 'source_path', 'error_message'});
        error_csv = fullfile(save_dir, [char(opt.summary_prefix), '_errors.csv']);
        writetable(T_errors, error_csv);
    else
        T_errors = table();
        error_csv = '';
    end

    if opt.show_progress
        fprintf('Completed %d files: %d succeeded, %d failed.\n', ...
            n_file, sum(valid_best), numel(error_file));
    end

    summary = struct();
    summary.save_dir = save_dir;
    summary.best_table = T_best;
    summary.all_candidates_table = T_candidates;
    summary.error_table = T_errors;
    summary.best_csv = best_csv;
    summary.all_candidates_csv = all_candidates_csv;
    summary.error_csv = error_csv;
end

function values = local_to_cellstr(x)
    values = cellstr(string(x(:)));
end

function pv = local_unmatched_to_pv(unmatched, exclude_names)
    if nargin < 2
        exclude_names = {};
    end

    names = fieldnames(unmatched);
    if ~isempty(exclude_names)
        keep_mask = ~ismember(lower(names), lower(cellstr(string(exclude_names(:)))));
        names = names(keep_mask);
    end

    pv = cell(1, 2 * numel(names));
    for i = 1:numel(names)
        pv{2 * i - 1} = names{i};
        pv{2 * i} = unmatched.(names{i});
    end
end
