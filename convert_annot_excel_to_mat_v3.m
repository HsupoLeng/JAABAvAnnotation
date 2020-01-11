% Convert from human annotation EXCEL to FLYMAT fashioned MATLAB .mat file
% Parameters: where to extract data, and what behvaiours to look at
load('common-params-annot-analysis.mat');

movie_name_pairs = cellfun(@(c) strsplit(c, '\'), groundtruth_movie_list, 'UniformOutput', false);

% Create structure conforming with the format of flymatAll
init_fields = {'movie', 'fly', 'L_t0s', 'L_t1s', 'L_scores', 'L_combined_score', 'L_annotators', 'L_is_duplicate'...
    'HB_t0s', 'HB_t1s', 'HB_scores', 'HB_combined_score', 'HB_annotators', 'HB_is_duplicate', 'WE_t0s', 'WE_t1s', ...
    'WE_scores', 'WE_combined_score', 'WE_annotators', 'WE_is_duplicate'};
init_flymat_args = [init_fields; cell(1, length(init_fields))];
flymatHumanAnnot = struct(init_flymat_args{:});
flymatHumanAnnot(1) = []; % Remove the first, dummy entry in flymatHumanAnnot

% Iterate over behaviours, corresponding to sheets in the original excel
for i=1:length(behav_list)
    [~,~,raw_xls] = xlsread(annot_file, behav_list{i});
    raw_xls = cellfun(@(x) num2str(x), raw_xls, 'UniformOutput', false); % Convert all fields to str 
    raw_xls = raw_xls(sum(cellfun(@(c) ~strcmp(c, 'NaN'), raw_xls), 2)>0, :); % Remove rows where entries are all NaN
    
    % Create dynamic field name for this current behaviour
    t0s = strcat(behav_shorthands{i}, '_t0s');
    t1s = strcat(behav_shorthands{i}, '_t1s');
    scores = strcat(behav_shorthands{i}, '_scores');
    combined_score = strcat(behav_shorthands{i}, '_combined_score');
    annotators = strcat(behav_shorthands{i}, '_annotators');
    is_duplicate = strcat(behav_shorthands{i}, '_is_duplicate');
    
    % Iterate over movies
    for j=1:length(groundtruth_movie_list)
        % Extract fields for one movie
        movie_start_idx = find(contains(raw_xls, movie_name_pairs{j}{2}));
        if j<length(groundtruth_movie_list)
            movie_end_idx = find(contains(raw_xls, movie_name_pairs{j+1}{2}));
        else
            movie_end_idx = numel(raw_xls);
        end
        [movie_start_row, movie_start_col] = ind2sub(size(raw_xls), movie_start_idx);
        [~, movie_end_col] = ind2sub(size(raw_xls), movie_end_idx);
        movie_fields = raw_xls(movie_start_row:size(raw_xls,1), movie_start_col:movie_end_col-1);
        
        % Extract the number and indices of flies annotated in one movie
        fly_start_indices = find(contains(movie_fields, 'Fly'));
        fly_end_indices = find(contains(movie_fields, 'Human cobmbined score'));
        num_of_flies = {movie_fields{fly_start_indices}};
        fly_numbers = regexp(num_of_flies, '[0-9]+', 'match');
        fly_numbers = cellfun(@(c) str2double(c), fly_numbers);
        
        % Iterate over annotated flies in this movie
        for k=1:length(fly_numbers)
            [fly_start_row, fly_start_col] = ind2sub(size(movie_fields), fly_start_indices(k));
            [~, fly_end_col] = ind2sub(size(movie_fields), fly_end_indices(k));
            fly_fields = movie_fields(fly_start_row:size(movie_fields,1), fly_start_col:fly_end_col);
            valid_fly_fields = fly_fields(sum(cellfun(@(c) ~strcmp(c, 'NaN'), fly_fields), 2)>0, :);
            annot_start_indices = find(contains(valid_fly_fields, {'Start', 'start'}), 1, 'first');
            annotator_names = valid_fly_fields(annot_start_indices(1)-1, 1:end-1);
            annotator_names = annotator_names(~cellfun(@(s) strcmp(s, 'NaN'), annotator_names));
            if length(annotator_names) > 2
                annotator_names = annotator_names(1:2); % We ignore the third annotator
            end
            annot_fields = valid_fly_fields(annot_start_indices(1):end, :);   
            
            rows_nonzero_total_score_mask = ~isnan(cellfun(@(c) str2double(c), annot_fields(:, end)));
            true_positive_fields = annot_fields(rows_nonzero_total_score_mask, :);
            cols_to_keep = boolean(ones(1, size(true_positive_fields, 2)));
            if ~isempty(annot_fields) && any(strcmp(annot_fields(1,:), 'Duration'))
                cols_to_keep(4:4:end) = 0; % Discard duration column, if present
            end
            true_positive_fields = cellfun(@(c) str2double(c), true_positive_fields(:, cols_to_keep));            
            true_positive_fields(isnan(true_positive_fields)) = 0;
            
            % Fill in cells where one bout annotated by annotator A is
            % labelled as two bouts by another annotator
            rows_to_discard = [];
            is_duplicate_mat = false(size(true_positive_fields,1), floor(size(true_positive_fields, 2)/3));
            for l=1:size(true_positive_fields, 1)
                if sum(true_positive_fields(l, 3:3:end-1)) ~= true_positive_fields(l, end)
                    scores_col = 3:3:size(true_positive_fields, 2)-1;
                    for m=1:length(scores_col)
                        if true_positive_fields(l, scores_col(m))
                            continue
                        end
                        history_offset = 1;
                        while ~true_positive_fields(l, scores_col(m))
                            true_positive_fields(l, scores_col(m)-2:scores_col(m)) = ...
                                true_positive_fields(l-history_offset, scores_col(m)-2:scores_col(m));
                            history_offset = history_offset + 1;
                        end
                        is_duplicate_mat(l, scores_col(m)/3) = true; 
                        if sum(true_positive_fields(l, 3:3:end-1)) == true_positive_fields(l, end)
                            break;
                        end
                    end
                    bout_ends = true_positive_fields(l, 2:3:end-2);
                    bout_ends = sort(bout_ends, 'descend');
                    if max(true_positive_fields(l, 1:3:end-3)) > bout_ends(2)
                        fprintf('Bouts dont overlap: behaviour %s, contents from that row\n', behav_shorthands{i});
                        true_positive_fields(l, :)
                        % Deal with two instances where we should pull from
                        % one row downward
                        history_offset = -1;
                        true_positive_fields(l, scores_col(m)-2:scores_col(m)) = ...
                                true_positive_fields(l-history_offset, scores_col(m)-2:scores_col(m));
                        is_duplicate_mat(l, scores_col(m)/3) = false;
                        is_duplicate_mat(l+1, scores_col(m)/3) = true;
                        fprintf('Fixed by pulling from the immediate next row\n');
                    end
                    
                    if sum(true_positive_fields(l, 3:3:end-1)) ~= true_positive_fields(l, end)
                        fprintf('Things still dont add up: contents from that row\n');
                        true_positive_fields(l, :)
                        rows_to_discard = [rows_to_discard, l];
                    end
                end
            end
            if ~isempty(rows_to_discard)
                for z=1:length(fliplr(rows_to_discard))
                    fprintf('Discarding this row\n');
                    true_positive_fields(rows_to_discard(z), :)
                    true_positive_fields(rows_to_discard(z), :) = [];
                    is_duplicate_mat(rows_to_discard(z)+1,:) = false(1, size(is_duplicate_mat, 2));
                    is_duplicate_mat(rows_to_discard(z), :) = [];
                    
                end
            end
            
            % Remove rows where only the third annotator labelled a bout
            if size(annot_fields, 2) > 9
                extra_scores = true_positive_fields(:, end-1);
                extra_scores(isnan(extra_scores)) = 0;
            else
                extra_scores = zeros(size(true_positive_fields, 1), 1);
            end
            only_third_annotator_mask = true_positive_fields(:, end) - extra_scores == 0; 
            true_positive_fields(:, end) = true_positive_fields(:, end) - extra_scores;
            true_positive_fields = true_positive_fields(~only_third_annotator_mask, :);
            is_duplicate_mat = is_duplicate_mat(~only_third_annotator_mask, :);
            % Change all zeros back to nan so we don't break the
            % analyze_human_jaaba_annot_corr_v2.m code
            true_positive_fields(true_positive_fields==0) = nan;
            
            % Find the entry in flymatHumanAnnot using (movie, fly) as ID
            fly_idx = find_fly_in_flymat(flymatHumanAnnot, movie_name_pairs{j}{1}, ...
                fly_numbers(k), false);
            if isempty(fly_idx)
                fly_idx = size(flymatHumanAnnot, 2) + 1;
                flymatHumanAnnot(fly_idx).movie = movie_name_pairs{j}{1};
                flymatHumanAnnot(fly_idx).fly = fly_numbers(k);
            end
            
            % Record the annotations into the entry
            % If there are more than two annotators, we take only the first
            % two annotators' annotation
            flymatHumanAnnot(fly_idx).(t0s) = true_positive_fields(:, 1:3:4);
            flymatHumanAnnot(fly_idx).(t1s) = true_positive_fields(:, 2:3:5);
            flymatHumanAnnot(fly_idx).(scores) = true_positive_fields(:,3:3:6);
            flymatHumanAnnot(fly_idx).(combined_score) = ...
                true_positive_fields(:,end);
            flymatHumanAnnot(fly_idx).(annotators) = annotator_names;
            flymatHumanAnnot(fly_idx).(is_duplicate) = is_duplicate_mat(:,[1,2]); 
        end
    end  
end

[~, annot_file_name, ~] = fileparts(annot_file);
annot_file_name_elements = strsplit(annot_file_name, '-');
annot_file_id = annot_file_name_elements{end};
save(sprintf('FLYMAT_HumanAnnotation_%s.mat', annot_file_id), 'flymatHumanAnnot');