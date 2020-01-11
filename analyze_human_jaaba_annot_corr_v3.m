function [recall_rates_per_cat_all_behav, score_post_analysis_struct] = analyze_human_jaaba_annot_corr_v3(flymat_id, behav_sel, override_jab_list, is_NoRel, bout_matches_mat, plot_bar, plot_violin, plot_box, img_formats)
    load('common-params-annot-analysis.mat', ...
         'annot_file', 'behav_list', 'behav_shorthands', 'jab_list', 'groundtruth_movie_list');
    [~, annot_file_name, ~] = fileparts(annot_file);
    annot_file_name_elements = strsplit(annot_file_name, '-');
    annot_file_id = annot_file_name_elements{end};
    load(sprintf('FLYMAT_HumanAnnotation_%s.mat', annot_file_id), 'flymatHumanAnnot');
    flymat_file = strcat('FLYMAT_MNL-KA JAABA training samples_', flymat_id, '.mat'); 
    load(flymat_file, 'flymatAll');

%     % For later use of flymat_id, shorten the names a bit
%     flymat_id_elems = strsplit(flymat_id, '_');
%     flymat_id_elems{1} = regexp(flymat_id_elems{1} ,'[^0-9]+', 'match');
%     flymat_id = strjoin([flymat_id_elems{:}], '_'); 
    
    movie_name_pairs = cellfun(@(c) strsplit(c, '\'), groundtruth_movie_list, 'UniformOutput', false);
    if is_NoRel
        behav_list = {'LungeNewNoRel', 'WingExtNoRel', 'HeadbuttNew'};
    end
    if ~isempty(override_jab_list)
        for i=1:length(override_jab_list)
            jab_list{override_jab_list(i).behav_sel} = override_jab_list(i).jab_name;
        end
    end
    
    init_bout_matches_all_args = [behav_list; cell(1, length(behav_list))];
    bout_matches_all = struct(init_bout_matches_all_args{:});

    if ~strcmp(bout_matches_mat, '')
        load(bout_matches_mat, 'bout_matches_all');
    else
        % Compile all the data we need 
        for i=1:length(behav_list)
            if ~ismember(i, behav_sel)
                continue;
            end
            % Create dynamic field name for this current behaviour
            t0s = strcat(behav_shorthands{i}, '_t0s');
            t1s = strcat(behav_shorthands{i}, '_t1s');
            scores = strcat(behav_shorthands{i}, '_scores');
            combined_score = strcat(behav_shorthands{i}, '_combined_score');
            t4s = strcat(behav_shorthands{i}, '_t4s');
            t5s = strcat(behav_shorthands{i}, '_t5s');
            is_duplicate = strcat(behav_shorthands{i}, '_is_duplicate');
            annotators = strcat(behav_shorthands{i}, '_annotators');

            % Create structure for bout matches in one type of behaviour
            % How we deal with no match situations:
            % 1) For a JAABA false negative, we average JAABA score in
            %    [annot_intsct_start, annot_intsct_end]. 
            % 2) For a JAABA false positive, set the annot_score as 0
            bout_matches_fields = {'movie', 'fly', ...
                'annot_union_start', 'annot_union_end', ...
                'annot_intsct_start', 'annot_intsct_end', 'annot_score', 'annot_score_pair', 'annotators'...
                'jaaba_bout_start', 'jaaba_bout_end', 'jaaba_score_avg', ...
                'jaaba_score_min', 'jaaba_score_max', ...
                'jaaba_score_avg_normed', 'jaaba_score_min_normed', 'jaaba_score_max_normed', ...
                'multi_match', 'virtual_jaaba_match', 'invalid_match'};
            init_bout_matches_args = [bout_matches_fields; cell(1,length(bout_matches_fields))];
            bout_matches = struct(init_bout_matches_args{:});
            bout_matches(1) = [];

            % Extract annotated bouts
            % Currently, we take the union of two or more annotations to calculate
            % the start and end of the annotated bout. This is rather liberal, and
            % is aimed to make matching with JAABA bouts eaiser. 
            % When one annotator did not label a bout, while another did, NaN
            % values will appear, but they do not affect the min/max operation
            bout_idx = uint16(1);
            nan_init_cell_arr = cell(1,length(bout_matches_fields));
            nan_init_cell_arr(:) = {nan};
            init_bout_matches_args = [bout_matches_fields; nan_init_cell_arr];
            bout_match = struct(init_bout_matches_args{:});
            for j=1:length(flymatHumanAnnot)
                for k=1:length(flymatHumanAnnot(j).(combined_score))
                    bout_match.movie = flymatHumanAnnot(j).movie;
                    bout_match.fly = flymatHumanAnnot(j).fly;
                    bout_match.annot_intsct_start = max(flymatHumanAnnot(j).(t0s)(k,:));
                    bout_match.annot_intsct_end = min(flymatHumanAnnot(j).(t1s)(k,:));
                    bout_match.annot_union_start = min(flymatHumanAnnot(j).(t0s)(k,:));
                    bout_match.annot_union_end = max(flymatHumanAnnot(j).(t1s)(k,:));
                    bout_match.annot_score = flymatHumanAnnot(j).(combined_score)(k);
                    bout_match.annot_score_pair = flymatHumanAnnot(j).(scores)(k, :);
                    bout_match.annotators = flymatHumanAnnot(j).(annotators);
                    bout_match_struct = bout_match;
                    if bout_match.annot_intsct_start > bout_match.annot_intsct_end % If two non-overlapping bouts are recorded as overlapping
                        bout_match_struct(2) = bout_match;
                        for l=1:2
                            bout_match_struct(l).annot_intsct_start = flymatHumanAnnot(j).(t0s)(k,l);
                            bout_match_struct(l).annot_intsct_end = flymatHumanAnnot(j).(t1s)(k,l);
                            bout_match_struct(l).annot_union_start = bout_match_struct(l).annot_intsct_start;
                            bout_match_struct(l).annot_union_end = bout_match_struct(l).annot_intsct_end;
                            bout_match_struct(l).annot_score = flymatHumanAnnot(j).(scores)(k, l);
                            if l == 1
                                bout_match_struct(l).annot_score_pair = [flymatHumanAnnot(j).(scores)(k, l), nan];
                            else
                                bout_match_struct(l).annot_score_pair = [nan, flymatHumanAnnot(j).(scores)(k, l)];
                            end
                        end
                    end
                    for l=1:length(bout_match_struct)
                        bout_match = bout_match_struct(l);
                        if ~any(flymatHumanAnnot(j).(is_duplicate)(k,:))
                            bout_matches(bout_idx) = bout_match;
                            bout_idx = bout_idx + 1;
                        elseif bitand(flymatHumanAnnot(j).(is_duplicate)(k,1), flymatHumanAnnot(j).(is_duplicate)(k,2))
                            continue;
                        else
                            % If one annotation corresponds to multiple rows in
                            % the annotation matrix/EXCEL
                            for m=1:2 % Deal additionally with the rare possibility that immediate previous row is faulty and split into 2
                                if ~isempty(intersect(bout_matches(bout_idx-m).annot_union_start:bout_matches(bout_idx-m).annot_union_end, ...
                                        bout_match.annot_union_start:bout_match.annot_union_end))
                                    annot_intsct_prev = arrayfun(@(s, e) s:e, bout_matches(bout_idx-m).annot_intsct_start, bout_matches(bout_idx-m).annot_intsct_end-1, 'UniformOutput', false);
                                    annot_intsct_prev = horzcat(annot_intsct_prev{:});
                                    intsct_overlap = intersect(annot_intsct_prev, ...
                                        bout_match.annot_intsct_start:bout_match.annot_intsct_end-1); 
                                    % If one annotator's annotation is broken into two
                                    % by another annotator
                                    if isempty(intsct_overlap)
                                        annot_intsct_frames = unique(horzcat(annot_intsct_prev, ...
                                            bout_match.annot_intsct_start:bout_match.annot_intsct_end-1));
                                        annot_intsct_boundaries = find(diff(annot_intsct_frames)>1);
                                        annot_intsct_starts = annot_intsct_frames([1, annot_intsct_boundaries+1]);
                                        annot_intsct_ends = annot_intsct_frames([annot_intsct_boundaries, end]);
                                        annot_intsct_ends = annot_intsct_ends + 1;
                                        %bout_matches(bout_idx-m).annot_score = ...
                                        %    max(bout_matches(bout_idx-m).annot_score, bout_match.annot_score);
                                        bout_matches(bout_idx-m).annot_score = [bout_matches(bout_idx-m).annot_score, bout_match.annot_score];
                                        bout_matches(bout_idx-m).annot_score_pair = [bout_matches(bout_idx-m).annot_score_pair; bout_match.annot_score_pair];
                                        bout_matches(bout_idx-m).annot_intsct_start = annot_intsct_starts; 
                                        bout_matches(bout_idx-m).annot_intsct_end = annot_intsct_ends;  
                                    % If one overlapping annotation is recorded in
                                    % two rows
                                    elseif any(arrayfun(@(l) l == length(intsct_overlap), ...
                                            [bout_matches(bout_idx-m).annot_intsct_end - bout_matches(bout_idx-m).annot_intsct_start, ...
                                            bout_match.annot_intsct_end - bout_match.annot_intsct_start]))
                                        bout_matches(bout_idx-m).annot_score = max([bout_matches(bout_idx-m).annot_score, bout_match.annot_score]);
                                        bout_matches(bout_idx-m).annot_score_pair = [bout_matches(bout_idx-m).annot_score_pair; bout_match.annot_score_pair];
                                        bout_matches(bout_idx-m).annot_score_pair(any(isnan(bout_matches(bout_idx-m).annot_score_pair), 2), :) = [];
                                        bout_matches(bout_idx-m).annot_intsct_start = max([bout_matches(bout_idx-m).annot_intsct_start, bout_match.annot_intsct_start]); 
                                        bout_matches(bout_idx-m).annot_intsct_end = min([bout_matches(bout_idx-m).annot_intsct_end, bout_match.annot_intsct_end]);  
                                        bout_matches(bout_idx-m).annot_union_start = min([bout_matches(bout_idx-m).annot_union_start, bout_match.annot_union_start]);
                                        bout_matches(bout_idx-m).annot_union_end = max([bout_matches(bout_idx-m).annot_union_end, bout_match.annot_union_end]); 
                                    else
                                        fprintf('Undocumented situation for bout overlap appears\n');
                                        pause;
                                    end
                                    break;
                                end
                            end
                        end
                    end
                end
            end

            % Iterate over annotated flies, find them in the JAABA flymat, match
            % JAABA bouts one-by-one with annotated bouts. If no match is found,
            % add the JAABA bout as a new entry in bout_matches
            for j=1:length(flymatHumanAnnot)
                fly_idx = find_fly_in_flymat(flymatAll, flymatHumanAnnot(j).movie, ...
                    flymatHumanAnnot(j).fly, true);
                candidates = find_fly_in_flymat(bout_matches, ...
                        flymatHumanAnnot(j).movie, flymatHumanAnnot(j).fly, false);

                for k=1:length(flymatAll(fly_idx).(t4s))
                    jaaba_bout_start = flymatAll(fly_idx).(t4s)(k);
                    jaaba_bout_end = flymatAll(fly_idx).(t5s)(k);     
                    matched = 0;
                    for m=1:length(candidates)
                        % Declare a match if there is any overlap between annotated
                        % bout and JAABA bout
                        if (jaaba_bout_start <= bout_matches(candidates(m)).annot_union_start ...
                                && jaaba_bout_end > bout_matches(candidates(m)).annot_union_start) ...
                            || ...
                            (jaaba_bout_start > bout_matches(candidates(m)).annot_union_start ...
                                && jaaba_bout_start < bout_matches(candidates(m)).annot_union_end)
                            if isnan(bout_matches(candidates(m)).jaaba_bout_start)
                                bout_matches(candidates(m)).jaaba_bout_start = ...
                                    jaaba_bout_start; 
                                bout_matches(candidates(m)).jaaba_bout_end = ...
                                    jaaba_bout_end;
                            else
                                bout_matches(candidates(m)).jaaba_bout_start = ...
                                    [bout_matches(candidates(m)).jaaba_bout_start, jaaba_bout_start]; 
                                bout_matches(candidates(m)).jaaba_bout_end = ...
                                    [bout_matches(candidates(m)).jaaba_bout_end, jaaba_bout_end];
                            end
                            bout_matches(candidates(m)).virtual_jaaba_match = false;
                            bout_matches(candidates(m)).invalid_match = false;
                            if ~matched
                                matched = candidates(m);
                                bout_matches(candidates(m)).multi_match = false;
                            else
                                % bout_matches(matched).multi_match = true;
                                bout_matches(candidates(m)).multi_match = true;
                            end
                        end
                    end
                    % If no match is found for a JAABA bout, add it to bout_matches
                    if ~matched
                        bout_match = struct(init_bout_matches_args{:});
                        bout_match.movie = flymatHumanAnnot(j).movie;
                        bout_match.fly = flymatHumanAnnot(j).fly;
                        bout_match.annot_union_start = nan;
                        bout_match.annot_union_end = nan;
                        bout_match.annot_score = 0;
                        bout_match.annot_score_pair = [0, 0];
                        bout_match.jaaba_bout_start = jaaba_bout_start;
                        bout_match.jaaba_bout_end = jaaba_bout_end;
                        bout_match.virtual_jaaba_match = false;
                        bout_match.invalid_match = false;
                        bout_matches(bout_idx) = bout_match;
                        bout_idx = bout_idx + 1;
                    end
                end
            end

            % If no match is found for an existing human annotation, use
            % intersection of human annotation (consensus) as a virtual JAABA bout 
            for j=1:length(bout_matches)
                if isnan(bout_matches(j).jaaba_bout_start)
                    bout_matches(j).jaaba_bout_start = bout_matches(j).annot_intsct_start;
                    bout_matches(j).jaaba_bout_end = bout_matches(j).annot_intsct_end;
                    bout_matches(j).virtual_jaaba_match = true;
                end
            end

            % Sort the bouts by movie field
            bout_matches_cell = struct2cell(bout_matches);
            bout_matches_cell = permute(squeeze(bout_matches_cell), [2,1]);
            bout_matches_cell = sortrows(bout_matches_cell, 1:4);
            bout_matches_cell = permute(bout_matches_cell, [2,1]);
            bout_matches = cell2struct(bout_matches_cell, bout_matches_fields, 1)';

            % Merge human annotated bouts if their unions intersect. 
            % This deals with cases that were missed by the
            % duplication process in human annotation creation, because the
            % combined annotation score is not recorded properly
            rerun_merge_bout = true; 
            while rerun_merge_bout
                bouts_to_remove = [];
                for j=2:length(bout_matches)
                    annot_union_intsct = intersect(bout_matches(j-1).annot_union_start:bout_matches(j-1).annot_union_end-1, bout_matches(j).annot_union_start:bout_matches(j).annot_union_end-1);
                    if ~isempty(annot_union_intsct)
                        % Fix the combined score that is not correctly
                        % recorded
                        if any(isnan(bout_matches(j-1).annot_score_pair))
                            bout_matches(j-1).annot_score_pair(isnan(bout_matches(j-1).annot_score_pair)) = ...
                                bout_matches(j).annot_score_pair(isnan(bout_matches(j-1).annot_score_pair));
                        end
                        if any(isnan(bout_matches(j).annot_score_pair))
                            bout_matches(j).annot_score_pair(isnan(bout_matches(j).annot_score_pair)) = ...
                                bout_matches(j-1).annot_score_pair(isnan(bout_matches(j).annot_score_pair));    
                        end
                        bout_matches(j-1).annot_score = [sum(bout_matches(j-1).annot_score_pair, 2)', sum(bout_matches(j).annot_score_pair, 2)];
                        bout_matches(j-1).annot_score_pair = [bout_matches(j-1).annot_score_pair; bout_matches(j).annot_score_pair];
                       
                        
                        bout_matches(j-1).annot_union_start = min(bout_matches(j-1).annot_union_start, bout_matches(j).annot_union_start);
                        bout_matches(j-1).annot_union_end = max(bout_matches(j-1).annot_union_end, bout_matches(j).annot_union_end);
                        
                        annot_intsct_prev = arrayfun(@(s, e) s:e, bout_matches(j-1).annot_intsct_start, bout_matches(j-1).annot_intsct_end-1, 'UniformOutput', false);
                        annot_intsct_prev = horzcat(annot_intsct_prev{:});
                        intsct_overlap = intersect(annot_intsct_prev, ...
                            bout_matches(j).annot_intsct_start:bout_matches(j).annot_intsct_end-1); 
                        % If one annotator's annotation is broken into two
                        % by another annotator
                        if isempty(intsct_overlap)
                            annot_intsct_frames = unique(horzcat(annot_intsct_prev, ...
                                        bout_matches(j).annot_intsct_start:bout_matches(j).annot_intsct_end-1));
                            annot_intsct_boundaries = find(diff(annot_intsct_frames)>1);
                            annot_intsct_starts = annot_intsct_frames([1, annot_intsct_boundaries+1]);
                            annot_intsct_ends = annot_intsct_frames([annot_intsct_boundaries, end]);
                            annot_intsct_ends = annot_intsct_ends + 1;
                            bout_matches(j-1).annot_intsct_start = annot_intsct_starts; 
                            bout_matches(j-1).annot_intsct_end = annot_intsct_ends; 
                        % If one overlapping annotation is recorded in
                        % two rows
                        elseif any(arrayfun(@(l) l == length(intsct_overlap), ...
                                [bout_matches(j-1).annot_intsct_end - bout_matches(j-1).annot_intsct_start, ...
                                bout_matches(j).annot_intsct_end - bout_matches(j).annot_intsct_start]))
                            bout_matches(j-1).annot_score = bout_matches(j-1).annot_score(1);
                            bout_matches(j-1).annot_score_pair = bout_matches(j-1).annot_score_pair(1, :);
                            bout_matches(j-1).annot_intsct_start = max([bout_matches(j-1).annot_intsct_start, bout_matches(j).annot_intsct_start]); 
                            bout_matches(j-1).annot_intsct_end = min([bout_matches(j-1).annot_intsct_end, bout_matches(j).annot_intsct_end]);
                        else
                            fprintf('Undocumented situation for bout overlap appears\n');
                            pause;
                        end
                        bouts_to_remove = [bouts_to_remove, j];
                    end
                    
                    % Keep only the max combined score if we have only one
                    % intersection segment. 
                    if length(bout_matches(j-1).annot_intsct_start) < length(bout_matches(j-1).annot_score)
                        bout_matches(j-1).annot_score = max(bout_matches(j-1).annot_score);
                    end
                end
                if isempty(bouts_to_remove)
                    rerun_merge_bout = false;
                end
                bout_matches(bouts_to_remove) = [];
                fprintf('Removed or combined %d intersecting annotation bouts for behavior %s\n', length(bouts_to_remove), behav_list{i});
            end
            
            % Iterate over the movies, load scores.mat, 
            % calculate the avg, min and max score for each JAABA bout entry 
            % in bout_matches
            for j=1:length(groundtruth_movie_list)
                load(fullfile(groundtruth_movie_list{j}, ...
                    strcat(movie_name_pairs{j}{2},'_JAABA'), ...
                    strcat('scores_', erase(behav_list{i}, {'_','s'}),'.mat')), 'allScores');
                bouts = uint16(find_fly_in_flymat(bout_matches, movie_name_pairs{j}{1}, [], false));
                for k=1:length(bouts)
                    jaaba_bout_start = bout_matches(bouts(k)).jaaba_bout_start;
                    jaaba_bout_end = bout_matches(bouts(k)).jaaba_bout_end;

                    if (~isempty(jaaba_bout_start) && ~isempty(jaaba_bout_end)) || ...
                            (~isnan(jaaba_bout_start) && ~isnan(jaaba_bout_end))
                        jaaba_bout_montage = arrayfun(@(j_bout_start, j_bout_end) j_bout_start:j_bout_end-1, jaaba_bout_start, jaaba_bout_end, 'UniformOutput', false);
                        scores_slice = allScores.scores{bout_matches(bouts(k)).fly}(cell2mat(jaaba_bout_montage));
                        bout_matches(bouts(k)).jaaba_score_avg = mean(scores_slice);
                        bout_matches(bouts(k)).jaaba_score_min = min(scores_slice);
                        bout_matches(bouts(k)).jaaba_score_max = max(scores_slice);
                        bout_matches(bouts(k)).jaaba_score_avg_normed = bout_matches(bouts(k)).jaaba_score_avg/allScores.scoreNorm;
                        bout_matches(bouts(k)).jaaba_score_min_normed = bout_matches(bouts(k)).jaaba_score_min/allScores.scoreNorm;
                        bout_matches(bouts(k)).jaaba_score_max_normed = bout_matches(bouts(k)).jaaba_score_max/allScores.scoreNorm;
                    end
                end
            end

            bout_matches_all.(behav_list{i}) = bout_matches;
        end

        save(strcat('bout_matches_', flymat_id, '.mat'), 'bout_matches_all');
    end
    
    % Plot human annotation score vs. JAABA score for valid bout matches, for
    % each behaviour type
    stats_to_plot = {'avg', 'min', 'max'};
    human_annot_scores = 0:6;
    xtick_label_cell_arr = strtrim(cellstr(num2str(human_annot_scores')));
    xtick_label_cell_arr{end+1} = 'training';
    recall_prec_args = [behav_list; cell(1, length(behav_list))];
    recalls = struct(recall_prec_args{:});
    precisions = struct(recall_prec_args{:});
    recall_rates_per_cat_all_behav = zeros(length(behav_list), length(human_annot_scores) - 1);
    score_post_analysis_struct = struct('annot_score', [], 'jaaba_score', [], 'false_negat_idxs', []);
    score_post_analysis_struct(1) = [];
    jaaba_score_quantile_struct_all = struct(recall_prec_args{:});
    for i=1:length(behav_list)
        if ~ismember(i, behav_sel)
            continue;
        end
        score_benchmrk_train = find_train_bout_jaaba_score_stats(...
            behav_list{i}, jab_list{i});
%         annot_score = uint8([bout_matches_all.(behav_list{i})(bout_masks{i}).annot_score]);
%         annot_score = cellfun(@(scores) max(scores), {bout_matches_all.(behav_list{i})(bout_masks{i}).annot_score});
%         false_negat_idxs = [bout_matches_all.(behav_list{i})(bout_masks{i}).virtual_jaaba_match];
        annot_score_max = cellfun(@(scores) max(scores), {bout_matches_all.(behav_list{i}).annot_score});
        annot_score_all = horzcat(bout_matches_all.(behav_list{i}).annot_score);
        associated_annot_per_bout_match = arrayfun(@(s) length(s.annot_score), bout_matches_all.(behav_list{i}));
        false_negat_idxs = [bout_matches_all.(behav_list{i}).virtual_jaaba_match];
        false_negat_idxs_w_duplicate = repelem(false_negat_idxs, associated_annot_per_bout_match);
        
%         % Calculate recalls without counting bout matches with
%         % multi_match=1
%         multi_match_idxs = [bout_matches_all.(behav_list{i})(bout_masks{i}).multi_match];
%         multi_match_idxs(isnan(multi_match_idxs)) = 0;
%         virt_annot_score_no_mm = annot_score(bitand(false_negat_idxs, ~multi_match_idxs)); 
%         true_annot_score_no_mm = annot_score(bitand(~false_negat_idxs, ~multi_match_idxs));
%         false_negat_accum_no_mm = hist(virt_annot_score_no_mm, human_annot_scores);
%         positive_accum_no_mm = hist(true_annot_score_no_mm, human_annot_scores);
%         recalls.(behav_list{i}) = sum(positive_accum_no_mm(2:end))/(sum(false_negat_accum_no_mm(2:end))+sum(positive_accum_no_mm(2:end)));
%         collect_accum_no_mm = [positive_accum_no_mm(2:end)', false_negat_accum_no_mm(2:end)'];
%         recall_rates = collect_accum_no_mm(:,1)./sum(collect_accum_no_mm, 2);
                   
        % Calculate precision and recall
        virt_annot_score_max = annot_score_max(false_negat_idxs); 
        true_annot_score_max = annot_score_max(~false_negat_idxs);
        virt_annot_score_all = annot_score_all(false_negat_idxs_w_duplicate); 
        true_annot_score_all = annot_score_all(~false_negat_idxs_w_duplicate);
        false_negat_accum = hist(virt_annot_score_all, human_annot_scores);
        positive_accum = hist(true_annot_score_all, human_annot_scores);
        collect_accum = [positive_accum(2:end)', false_negat_accum(2:end)'];
        test_violin_norm = max([false_negat_accum, positive_accum]);
        positive_train_size = length(score_benchmrk_train.behav.avg);
        negat_train_size = length(score_benchmrk_train.none_behav.avg);
        train_violin_norm = max([positive_train_size, negat_train_size]);
        
        recalls.(behav_list{i}) = sum(positive_accum(2:end))/(sum(false_negat_accum(2:end))+sum(positive_accum(2:end)));
        recall_rates = collect_accum(:,1)./sum(collect_accum, 2);
        precisions.(behav_list{i}) = sum(positive_accum(2:end))/sum(positive_accum);
        
        for k=1:length(stats_to_plot)
            field_to_plot = strcat('jaaba_score_', stats_to_plot{k}, '_normed');
            jaaba_score = [bout_matches_all.(behav_list{i}).(field_to_plot)];
            virt_jaaba_score = jaaba_score(false_negat_idxs);
            true_jaaba_score = jaaba_score(~false_negat_idxs);
            positive_train_jaaba_score = cell2mat(score_benchmrk_train.behav.(stats_to_plot{k}));
            negative_train_jaaba_score = cell2mat(score_benchmrk_train.none_behav.(stats_to_plot{k}));
            
            if k == 1
                score_post_analysis_struct(i).annot_score = annot_score_max; 
                score_post_analysis_struct(i).jaaba_score = jaaba_score; 
                score_post_analysis_struct(i).false_negat_idxs = false_negat_idxs; 
            end
            % =====
            % Box plot
            % =====
            if plot_box
                figure();
                hold on;

                boxplot(true_jaaba_score, true_annot_score_max, 'Positions', double(human_annot_scores)-0.15, 'PlotStyle', 'compact', 'Widths', 0.2, 'Colors', 'b');
                boxplot(virt_jaaba_score, virt_annot_score_max, 'Positions', double(human_annot_scores(2:end))+0.15, 'PlotStyle', 'compact', 'Widths', 0.2, 'Colors', 'm');
                boxplot(positive_train_jaaba_score, 'Positions', 7-0.15, 'PlotStyle', 'compact', 'Widths', 0.2, 'Colors', [0, 0.5, 0.5]);
                boxplot(negative_train_jaaba_score, 'Positions', 7+0.15, 'PlotStyle', 'compact', 'Widths', 0.2, 'Colors', [0.75, 0.25, 0]);
                boxs = findobj(gca, 'Tag', 'Box');
                legend(boxs([10, 3, 2, 1]),  {'true or false positive', 'false negative', 'positive training', 'negative training'},...
                'NumColumns', 2, 'Location', 'southwest');
                offset_hack = 0;
                
                all_plotted_scores = [true_jaaba_score, virt_jaaba_score, positive_train_jaaba_score, negative_train_jaaba_score];
                axis([-0.5+offset_hack, 7.5+offset_hack, min(all_plotted_scores)-1.2, ...
                    max(all_plotted_scores)+0.3]);
                line([-0.5, 6.5]+offset_hack, repmat(mean(cell2mat(score_benchmrk_train.behav.(stats_to_plot{k}))), [1, 2]), ...
                    'Color', 'k', 'LineStyle', '--', 'DisplayName', 'avg. over all positive{\it training} bouts');
                line([-0.5, 6.5]+offset_hack, repmat(mean(cell2mat(score_benchmrk_train.none_behav.(stats_to_plot{k}))), [1, 2]), ...
                    'Color', 'k', 'LineStyle', '-.', 'DisplayName', 'avg. over all negative{\it training} bouts');
                ylim_tmp = ylim; 
                train_bkg = area([6.5+offset_hack, 7.5+offset_hack], [ylim_tmp(2),ylim_tmp(2)], ylim_tmp(1), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                train_bkg.Annotation.LegendInformation.IconDisplayStyle = 'off';
                title({strcat(stats_to_plot{k}, '. of JAABA score vs. human combined score for test bouts'), ...
                    strcat('behaviour:', erase(behav_list{i},'_'))});
                xlabel('Human annotation combined score \in [0,6]');
                ylabel(strcat(stats_to_plot{k}, '. of JAABA score in a bout'));
                set(gca, 'xtick', [human_annot_scores, 7]++offset_hack);
                set(gca, 'xticklabel', xtick_label_cell_arr);
                hold off;
                set(gcf,'renderer','Painters');
                for j=1:length(img_formats)
                    if contains(img_formats{j}, 'eps')
                        img_format = 'epsc';
                    else
                        img_format = 'png';
                    end
                    saveas(gcf, strcat(strcat('jaaba_score_', stats_to_plot{k}, '_normed'), ...
                    '_vs_annot_combined-',erase(behav_list{i},'_'), '-box-', flymat_id, '.', erase(img_format, 'c')), img_format);
                end
            end 
            % =====
            % Violin plot
            % =====
            if plot_violin
                figure();
                hold on;
                v1 = violinplot(true_jaaba_score, double(true_annot_score_max)-0.3, 'ViolinColor', [0, 0, 1], 'ShowData', false, 'ViolinNorm', test_violin_norm);
                v2 = violinplot([-10, virt_jaaba_score], double([0, virt_annot_score_max]+0.3), 'ViolinColor', [1, 0, 1], 'ShowData', false, 'ViolinNorm', test_violin_norm);
                v3 = violinplot([repmat(-10, 1, 7), positive_train_jaaba_score], ...
                    [0:6, repmat(7, 1, length(positive_train_jaaba_score))], 'ViolinColor', [0, 0.5, 0.5], 'ShowData', false, 'ViolinNorm', train_violin_norm);
                v4 = violinplot([repmat(-10, 1, 7), negative_train_jaaba_score], ...
                    [0:6, repmat(7, 1, length(negative_train_jaaba_score))], 'ViolinColor', [0.75, 0.25, 0], 'ShowData', false, 'ViolinNorm', train_violin_norm);
                set(gca, 'xtick', [human_annot_scores, 7]+1);
                set(gca, 'xticklabel', xtick_label_cell_arr);
                legend([v1(end).ViolinPlot, v2(2).ViolinPlot, v3(end).ViolinPlot, v4(end).ViolinPlot], {'true or false positive', 'false negative', 'positive training', 'negative training'},...
                'NumColumns', 2, 'Location', 'southwest');
                offset_hack = 1;
                
                jaaba_score_quantile_struct = struct();
                for m=human_annot_scores
                    jaaba_score_quantile_struct(m+1).data_type = 'test';
                    jaaba_score_quantile_struct(m+1).annot_score = m;
                    scores_one_cat = true_jaaba_score(true_annot_score_max == m);
                    if ~isempty(scores_one_cat)
                        jaaba_score_quantile_struct(m+1).jaaba_posit_quantiles = quantile(scores_one_cat, [0.1, 0.25, 0.5, 0.75, 0.9]);
                    else
                        jaaba_score_quantile_struct(m+1).jaaba_posit_quantiles = nan(5, 1);
                    end
                    scores_one_cat = virt_jaaba_score(virt_annot_score_max == m);
                    if ~isempty(scores_one_cat)
                        jaaba_score_quantile_struct(m+1).jaaba_negat_quantiles = quantile(scores_one_cat, [0.1, 0.25, 0.5, 0.75, 0.9]);
                    else
                        jaaba_score_quantile_struct(m+1).jaaba_negat_quantiles = nan(5, 1);
                    end
                    scores_one_cat = horzcat(true_jaaba_score(true_annot_score_max == m), virt_jaaba_score(virt_annot_score_max == m));
                    if ~isempty(scores_one_cat)
                        jaaba_score_quantile_struct(m+1).jaaba_combined_quantiles = quantile(scores_one_cat, [0.1, 0.25, 0.5, 0.75, 0.9]);
                    else
                        jaaba_score_quantile_struct(m+1).jaaba_combined_quantiles = nan(5, 1);
                    end
                    
                end
                jaaba_score_quantile_struct(8).data_type = 'train';
                jaaba_score_quantile_struct(8).annot_score = nan;
                jaaba_score_quantile_struct(8).jaaba_posit_quantiles = quantile(positive_train_jaaba_score, [0.1, 0.25, 0.5, 0.75, 0.9]);
                jaaba_score_quantile_struct(8).jaaba_negat_quantiles = quantile(negative_train_jaaba_score, [0.1, 0.25, 0.5, 0.75, 0.9]);
                jaaba_score_quantile_struct(8).jaaba_combined_quantiles = quantile(horzcat(positive_train_jaaba_score, negative_train_jaaba_score), [0.1, 0.25, 0.5, 0.75, 0.9]);
                jaaba_score_quantile_struct_all.(behav_list{i}).(stats_to_plot{k}) = jaaba_score_quantile_struct; 
                
                all_plotted_scores = [true_jaaba_score, virt_jaaba_score, positive_train_jaaba_score, negative_train_jaaba_score];
                axis([-0.5+offset_hack, 7.5+offset_hack, min(all_plotted_scores)-1.2, ...
                    max(all_plotted_scores)+0.3]);
                line([-0.5, 6.5]+offset_hack, repmat(mean(cell2mat(score_benchmrk_train.behav.(stats_to_plot{k}))), [1, 2]), ...
                    'Color', 'k', 'LineStyle', '--', 'DisplayName', 'avg. over all positive{\it training} bouts');
                line([-0.5, 6.5]+offset_hack, repmat(mean(cell2mat(score_benchmrk_train.none_behav.(stats_to_plot{k}))), [1, 2]), ...
                    'Color', 'k', 'LineStyle', '-.', 'DisplayName', 'avg. over all negative{\it training} bouts');
                ylim_tmp = ylim; 
                train_bkg = area([6.5+offset_hack, 7.5+offset_hack], [ylim_tmp(2),ylim_tmp(2)], ylim_tmp(1), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                train_bkg.Annotation.LegendInformation.IconDisplayStyle = 'off';
                title({strcat(stats_to_plot{k}, '. of JAABA score vs. human combined score for test bouts'), ...
                    strcat('behaviour:', erase(behav_list{i},'_'))});
                xlabel('Human annotation combined score \in [0,6]');
                ylabel(strcat(stats_to_plot{k}, '. of JAABA score in a bout'));
                set(gca, 'xtick', [human_annot_scores, 7]++offset_hack);
                set(gca, 'xticklabel', xtick_label_cell_arr);
                hold off;
                set(gcf,'renderer','Painters');
                for j=1:length(img_formats)
                    if contains(img_formats{j}, 'eps')
                        img_format = 'epsc';
                    else
                        img_format = 'png';
                    end
                    saveas(gcf, strcat(strcat('jaaba_score_', stats_to_plot{k}, '_normed'), ...
                        '_vs_annot_combined-',erase(behav_list{i},'_'), '-violin-', flymat_id, '.', erase(img_format, 'c')), img_format);
                end
            end
        end 

        % =====
        % Bar plot
        % =====
        if plot_bar
            figure();
            hold on;
            yyaxis left
            bar_plt = bar(1:6, collect_accum, 'stacked', 'FaceAlpha', 0.5);
            false_positive_bar = bar(0, positive_accum(1), 'FaceColor', 'b', 'FaceAlpha', 0.5);
            yyaxis right
            plot(1:6, recall_rates*100, 'LineWidth', 1.5);   
            yyaxis left
            sum_collect_accum = sum(collect_accum, 2);
            % Put annotations of the binned recall rate
            for j=0:1
                binned_recall = uint8(sum(collect_accum(1+j*3:3+j*3, 1))/sum(sum_collect_accum(1+j*3:3+j*3))*100);
                curr_lim = ylim;
                bin_label_yloc = max(sum_collect_accum(1+j*3:3+j*3))+20;
                axes_colors = get(gcf, 'defaultAxesColorOrder');
                drawbrace([1+j*3, bin_label_yloc], [3+j*3, bin_label_yloc], 7, 'Color', axes_colors(2,:));
                text(2+j*3, bin_label_yloc+40*(curr_lim(2)/500), strcat(num2str(binned_recall), '%'), 'Color', axes_colors(2,:));
            end
            title({'JAABA bout classification count', ...
                strcat('behaviour:', erase(behav_list{i},'_')), ...
                sprintf('recall %d%%, precision %d%%', uint8(recalls.(behav_list{i})*100), uint8(precisions.(behav_list{i})*100))});
            xlabel('Human annotation combined score \in [1,6]');
            yyaxis left
            ylabel('Human annotation bout classification count');
            curr_lim = ylim;
            ylim([curr_lim(1), curr_lim(2)*1.2]);
            yyaxis right
            ylabel('recall in %');
            ylim([0, 100]);
            curr_lim = xlim;
            xlim([curr_lim(1)+0.3, curr_lim(2)-0.3])
            set(bar_plt, {'FaceColor'}, {'b'; 'm'});
            legend({'true or false positives', 'false negatives'}, 'Location', 'northwest');
            hold off;
            set(gcf,'renderer','Painters');
            for j=1:length(img_formats)
                if contains(img_formats{j}, 'eps')
                    img_format = 'epsc';
                else
                    img_format = 'png';
                end
                saveas(gcf, strcat('human_bout_classification_count-',erase(behav_list{i},'_'), '-bar-', flymat_id, '.', erase(img_format, 'c')), ...
                    img_format);
            end
        end
        recall_rates_per_cat_all_behav(i, :) = recall_rates; 
        save(sprintf('jaaba_score_quantile_struct_all-%s.mat', flymat_id), 'jaaba_score_quantile_struct_all');
    end
end