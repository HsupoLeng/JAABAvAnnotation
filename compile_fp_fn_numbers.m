function compile_fp_fn_numbers(behav_sel, match_sel_str, is_frame_matches)
    if is_frame_matches
        load(strcat('frame_matches_', match_sel_str, '.mat'), 'frame_matches_all');
        matches_all = frame_matches_all; 
    else
        load(strcat('bout_matches_', match_sel_str, '.mat'), 'bout_matches_all');
        matches_all = bout_matches_all; 
    end
    load('common-params-annot-analysis.mat', 'behav_list');
    
    if ~isempty(regexp(match_sel_str, '.*woRel.*', 'once'))
        behav_list = {'LungeNewNoRel', 'WingExtNoRel', 'HeadbuttNew'};
    end
    true_false_count_all_args = [behav_list; cell(1, length(behav_list))];
    true_false_count_all = struct(true_false_count_all_args{:});
    for i=1:length(behav_list)
        if ~ismember(i, behav_sel)
            continue;
        end
        
%         annot_score = uint8(full(cellfun(@(scores) max(scores), {matches_all.(behav_list{i}).annot_score})));
        annot_score_all = horzcat(matches_all.(behav_list{i}).annot_score);
        associated_annot_per_bout_match = arrayfun(@(s) length(s.annot_score), matches_all.(behav_list{i}));
        false_positive_idxs = cellfun(@(scores) max(scores), {matches_all.(behav_list{i}).annot_score}) == 0;
        false_positive_idxs_w_duplicate = repelem(false_positive_idxs, associated_annot_per_bout_match);
        if is_frame_matches
            false_negat_idxs = [matches_all.(behav_list{i}).is_false_negative];
        else
            false_negat_idxs = [matches_all.(behav_list{i}).virtual_jaaba_match];    
        end
        false_negat_idxs_w_duplicate = repelem(false_negat_idxs, associated_annot_per_bout_match);
        
        human_annot_scores = 0:6; 
        virt_annot_score = annot_score_all(false_negat_idxs_w_duplicate); 
        true_annot_score = annot_score_all(~false_negat_idxs_w_duplicate);
        false_negat_accum = hist(virt_annot_score, human_annot_scores);
        positive_accum = hist(true_annot_score, human_annot_scores);
        true_counts = [0, positive_accum(2:end)];
        false_counts = [positive_accum(1), false_negat_accum(2:end)];
        per_category_struct = struct('score', nan, 'false_positives', nan, 'false_negatives', nan, 'true_positives', nan, 'true_total', nan);
        for j=1:length(false_counts)+1
            if j==1
                per_category_struct(j).score = j-1;
                per_category_struct(j).false_positives = false_counts(j);
                per_category_struct(j).false_negatives = 0; 
                per_category_struct(j).true_positives = 0;
                per_category_struct(j).true_total = 0;
                per_category_struct(j).recall = nan;
            elseif j==length(false_counts)+1
                per_category_struct(j).score = 'summary';
                per_category_struct(j).false_positives = sum([per_category_struct(1:j-1).false_positives]);
                per_category_struct(j).false_negatives = sum([per_category_struct(1:j-1).false_negatives]); 
                per_category_struct(j).true_positives = sum([per_category_struct(1:j-1).true_positives]);
                per_category_struct(j).true_total = sum([per_category_struct(1:j-1).true_total]);
                per_category_struct(j).recall = nan;
            else
                per_category_struct(j).score = j-1;
                per_category_struct(j).false_positives = 0; 
                per_category_struct(j).false_negatives = false_counts(j); 
                per_category_struct(j).true_positives = true_counts(j);
                per_category_struct(j).true_total = per_category_struct(j).false_negatives ...
                    + per_category_struct(j).true_positives;
                per_category_struct(j).recall = per_category_struct(j).true_positives / per_category_struct(j).true_total;
            end
        end
        true_false_count_all.(behav_list{i}).per_category_count = per_category_struct;
        true_false_count_all.(behav_list{i}).recall = sum([per_category_struct(1:end-1).true_positives])/sum([per_category_struct(1:end-1).true_total]);
        true_false_count_all.(behav_list{i}).precision = sum([per_category_struct(1:end-1).true_positives])/...
            (sum([per_category_struct(1:end-1).true_positives]) + sum([per_category_struct(1:end-1).false_positives]));

        
        all_sources = {matches_all.(behav_list{i})(:).movie};
        all_sources_w_duplicate = repelem(all_sources, associated_annot_per_bout_match);
        false_positive_sources_w_duplicate = all_sources_w_duplicate(false_positive_idxs_w_duplicate); 
        false_negative_sources_w_duplicate = all_sources_w_duplicate(false_negat_idxs_w_duplicate); 
        all_sources_set = unique(all_sources);
        
        per_movie_count = struct('movie', '', 'false_positives', nan, 'false_negatives', nan, 'true_positives', nan, 'true_total', nan, ...
            'false_positive_fly_frame', [], 'false_negative_fly_frame', []);
        for j=1:length(all_sources_set)+1
            if j<length(all_sources_set)+1
                per_movie_count(j).movie = all_sources_set{j};
                this_movie_mask = strcmp(all_sources, all_sources_set(j));
                per_movie_count(j).false_positives = nnz(strcmp(false_positive_sources_w_duplicate, all_sources_set{j}));
                per_movie_count(j).false_negatives = nnz(strcmp(false_negative_sources_w_duplicate, all_sources_set{j})); 
                per_movie_count(j).true_total = nnz(strcmp(all_sources_w_duplicate, all_sources_set{j})) - per_movie_count(j).false_positives;
                per_movie_count(j).true_positives = per_movie_count(j).true_total - per_movie_count(j).false_negatives; 
            else
                per_movie_count(j).movie = 'summary';
                per_movie_count(j).false_positives = sum([per_movie_count(1:j-1).false_positives]);
                per_movie_count(j).false_negatives = sum([per_movie_count(1:j-1).false_negatives]); 
                per_movie_count(j).true_total = sum([per_movie_count(1:j-1).true_total]);
                per_movie_count(j).true_positives = sum([per_movie_count(1:j-1).true_positives]); 
            end
            
            if is_frame_matches
                per_movie_count(j).false_positive_fly_frame = [[matches_all.(behav_list{i})(bitand(this_movie_mask, false_positive_idxs)).fly]; [matches_all.(behav_list{i})(bitand(this_movie_mask, false_positive_idxs)).frame_num]]';
                per_movie_count(j).false_negative_fly_frame = [[matches_all.(behav_list{i})(bitand(this_movie_mask, false_negat_idxs)).fly]; [matches_all.(behav_list{i})(bitand(this_movie_mask, false_negat_idxs)).frame_num]]';
            else
                per_movie_count(j).false_positive_fly_frame = [];
                per_movie_count(j).false_negative_fly_frame = [];
            end
        end
        
        true_false_count_all.(behav_list{i}).per_movie_count = per_movie_count;
        
        all_flies = cellfun(@(m, f) sprintf('Movie %s Fly %d', m, f), ...
            {matches_all.(behav_list{i})(:).movie}, {matches_all.(behav_list{i})(:).fly}, ...
            'UniformOutput', false);
        all_flies_w_duplicate = repelem(all_flies, associated_annot_per_bout_match);
        false_positive_flies_w_duplicate = all_flies_w_duplicate(false_positive_idxs_w_duplicate); 
        false_negative_flies_w_duplicate = all_flies_w_duplicate(false_negat_idxs_w_duplicate); 
        all_flies_set = unique(all_flies);
        
        per_fly_count = struct('movie', '', 'fly', nan, 'false_positives', nan, 'false_negatives', nan, 'true_positives', nan, 'true_total', nan, ...
            'false_positive_fly_frame', [], 'false_negative_fly_frame', []);
        for j=1:length(all_flies_set)+1
            if j<length(all_flies_set)+1
                tokens = regexp(all_flies_set{j}, 'Movie ([\w-]+) Fly (\d+)', 'tokens');
                per_fly_count(j).movie = tokens{1}{1};
                per_fly_count(j).fly = str2double(tokens{1}{2});
                this_fly_mask = strcmp(all_flies, all_flies_set(j));
                per_fly_count(j).false_positives = nnz(strcmp(false_positive_flies_w_duplicate, all_flies_set{j}));
                per_fly_count(j).false_negatives = nnz(strcmp(false_negative_flies_w_duplicate, all_flies_set{j})); 
                per_fly_count(j).true_total = nnz(strcmp(all_flies_w_duplicate, all_flies_set{j})) - per_fly_count(j).false_positives;
                per_fly_count(j).true_positives = per_fly_count(j).true_total - per_fly_count(j).false_negatives; 
            else
                per_fly_count(j).movie = 'summary';
                per_fly_count(j).false_positives = sum([per_fly_count(1:j-1).false_positives]);
                per_fly_count(j).false_negatives = sum([per_fly_count(1:j-1).false_negatives]); 
                per_fly_count(j).true_total = sum([per_fly_count(1:j-1).true_total]);
                per_fly_count(j).true_positives = sum([per_fly_count(1:j-1).true_positives]); 
            end
            
            if is_frame_matches
                per_fly_count(j).false_positive_fly_frame = [[matches_all.(behav_list{i})(bitand(this_fly_mask, false_positive_idxs)).fly]; [matches_all.(behav_list{i})(bitand(this_fly_mask, false_positive_idxs)).frame_num]]';
                per_fly_count(j).false_negative_fly_frame = [[matches_all.(behav_list{i})(bitand(this_fly_mask, false_negat_idxs)).fly]; [matches_all.(behav_list{i})(bitand(this_fly_mask, false_negat_idxs)).frame_num]]';
            else
                per_fly_count(j).false_positive_fly_frame = [];
                per_fly_count(j).false_negative_fly_frame = [];
            end
        end
        
        true_false_count_all.(behav_list{i}).per_fly_count = per_fly_count;
    end
    
    if is_frame_matches
        prefix = 'frame';
    else
        prefix = 'bout';
    end
    save(strcat(prefix, '_true_false_count_', match_sel_str, '.mat'), 'true_false_count_all');
end