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
        
        annot_score = uint8(full(cellfun(@(scores) max(scores), {matches_all.(behav_list{i}).annot_score})));
        false_positive_idxs = cellfun(@(scores) max(scores), {matches_all.(behav_list{i}).annot_score}) == 0;
        if is_frame_matches
            false_negat_idxs = [matches_all.(behav_list{i}).is_false_negative];
        else
            false_negat_idxs = [matches_all.(behav_list{i}).virtual_jaaba_match];
        end
        
        human_annot_scores = 0:6; 
        virt_annot_score = annot_score(false_negat_idxs); 
        true_annot_score = annot_score(~false_negat_idxs);
        false_negat_accum = hist(virt_annot_score, human_annot_scores);
        positive_accum = hist(true_annot_score, human_annot_scores);
        true_counts = [0, positive_accum(2:end)];
        false_counts = [positive_accum(1), false_negat_accum(2:end)];
        per_category_struct = struct('score', nan, 'false_positives', nan, 'false_negatives', nan, 'true_positives', nan, 'true_total', nan);
        for j=1:length(false_counts)
            per_category_struct(j).score = j-1;
            if j==1
                per_category_struct(j).false_positives = false_counts(j);
                per_category_struct(j).false_negatives = 0; 
                per_category_struct(j).true_positives = 0;
                per_category_struct(j).true_total = 0;
            else
                per_category_struct(j).false_positives = 0; 
                per_category_struct(j).false_negatives = false_counts(j); 
                per_category_struct(j).true_positives = true_counts(j);
                per_category_struct(j).true_total = per_category_struct(j).false_negatives ...
                    + per_category_struct(j).true_positives;
            end
        end
        true_false_count_all.(behav_list{i}).per_category_count = per_category_struct;
        
        false_positive_sources = {matches_all.(behav_list{i})(false_positive_idxs).movie}; 
        false_negative_sources = {matches_all.(behav_list{i})(false_negat_idxs).movie}; 
        all_sources = {matches_all.(behav_list{i})(:).movie};
        all_sources_set = unique(all_sources);
        
        per_movie_count = struct('movie', '', 'false_positives', nan, 'false_negatives', nan, 'true_positives', nan, 'true_total', nan, ...
            'false_positive_fly_frame', [], 'false_negative_fly_frame', []);
        for j=1:length(all_sources_set)
            per_movie_count(j).movie = all_sources_set{j};
            this_movie_mask = strcmp(all_sources, all_sources_set(j));
            per_movie_count(j).false_positives = nnz(strcmp(false_positive_sources, all_sources_set{j}));
            per_movie_count(j).false_negatives = nnz(strcmp(false_negative_sources, all_sources_set{j})); 
            per_movie_count(j).true_total = nnz(strcmp(all_sources, all_sources_set{j})) - per_movie_count(j).false_positives;
            per_movie_count(j).true_positives = per_movie_count(j).true_total - per_movie_count(j).false_negatives; 
            
            if is_frame_matches
                per_movie_count(j).false_positive_fly_frame = [[matches_all.(behav_list{i})(bitand(this_movie_mask, false_positive_idxs)).fly]; [matches_all.(behav_list{i})(bitand(this_movie_mask, false_positive_idxs)).frame_num]]';
                per_movie_count(j).false_negative_fly_frame = [[matches_all.(behav_list{i})(bitand(this_movie_mask, false_negat_idxs)).fly]; [matches_all.(behav_list{i})(bitand(this_movie_mask, false_negat_idxs)).frame_num]]';
            else
                per_movie_count(j).false_positive_fly_frame = [];
                per_movie_count(j).false_negative_fly_frame = [];
            end
        end
        
        true_false_count_all.(behav_list{i}).per_movie_count = per_movie_count;
    end
    
    if is_frame_matches
        prefix = 'frame';
    else
        prefix = 'bout';
    end
    save(strcat(prefix, '_true_false_count_', match_sel_str, '.mat'), 'true_false_count_all');
end