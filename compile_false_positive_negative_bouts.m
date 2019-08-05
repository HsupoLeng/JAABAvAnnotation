clear; 
bout_match_dir = 'D:\xubo\code\annot-analysis\bout_frame_matches';
bout_frame_match_filelist = dir(bout_match_dir);
bout_frame_match_filelist = {bout_frame_match_filelist(:).name};
bout_match_filelist = bout_frame_match_filelist(contains(bout_frame_match_filelist, 'bout_matches'));

for i=1:length(bout_match_filelist)
    load(bout_match_filelist{i}, 'bout_matches_all');
    fields = fieldnames(bout_matches_all);
    
    false_positive_bouts_all = struct();
    false_negative_bouts_all = struct();
    for j=1:length(fields)
        bout_matches = bout_matches_all.(fields{j});
        if ~isempty(bout_matches)
            false_positive_mask = ~arrayfun(@(s) max(s.annot_score), bout_matches);
            false_negative_mask = [bout_matches(:).virtual_jaaba_match];
            false_positive_bouts_all.(fields{j}) = bout_matches(false_positive_mask);
            false_negative_bouts_all.(fields{j}) = bout_matches(false_negative_mask);

            % Remove fields that are unneccessary
            false_positive_bouts_all.(fields{j}) = rmfield(false_positive_bouts_all.(fields{j}), {'annot_union_start', 'annot_union_end', ...
                'annot_intsct_start', 'annot_intsct_end', 'annot_score', 'multi_match', 'virtual_jaaba_match', 'invalid_match'});
            false_negative_bouts_all.(fields{j}) = rmfield(false_negative_bouts_all.(fields{j}), ...
                {'jaaba_bout_start', 'jaaba_bout_end', 'jaaba_score_avg', 'jaaba_score_min', 'jaaba_score_max', ...
                'jaaba_score_avg_normed', 'jaaba_score_min_normed', 'jaaba_score_max_normed', 'multi_match', 'virtual_jaaba_match', 'invalid_match'});
        else
            false_positive_bouts_all.(fields{j}) = [];
            false_negative_bouts_all.(fields{j}) = [];
        end
    end
    
    save(strcat('false_positive_', bout_match_filelist{i}), 'false_positive_bouts_all');
    save(strcat('false_negative_', bout_match_filelist{i}), 'false_negative_bouts_all');
end


