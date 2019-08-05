% This function is written for sanity check purposes. It compiles all false
% positives in bout_matches_all, so that a comparison can be made with the
% annotated false positives in FLYMAT_NewTrainingFiles-OriginalCopy.mat
function false_positives_all = compile_false_positives(thres)
    load(strcat('bout_matches_all-thres', num2str(thres), '.mat'), 'bout_matches_all');
    load('common-params-annot-analysis.mat', ...
        'annot_file', 'behav_list', 'behav_shorthands', 'jab_list', 'movie_name_list');

    false_positives_all_init_args = [behav_list; cell(1, length(behav_list))];
    false_positives_all = struct(false_positives_all_init_args{:});

    for i=1:length(behav_list)
        bout_matches = bout_matches_all.(behav_list{i});
        false_positive_idxs = [bout_matches(:).annot_score] == 0;
        fields_to_extract = {'movie', 'fly', 'jaaba_bout_start', 'jaaba_bout_end'};
        struct_args = cell(2, length(fields_to_extract));
        for j=1:length(fields_to_extract)
            struct_args{1, j} = fields_to_extract{j};
            struct_args{2, j} = {bout_matches(false_positive_idxs).(struct_args{1, j})}';
        end
        false_positives = struct(struct_args{:});
        false_positives_all.(behav_list{i}) = false_positives;
    end
end
