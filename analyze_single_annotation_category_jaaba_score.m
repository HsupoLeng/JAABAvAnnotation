load('bout_matches_all-thres0.1.mat');
load('FLYMAT_HumanAnnotation_new.mat');

behaviour_to_analyze = 'Lunges_New';
behaviour_shorthand = 'L';
score_to_analyze = 3;
score_to_compare = 4; 

bout_mask = [bout_matches_all.(behaviour_to_analyze).annot_score] == score_to_analyze;
bouts_to_analyze = bout_matches_all.(behaviour_to_analyze)(bout_mask);
bout_mask = [bout_matches_all.(behaviour_to_analyze).annot_score] == score_to_compare;
bouts_to_compare_1 = bout_matches_all.(behaviour_to_analyze)(bout_mask);
bout_mask = [bout_matches_all.(behaviour_to_analyze).annot_score] == 5;
bouts_to_compare_2 = bout_matches_all.(behaviour_to_analyze)(bout_mask);

curr_movie_and_fly = cell(1, 2);
for i=1:length(bouts_to_analyze)
    if ~isequal(curr_movie_and_fly, {bouts_to_analyze(i).movie, bouts_to_analyze(i).fly})
        curr_movie_and_fly = {bouts_to_analyze(i).movie, bouts_to_analyze(i).fly};
        bouts_to_analyze(i).bout_cnt = 1;
    else
        bouts_to_analyze(i).bout_cnt = bouts_to_analyze(i-1).bout_cnt + 1;
    end
end

for i=1:length(bouts_to_analyze)
    bout = bouts_to_analyze(i);
    fly_idx = find_fly_in_flymat(flymatHumanAnnot, bout.movie, bout.fly, false);
    annot_score_pair = flymatHumanAnnot(fly_idx).(strcat(behaviour_shorthand, '_scores'))(bout.bout_cnt, :);
    if any(isnan(annot_score_pair))
        bouts_to_analyze(i).is_single_labeled = true;
    else
        bouts_to_analyze(i).is_single_labeled = false; 
    end
end

bouts_single_labeled = bouts_to_analyze([bouts_to_analyze.is_single_labeled]);
bouts_both_labeled = bouts_to_analyze(~[bouts_to_analyze.is_single_labeled]);
figure();
hold on;
boxplot([[bouts_single_labeled.jaaba_score_avg_normed]'; [bouts_both_labeled.jaaba_score_avg_normed]'; ...
    [bouts_to_compare_1.jaaba_score_avg_normed]'; [bouts_to_compare_2.jaaba_score_avg_normed]'], ...
    [ones(length(bouts_single_labeled), 1); 2.*ones(length(bouts_both_labeled), 1); ...
    3.*ones(length(bouts_to_compare_1), 1); 4.*ones(length(bouts_to_compare_2), 1)]);
xlim([0.5, 4.5]);
xticks(1:4);
xticklabels({'Bouts w. score 3 labeled by one annotator', ...
    'Bouts w. score 3 labeled by two annotators', ...
    'Bouts w. score 4', 'Bouts w. score 5'});
hold off; 