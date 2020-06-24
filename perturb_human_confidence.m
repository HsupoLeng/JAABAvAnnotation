bout_match_str_ori = 'bout_matches_ALL_All_thres0.1'; 
frame_match_str_ori = 'frame_matches_ALL_All_thres0.1';
bout_match_str_new = strcat(bout_match_str_ori, '_rand_persv_ratio');
s_temp = load(strcat(bout_match_str_ori, '.mat'));
bout_matches_all_ori = s_temp.bout_matches_all; 
load('common-params-annot-analysis.mat');

flymat_id = 'ALL_All_thres0.1';
behav_sel = [1, 2, 3];
override_jab_list = {};
is_NoRel = false;
plot_bar = false;
plot_violin = true;
plot_box = false;
img_format = {'png', 'eps'};
num_repeat = 50;
gen_new_randomization = false;

%% Run statistical significance test on original data
% w.r.t. the hypothesis that all categories come from the same distribution
% then run pair-wise multi-comparison Bofferroni test
[~, score_post_analysis_struct] = analyze_human_jaaba_annot_corr_v3(flymat_id, behav_sel, override_jab_list, is_NoRel, strcat(bout_match_str_ori, '.mat'), plot_bar, plot_violin, plot_box, img_format);
[p_val_arr_ori, pairwise_test_p_val_ori] = run_stat_signif_jaaba_scores(behav_list, score_post_analysis_struct, 'on');

score_post_analysis_struct_frame_wise = analyze_human_jaaba_annot_corr_frame_wise(flymat_id, behav_sel, override_jab_list, is_NoRel, strcat(frame_match_str_ori, '.mat'), plot_bar, plot_violin, img_format);
[p_val_arr_ori_frame_wise, pairwise_test_p_val_ori_frame_wise] = run_stat_signif_jaaba_scores(behav_list, score_post_analysis_struct_frame_wise, 'on');

%% Randomize human annotation score for each bout multiple times
recall_rate_all = zeros(num_repeat, 3, 6);
ori_human_conf_prob_cell = cell(size(behav_list));
p_val_all_rand = zeros(num_repeat, 3);
for i=1:length(behav_list)
    ori_human_conf = [bout_matches_all_ori.(behav_list{i})(:).annot_score];
    ori_human_conf_count = histcounts(ori_human_conf, 'BinMethod', 'integers', 'BinLimits', [1, 6]);
    ori_human_conf_count_accum = cumsum(ori_human_conf_count);
    ori_human_conf_prob_cell{i} = ori_human_conf_count_accum./ori_human_conf_count_accum(end);
end

for k=1:num_repeat
    if gen_new_randomization
        % Randomize human combined confidence score for each bout, preserving ratio
        % of combined scores
        bout_matches_all = bout_matches_all_ori; 
        for i=1:length(behav_list)
            ori_human_conf = [bout_matches_all_ori.(behav_list{i})(:).annot_score];
            rand_human_conf_src = rand(size(ori_human_conf));
            rand_human_conf_mask = arrayfun(@(prob) rand_human_conf_src < prob, ori_human_conf_prob_cell{i}, 'UniformOutput', false);
            rand_human_conf = 7 - sum(vertcat(rand_human_conf_mask{:}));
            rand_human_conf(~ori_human_conf) = 0;
            rand_human_conf_cell = mat2cell(rand_human_conf, 1, arrayfun(@(s) length(s.annot_score), bout_matches_all.(behav_list{i})));
            for j=1:length(rand_human_conf_cell)
                if ~bout_matches_all.(behav_list{i})(j).annot_score 
                    continue; 
                else
                    bout_matches_all.(behav_list{i})(j).annot_score = rand_human_conf_cell{j};
                end
            end
        end

        save(sprintf('%s_%d.mat', bout_match_str_new, k), 'bout_matches_all');
    end
        
    % Compute and store recall rate for each randomization
    [recall_rate_all(k, :, :), score_post_analysis_struct] = analyze_human_jaaba_annot_corr_v3(flymat_id, behav_sel, override_jab_list, is_NoRel, sprintf('%s_%d.mat', bout_match_str_new, k), plot_bar, plot_violin, plot_box, img_format);
    
    % Run statistical significance test w.r.t. null hypothesis that all
    % categories come from the same distribution
    [p_val_all_rand(k, :), ~] = run_stat_signif_jaaba_scores(behav_list, score_post_analysis_struct, 'off');
end

% Compute recall rate distribution over the randomizations
for i=1:length(behav_sel)
    figure();
    hold on;
    plot(1:6, squeeze(recall_rate_all(:, i, :)).*100, 'Color', [0.7, 0.7, 0.7]);
    errorbar(1:6, mean(squeeze(recall_rate_all(:, i, :)).*100, 'omitnan'), 2*std(squeeze(recall_rate_all(:, i, :)).*100, 'omitnan'), 'LineWidth', 2);
    xlim([0.5, 6.5]);
    title({'Recall rate per confidence category', sprintf('Behavior: %s Randomization: %d times', strrep(behav_list{i}, '_', '-'), num_repeat)});   
    xlabel('Human annotation combined score \in [1,6]');
    ylabel('recall in %');
    hold off;
    set(gcf,'renderer','Painters');
    saveas(gcf, sprintf('perturbed_recall_rate_curve-%s-random_%d_times.png', behav_list{i}, num_repeat));
    saveas(gcf, sprintf('perturbed_recall_rate_curve-%s-random_%d_times.eps', behav_list{i}, num_repeat), 'epsc');
end

%% Compute confidence interal for the means
recall_rate_all_mean = squeeze(mean(recall_rate_all));
recall_rate_all_std = squeeze(std(recall_rate_all));

recall_rate_all_dist_fit = squeeze(cellfun(@(x) fitdist(x, 'normal'), num2cell(recall_rate_all, 1), 'UniformOutput', false));
recall_rate_all_param_conf_intv = cellfun(@paramci, recall_rate_all_dist_fit, 'UniformOutput', false);
recall_rate_all_mean_conf_intv = cellfun(@(a) a(:, 1), recall_rate_all_param_conf_intv, 'UniformOutput', false);
recall_rate_all_std_conf_intv = cellfun(@(a) a(:, 2), recall_rate_all_param_conf_intv, 'UniformOutput', false);
%recall_rate_all_mean_conf_range = 2.021 .* recall_rate_all_std./sqrt(size(recall_rate_all, 1));
%recall_rate_all_mean_conf_intv = arrayfun(@(m, r) [m-r, m+r], recall_rate_all_mean, recall_rate_all_mean_conf_range, 'UniformOutput', false); 

%% Draw the plots
color_opts = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};
figure();
hold on;

for i=1:length(behav_sel)
    scatter(i+0.4.*rand(num_repeat, 1)-0.2, p_val_all_rand(:, i), [], color_opts{i});
    scatter(i, p_val_arr_ori(i), [], color_opts{i}, 'filled');
end
set(gca, 'YScale', 'log');
x_lims = xlim;
line([x_lims(1), x_lims(2)], [-2, -2], 'Linestyle', '--', 'Color', 'k');
hold off; 
set(gcf,'renderer','Painters');
% yticks(fliplr([log10(0.999), log10(0.99), log10(0.9), log10(0.5), -1, -2, -3, -10, -50, -100, -150]));
% yticklabels(fliplr([0.999, 0.99, 0.9, 0.5, 10.^[-1, -2, -3, -10, -50, -100, -150]]));
yticks(fliplr(10.^[-1, -2, -3, -10, -50, -100, -150]));
% yticklabels(fliplr(10.^[-1, -2, -3, -10, -50, -100, -150]));
ylabel('p value');
rand_p_signif_count = sum(p_val_all_rand < 0.01);
xticks(1:3);
xticklabels({sprintf('Lunge (N_{rand, p<0.01} = %d)', rand_p_signif_count(1)), ...
    sprintf('Wing extension (N_{rand, p<0.01} = %d)', rand_p_signif_count(2)), ...
    sprintf('Headbutt (N_{rand, p<0.01} = %d)', rand_p_signif_count(3))});
% set(gca, 'TickLabelInterpreter', 'tex');
legend({'Lunge randomized', 'Lunge original', 'Wing extension randomized', 'Wing extension original', 'Headbutt randomzied', 'Headbutt original'}, 'Location', 'northwest');
%% Save results
save('stat_signif_test_data.mat', 'p_val_all_rand', 'p_val_arr_ori', ...
    'pairwise_test_p_val_ori', 'rand_p_signif_count', 'recall_rate_all', ...
    'recall_rate_all_mean', 'recall_rate_all_std', ...
    'recall_rate_all_mean_conf_intv', 'recall_rate_all_std_conf_intv');
%% Test statistical significance that jaaba score come from different distributions (or not)
% Now using average jaaba score for each bout
function [p_val_arr, pairwise_test_p_val_struct] = run_stat_signif_jaaba_scores(behav_list, score_post_analysis_struct, disp)
    p_val_arr = zeros(size(score_post_analysis_struct));
    
    pairwise_test = [1, 2; 1, 3; 2, 3; 2, 4; 3, 4; 3, 5; 4, 5; 4, 6; 5, 6];
    pairwise_test_p_val_struct = struct(); 
    for i = 1:length(score_post_analysis_struct)
%         false_negat_idx = score_post_analysis_struct(i).false_negat_idxs; 
        annot_score_max = score_post_analysis_struct(i).annot_score;
        annot_score_p = annot_score_max(annot_score_max > 0);
        % randomly choose to make sample size equal in each category (to the category smallest in size)
        annot_score_unique = unique(annot_score_p);
        annot_score_count = arrayfun(@(s) sum(annot_score_p == s), annot_score_unique);
        sample_per_cat = min(annot_score_count);
        chosen_idx_all = zeros(sample_per_cat, length(annot_score_unique));
        for j=1:length(annot_score_unique)
            candidate_idx = find(annot_score_p == annot_score_unique(j));
            chosen_idx = candidate_idx(randi(length(candidate_idx), [sample_per_cat, 1]));
            chosen_idx_all(:, j) = chosen_idx;
        end

%         annot_score_chosen = annot_score_p(chosen_idx_all(:));
        jaaba_score = score_post_analysis_struct(i).jaaba_score;
        jaaba_score_p = jaaba_score(annot_score_max > 0);
%         jaaba_score_chosen = jaaba_score_max(chosen_idx_all(:));

        [p_val_arr(i), ~, stats] = kruskalwallis(jaaba_score_p, annot_score_p, disp);
%         pairwise_test_cell{i} = multcompare(stats,  'CType', 'bonferroni', 'Display', disp);
        for j=1:size(pairwise_test, 1)
            pairwise_test_p_val_struct.(behav_list{i})(j).categories = pairwise_test(j, :);
            pairwise_test_p_val_struct.(behav_list{i})(j).p_val = ranksum(...
                jaaba_score_p(annot_score_p == pairwise_test(j, 1)), ...
                jaaba_score_p(annot_score_p == pairwise_test(j, 2)));
        end
    end
end
