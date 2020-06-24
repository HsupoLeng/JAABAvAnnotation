%% Load model and data
data_root = 'D:\xubo\code\unsupervised-behavior-sequence-identification\';
analysis_root = 'D:\xubo\NewTrainingFiles-OriginalCopy\movies\gmm_analysis';
gt_root = 'D:\xubo\code\annot-analysis\bout_frame_matches';
load(fullfile(data_root, 'fitted_gmm-Mar_26.mat'), 'feature_names');
load(fullfile(gt_root, 'bout_matches_ALL_All_thres0.1.mat'));
load('common-params-annot-analysis.mat', 'behav_list');

load(fullfile(analysis_root, 'markov_transition-sort_by_feature-python_gmm_Mar_27_a-ana_Mar_28_214536-w_256_clusters-catch_bag_false.mat'), ...
    'movie_fly_sels', 'movie_list');

%% Prepare feature data
load(fullfile(data_root, 'feature_mat-all_data-Mar_26.mat'), 'feature_mat_cell');
feature_mat_cell = cellfun(@(f) squeeze(num2cell(f, [1, 2]))', feature_mat_cell, 'UniformOutput', false);
feature_mat_cell = vertcat(feature_mat_cell{:})';


%% Compute mean and variance for GT vs non-GT frames
% Compile features in GT frames and those in non-GT frames
feature_mat_gt_cell = cell(size(movie_fly_sels))';
feature_mat_non_gt_cell = cell(size(feature_mat_gt_cell));
feature_ranksum_p_val = struct();
feature_quantiles_gt = struct();
feature_quantiles_non_gt = struct();
feature_rule_weight_sum = struct();
feature_outliers_ratio = struct();
rule_mat_name = {'LungeNewRetrain2017-0219-Top_98_rules.mat', 'WingExtNew2017-0112b-Top_95_rules.mat', ...
    'HeadbuttNewRetrain2018-0905-Top_97_rules.mat'};

for k=1:length(behav_list)
    for i=1:size(movie_fly_sels, 1)
        movie_idx = movie_fly_sels{i, 1};
        flies_in_chamber = movie_fly_sels{i, 2};
        for j=1:length(flies_in_chamber)
            movie_name_str = movie_list{movie_fly_sels{i, 1}};
            movie_name_elems = strsplit(movie_name_str, '\');
            movie_name = movie_name_elems{1};
            fly_abs_idx = flies_in_chamber(j);

            bout_matches = bout_matches_all.(behav_list{k});
            bout_matches = bout_matches(bitand(contains({bout_matches(:).movie}, movie_name), [bout_matches(:).fly] == fly_abs_idx));
            bout_matches = bout_matches(~isnan([bout_matches(:).annot_union_start]));
            if ~isempty(bout_matches)
                gt_frames = arrayfun(@(s) s.annot_union_start:s.annot_union_end-1, bout_matches, 'UniformOutput', false);
                gt_frames = horzcat(gt_frames{:});
                gt_mask = boolean(accumarray(gt_frames', 1, [size(feature_mat_cell{j, i}, 1), 1]));
                feature_mat_gt_cell{j, i} = feature_mat_cell{j, i}(gt_mask, :);
                feature_mat_non_gt_cell{j, i} = feature_mat_cell{j, i}(~gt_mask, :);
            else
                feature_mat_gt_cell{j, i} = [];
                feature_mat_non_gt_cell{j, i} = feature_mat_cell{j, i};
            end
        end
    end

    feature_mat_gt = vertcat(feature_mat_gt_cell{:});
    feature_mat_non_gt = vertcat(feature_mat_non_gt_cell{:});
    
    ranksum_p_val_temp = nan(length(feature_names), 1);
    for m=1:size(feature_mat_gt, 2)
        ranksum_p_val_temp(m) = ranksum(feature_mat_gt(:, m), feature_mat_non_gt(:, m));
    end
    feature_ranksum_p_val.(behav_list{k}) = ranksum_p_val_temp;
    
    feature_outliers_ratio_temp = nan(length(feature_names), 2);
    for m=1:size(feature_mat_gt, 2)
        gt_qt = quantile(feature_mat_gt(:, m), [0.25, 0.75]);
        gt_q1 = gt_qt(1); gt_q3 = gt_qt(2);
        feature_outliers_ratio_temp(m, 1) = sum(bitor(feature_mat_gt(:, m) < gt_q1 - 1.5*(gt_q3-gt_q1), ...
            feature_mat_gt(:, m) > gt_q3 + 1.5*(gt_q3-gt_q1)))/size(feature_mat_gt, 1);
        non_gt_qt = quantile(feature_mat_non_gt(:, m), [0.25, 0.75]);
        non_gt_q1 = non_gt_qt(1); non_gt_q3 = non_gt_qt(2);
        feature_outliers_ratio_temp(m, 2) = sum(bitor(feature_mat_non_gt(:, m) < non_gt_q1 - 1.5*(non_gt_q3-non_gt_q1), ...
            feature_mat_non_gt(:, m) > non_gt_q3 + 1.5*(non_gt_q3-non_gt_q1)))/size(feature_mat_non_gt, 1);
    end
    feature_outliers_ratio.(behav_list{k}) = array2table(feature_outliers_ratio_temp, 'VariableNames', {'Behavior', 'Non_behavior'}, ...
        'RowNames', feature_names);
    
    figure((k-1)*2 + 2);
    feature_mat_gt = [feature_mat_gt; nan(size(feature_mat_non_gt, 1) - size(feature_mat_gt, 1), ...
        size(feature_mat_gt, 2))];
    b2 = iosr.statistics.boxPlot(cat(3, feature_mat_gt, feature_mat_non_gt), 'limit', [0.5, 99.5], 'showOutliers', false, ...
        'medianColor', {'k', 'w'}, 'boxColor', {'m', 'b'}, ...
        'groupLabels', {sprintf('Groundtruth frames (N=%d)', sum(~isnan(feature_mat_gt(:, 1)))), ...
        sprintf('Non-groundtruth frames (N=%d)', sum(~isnan(feature_mat_non_gt(:, 1))))}, ...
        'showLegend', true);
    xticks(1:size(feature_mat_gt, 2));
    xticklabels(strrep(feature_names, '_', '-'));
    xtickangle(45);
    ylabel('Feature z-score');
    set(gcf,'renderer','Painters');
    saveas(gcf, sprintf('feature_z_score-gt_vs_non_gt-%s-05_995_percentile.eps', behav_list{k}), 'epsc');
    saveas(gcf, sprintf('feature_z_score-gt_vs_non_gt-%s-05_995_percentile.png', behav_list{k})); 
    
    feature_quantiles_gt.(behav_list{k}) = quantile(feature_mat_gt, [0.75, 0.5, 0.25]);
    feature_quantiles_non_gt.(behav_list{k}) = quantile(feature_mat_non_gt, [0.75, 0.5, 0.25]);
    
    feature_var_gt = var(feature_mat_gt, 'omitnan');
    feature_var_non_gt = var(feature_mat_non_gt, 'omitnan');
    
%     figure();
%     hold on;
%     stem(2:2:2*size(feature_mat_non_gt, 2), feature_var_non_gt, 'filled');
%     stem(2:2:2*size(feature_mat_gt, 2), feature_var_gt, 'filled');
%     hold off; 
%     xticks(2:2:2*size(feature_mat_gt, 2));
%     xticklabels(strrep(feature_names, '_', '-'));
%     xtickangle(45);
%     legend({sprintf('Non-groundtruth frames (N=%d)', size(feature_mat_non_gt, 1)), ...
%         sprintf('Groundtruth frames (N=%d)', size(feature_mat_gt, 1))});
%     ylabel('Feature variance');
%     set(gcf,'renderer','Painters');
%     saveas(gcf, sprintf('feature_variance-gt_vs_non_gt-%s.eps', behav_list{k}), 'epsc');
    
    load(rule_mat_name{k});
    feature_rule_weight_sum_temp = nan(size(feature_names));
    for m=1:length(feature_names)
        feature_rule_weight_sum_temp(m) = sum(...
            [rules_struct_arr.rules(contains(...
            {rules_struct_arr.rules(:).feat_name}, feature_names{m})).rule_weight]);
    end
    feature_rule_weight_sum.(behav_list{k}) = feature_rule_weight_sum_temp;
    
%     figure();
%     bar(2:2:2*size(feature_mat_gt, 2), feature_rule_weight_sum_temp);
%     xticks(2:2:2*size(feature_mat_gt, 2));
%     xticklabels(strrep(feature_names, '_', '-'));
%     xtickangle(45);
%     ylabel('Summed feature rule weight');
%     set(gcf,'renderer','Painters');
%     saveas(gcf, sprintf('summed_feature_rule_weight-%s.eps', behav_list{k}), 'epsc');
end
