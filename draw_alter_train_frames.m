% Compile and draw percision and recall for classifiers trained with
% downsampled training set
load('common-params-annot-analysis.mat', 'behav_list', 'jab_list');
data_mat_prefix = 'frame_true_false_count_All_All_thres0.1_percent';
frames_ratio_opts_template = [0.05, 0.25, 0.5, 0.65, 0.75];
actual_basal_pairing_frames_ratio = [0.69, 0.62, 0.53];
ref = load('D:\xubo\code\annot-analysis\true_false_count\perframe\frame_true_false_count_ALL_All_thres0.1.mat');
ref = ref.true_false_count_all;
jab_root = 'D:\xubo\classifiers';
num_repeat = 10; 

%%
prec_recall_all = zeros(num_repeat, length(behav_list), 2, length(frames_ratio_opts_template));
prec_recall_by_cat_all = zeros(num_repeat, length(behav_list), 2, length(frames_ratio_opts_template), 6);
for iter=1:num_repeat
    for i=1:length(frames_ratio_opts_template)
        if i == length(frames_ratio_opts_template) - 1
            bout_true_false_count_file = sprintf('%s_%s_repeat_%d.mat', data_mat_prefix, 'vary', iter);
        else
            bout_true_false_count_file = sprintf('%s_%d_repeat_%d.mat', data_mat_prefix, floor(frames_ratio_opts_template(i)*100), iter);
        end
        temp = load(bout_true_false_count_file);
        for j=1:length(behav_list)
            prec_recall_all(iter, j, 1, i) = temp.true_false_count_all.(behav_list{j}).precision * 100;
            prec_recall_all(iter, j, 2, i) = temp.true_false_count_all.(behav_list{j}).recall * 100;
            
            prec_recall_by_cat_all(iter, j, 2, i, :) = [temp.true_false_count_all.(behav_list{j}).per_category_count(2:end-1).recall]*100;
        end
        
%         actual_train_frames_ratio = struct();
%         for j=1:length(behav_list)
%            [~, jab_id, ~] = fileparts(jab_list{j}); 
%            
%            [ref_classifier_p_sz, ref_classifier_n_sz] = report_training_set_size(jab_list{j});
%            
%            requested_frame_ratio = frames_ratio_opts_template(i);
%            if requested_frame_ratio == 0.65
%                requested_frame_ratio = actual_basal_pairing_frames_ratio(j);
%            end
%            new_classifier_path = fullfile(jab_root, '\rand_train_subset', sprintf('%s_percent_%d_repeat_%d.jab', jab_id, floor(requested_frame_ratio*100), iter));
%            [new_classifier_p_sz, new_classifier_n_sz] = report_training_set_size(new_classifier_path);
%            actual_train_frames_ratio(j).behavior = behav_list{j};
%            actual_train_frames_ratio(j).positive_train_frame_ratio = new_classifier_p_sz/ref_classifier_p_sz; 
%            actual_train_frames_ratio(j).negative_train_frame_ratio = new_classifier_n_sz/ref_classifier_n_sz; 
%            actual_train_frames_ratio(j).total_train_frame_ratio = (new_classifier_p_sz + new_classifier_n_sz)/(ref_classifier_p_sz + ref_classifier_n_sz);
%         end
%         true_false_count_all = temp.true_false_count_all;
%         save(bout_true_false_count_file, 'true_false_count_all', 'actual_train_frames_ratio');
    end
end
%%
mean_prec_recall_w_conf_all = struct();
for j=1:length(behav_list)
    frames_ratio_opts = [frames_ratio_opts_template(1:end-2), actual_basal_pairing_frames_ratio(j), frames_ratio_opts_template(end)];
    [frames_ratio_opts, sorted_idx] = sort(frames_ratio_opts);
    mean_prec_recall_w_conf = struct();
    
    figure();
    hold on; 
    data_temp = squeeze(prec_recall_all(:, j, 1, :));
    data_temp = data_temp(:, sorted_idx);
    std_err = std(data_temp, 0)./sqrt(num_repeat);
    CI95 = tinv(0.975, num_repeat-1);                   
    yCI95 = std_err.*CI95(:); 
    prec_mean = mean(data_temp); 
    errorbar(frames_ratio_opts*100, prec_mean, yCI95, 'LineWidth', 3);
    
    for k=1:length(frames_ratio_opts)
        mean_prec_recall_w_conf(k).train_frame_ratio = frames_ratio_opts(k);
        mean_prec_recall_w_conf(k).mean_precision = prec_mean(k);
        mean_prec_recall_w_conf(k).mean_precision_conf_interval = prec_mean(k) + [-yCI95(k), yCI95(k)];
    end
    
    set(gca, 'ColorOrderIndex', 1);
    plot(frames_ratio_opts*100, repmat(ref.(behav_list{j}).precision*100, 1, length(frames_ratio_opts)), 'LineStyle', '--');
    ylim([0, 100])
    ylabel('Precision (%)');
    
    for iter=1:num_repeat
        plot(frames_ratio_opts*100, data_temp, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '-');
    end
    
    title(sprintf('Behavior %s', behav_list{j}), 'Interpreter', 'none');
    xticks(frames_ratio_opts*100);
    xticklabels(cellstr(num2str(floor(frames_ratio_opts'*100))));
    xlabel('Percentage of training set size relative to the original (%)');
    xlim([0, 80]);
    
    legend({'Precision', 'Original precision'}, 'Location', 'southeast');
    hold off; 
    set(gcf, 'renderer', 'Painters');
    saveas(gcf, sprintf('frame_precision_vs_training_set_size_%s.eps', behav_list{j}), 'epsc');
    
    figure();
    hold on;
    data_temp = squeeze(prec_recall_all(:, j, 2, :));
    data_temp = data_temp(:, sorted_idx);
    std_err = std(data_temp, 0)./sqrt(num_repeat);
    CI95 = tinv(0.975, num_repeat-1);                   
    yCI95 = std_err.*CI95(:); 
    mean_recall = mean(data_temp);
    errorbar(frames_ratio_opts*100, mean_recall, yCI95, 'LineWidth', 3);
    
    set(gca, 'ColorOrderIndex', 1);
    plot(frames_ratio_opts*100, repmat(ref.(behav_list{j}).recall*100, 1, length(frames_ratio_opts)), 'LineStyle', '--');
    ylim([0, 100]);
    ylabel('Recall (%)');
    
    for k=1:length(frames_ratio_opts)
        mean_prec_recall_w_conf(k).mean_recall = mean_recall(k);
        mean_prec_recall_w_conf(k).mean_recall_conf_interval = mean_recall(k) + [-yCI95(k), yCI95(k)];
    end
    
    for iter=1:num_repeat
        plot(frames_ratio_opts*100, data_temp, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '-');
    end
    
    title(sprintf('Behavior %s', behav_list{j}), 'Interpreter', 'none');
    xticks(frames_ratio_opts*100);
    xticklabels(cellstr(num2str(floor(frames_ratio_opts'*100))));
    xlabel('Percentage of training set size relative to the original (%)');
    xlim([0, 80]);
    
    legend({'Mean recall', 'Original recall'}, 'Location', 'southeast');
    hold off; 
    set(gcf, 'renderer', 'Painters');
    saveas(gcf, sprintf('frame_recall_vs_training_set_size_%s.eps', behav_list{j}), 'epsc');
    
    figure();
    hold on;
    data_temp = squeeze(prec_recall_by_cat_all(:, j, 2, :, :));
    data_temp = data_temp(:, sorted_idx, :);
    std_err = std(data_temp, 0)./sqrt(num_repeat);
    CI95 = tinv(0.975, num_repeat-1);                   
    yCI95 = squeeze(std_err.*CI95(:)); 
    mean_recall = squeeze(mean(data_temp));
    
    for k=1:length(frames_ratio_opts)
        for l=1:size(mean_recall, 2)
            mean_prec_recall_w_conf(k).(sprintf('score_%d_mean_recall', l)) = mean_recall(k, l);
            mean_prec_recall_w_conf(k).(sprintf('score_%d_mean_recall_conf_interval', l)) = mean_recall(k, l) + [-yCI95(k, l), yCI95(k, l)];
        end
    end
    
    for k=1:size(mean_recall, 1)
        errorbar(1:6, mean_recall(k, :), yCI95(k, :), 'LineWidth', 3);
    end
    plot(1:6, [ref.(behav_list{j}).per_category_count(2:end-1).recall].*100, 'LineWidth', 3, 'LineStyle', '--');
    legend([arrayfun(@(r) sprintf('Training set size %d%% w. 95%% confidence interval (N=10)', floor(r*100)), frames_ratio_opts, 'UniformOutput', false), 'Original recall'], 'Location', 'southeast');
    xticks(1:6);
    xlabel('Human annotation combined score \in [0,6]');
    ylim([0, 100]);
    ylabel('Recall (%)');
    title({'Recall vs. human confidence for classifiers trained on downsized training set', sprintf('Behavior: %s', behav_list{j})}, 'Interpreter', 'none');
    hold off; 
    set(gcf, 'renderer', 'Painters');
    saveas(gcf, sprintf('frame_recall_vs_human_confidence_%s.eps', behav_list{j}), 'epsc');
    
    mean_prec_recall_w_conf_all(j).behavior = behav_list{j};
    mean_prec_recall_w_conf_all(j).data = mean_prec_recall_w_conf; 
end

%% Convert confidence interval to Kenta's EXCEL format
conf_intv_all = struct();
for j=1:length(behav_list)
    data_struct_temp = mean_prec_recall_w_conf_all(j).data;
    data_temp = struct2table(data_struct_temp);
    conf_intv_one_behav = struct();
    for i=1:size(data_temp, 1)
        conf_intv_temp = data_temp(i, (6+2*6-1):-2:6);
        conf_intv_temp = table2cell(conf_intv_temp);
        conf_intv_temp = vertcat(conf_intv_temp{:})./100;
        conf_intv_temp = conf_intv_temp(:, [2, 1]);
        conf_intv_one_behav.(sprintf('percentage_%d', floor(100*data_temp.train_frame_ratio(i)))) = ...
            conf_intv_temp;
        
        conf_intv_one_behav.(sprintf('precision_percentage_%d',  floor(100*data_temp.train_frame_ratio(i)))) = ...
            data_struct_temp(i).mean_precision_conf_interval([2, 1])./100;
    end
    conf_intv_all.(behav_list{j}) = conf_intv_one_behav;
end 