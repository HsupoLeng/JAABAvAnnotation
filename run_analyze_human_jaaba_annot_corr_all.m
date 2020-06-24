% Run analysis of correlation between human annotation and JAABA predition
% for different sets of analyses
load('common-params-annot-analysis.mat', 'all_jaaba_folder_list', 'jab_list');
% Initialize
plot_bar_bout_wise = true; 
plot_violin_bout_wise = false;
plot_bar_frame_wise = true; 
plot_violin_frame_wise = false; 

% ===== Original classifiers =====
% flymat_id_list = {{'ALL_All_thresn0.1'}, {'ALL_All_thres0'}, {'ALL_All_thres0.1'}, {'ALL_All_thres0.2'}, ...
%     {'All_All_thres0.3'}...
%     {'ALL_All_woRelFeat'}, ...
%     {'L2017-0219_MvsMonly2_thres0.1', 'WE2017-0112b_MvsFonly_thres0.1', ...
%     'HB2018-0905_FvsFonly_thres0.1'}, ...
%     {'L2017-0219_MvsMandMvsF_thres0.1', 'WE2017-0112b_MvsFandFvsF_thres0.1', ...
%     'HB2018-0905_FvsFandMvsF_thres0.1'}, ...
%     {'L2017-0219_MvsMandFvsF_thres0.1', 'WE2017-0112b_MvsFandMvsM_thres0.1', ...
%     'HB2018-0905_FvsFandMvsM_thres0.1', }...
%     {'L2017-0219_MvsMandFvsFandMvsF_thres0.1', 'WE2017-0112b_MvsMandMvsFandFvsF_thres0.1', ...
%     'HB2018-0905_FvsFandMvsFandMvsM_thres0.1'}
%     };
% threshold_list = [-0.1, 0, 0.1, 0.2, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1];
% jab_file_list = {jab_list, jab_list, jab_list, jab_list, jab_list, ...
%     {'LungeNewRetrain2017-0219-woRelFeat_v2.jab', 'WingExtNew2017-0112b_woRelFeat_v2.jab', ...
%     'HeadbuttNewRetrain2018-0905-woRelFeat.jab'}, ...
%     {'LungeNewRetrain2017-0219_MvsMonly.jab', 'WingExtNew2017-0112b_MvsFonly.jab', ...
%     'HeadbuttNewRetrain2018-0905_FvsFonly.jab'}, ...
%     {'LungeNewRetrain2017-0219_MvsMandMvsF2.jab', 'WingExtNew2017-0112b_MvsFandFvsF.jab', ...
%     'HeadbuttNewRetrain2018-0905_FvsFandMvsF.jab'}, ...
%     {'LungeNewRetrain2017-0219_MvsMandFvsF.jab', 'WingExtNew2017-0112b_MvsFandMvsM.jab', ...
%     'HeadbuttNewRetrain2018-0905_FvsFandMvsM.jab'}, ...
%     {'LungeNewRetrain2017-0219_MvsMandFvsFandMvsF.jab', 'WingExtNew2017-0112b_MvsMandMvsFandFvsF.jab', ...
%     'HeadbuttNewRetrain2018-0905_FvsFandMvsMandMvsF.jab'}};

% % ===== Top rules removed, re-trained classifiers =====
% flymat_id_list = {{'ALL_All_thres0.1_minusTop3'}, {'ALL_All_thres0.1_minusTop5'}, ...
%     {'ALL_All_thres0.1_minusTop10'}};
% threshold_list = [0.1, 0.1, 0.1];
% 
% jab_file_list = repmat(jab_list, 1, 3); 
% jab_file_list = reshape(jab_file_list, [], 3)';
% jab_file_list = mat2cell(jab_file_list, [1, 1, 1], 3);
% minus_rules = [3, 5, 10];
% for i=1:length(jab_file_list)
%     for j=1:length(jab_list)
%         [~, filename, ~] = fileparts(jab_file_list{i}{j});
%         jab_file_list{i}{j} = sprintf('%s_minusTop%dNew.jab', filename, minus_rules(i));
%     end
% end

% % ===== Classifiers trained using ratio-preserved partial training set
% % frames_ratio_opts = [0.05, 0.25, 0.5, 0.75];
% frames_ratio_opts = [0.69, 0.62, 0.53];
% jab_root = 'D:\xubo\classifiers\rand_train_subset';
% flymat_id_common = 'All_All_thres0.1';
% % flymat_id_list = arrayfun(@(r) {sprintf('%s_percent_%d', flymat_id_common, floor(r*100))}, frames_ratio_opts, 'UniformOutput', false);
% threshold_list = [0.1, 0.1, 0.1];
% num_repeat = 10; 
% use_existing_match = false; 
% output_dir = "D:\xubo\code\annot-analysis\rand_train_subset_repeat";

% ===== Top features removed, re-trained classifiers =====
flymat_id_list = {{'ALL_All_thres0.1_woTop1Feat'}, {'ALL_All_thres0.1_woTop3Feat'}, ...
    {'ALL_All_thres0.1_woBottom3Feat'}};
threshold_list = [0.1, 0.1, 0.1];

jab_root = 'D:\xubo\classifiers\classifiers_without_weighted_features';
jab_file_list = repmat(jab_list, 1, 3); 
jab_file_list = reshape(jab_file_list, [], 3)';
jab_file_list = mat2cell(jab_file_list, [1, 1, 1], 3);
minus_rules = [1, 3, 3];
for i=1:length(jab_file_list)
    for j=1:length(jab_list)
        [~, filename, ~] = fileparts(jab_file_list{i}{j});
        if i <= 2
            jab_file_list{i}{j} = fullfile(jab_root, sprintf('%s_woTop%dFeat.jab', filename, minus_rules(i)));
        else
            jab_file_list{i}{j} = fullfile(jab_root, sprintf('%s_woBottom%dFeat.jab', filename, minus_rules(i)));
        end
    end
end
num_repeat = 1; 
use_existing_match = false; 
output_dir = "D:\xubo\code\annot-analysis\remove_top_features";

% % ===== Only-top-features, re-trained classifiers =====
% flymat_id_list = {{'ALL_All_thres0.1_Top3FeatOnly'}};
% threshold_list = [0.1];
% 
% jab_root = 'D:\xubo\classifiers\classifiers_with_only_weighted_features';
% jab_file_list = {jab_list};
% for i=1:length(jab_file_list)
%     for j=1:length(jab_list)
%         [~, filename, ~] = fileparts(jab_file_list{i}{j});
%         jab_file_list{i}{j} = fullfile(jab_root, sprintf('%s_Top3FeatOnly.jab', filename));
%     end
% end
% num_repeat = 1; 
% use_existing_match = false; 
% output_dir = "D:\xubo\classifiers\classifiers_with_only_weighted_features";


for iter=1:num_repeat
    % For random repeat analysis
% %     flymat_id_list = arrayfun(@(r) {sprintf('%s_percent_%d_repeat_%d', flymat_id_common, floor(r*100), iter)}, frames_ratio_opts, 'UniformOutput', false);
% %     jab_file_list = repmat(jab_list, 1, length(frames_ratio_opts)); 
% %     jab_file_list = reshape(jab_file_list, [], length(frames_ratio_opts))';
% %     jab_file_list = mat2cell(jab_file_list, ones(size(frames_ratio_opts)), 3);
%     flymat_id_list = {{sprintf('%s_percent_%s_repeat_%d', flymat_id_common, 'vary', iter)}};
%     jab_file_list = {jab_list};
%     for i=1:length(jab_file_list)
%         for j=1:length(jab_list)
%             [~, filename, ~] = fileparts(jab_file_list{i}{j});
% %             jab_file_list{i}{j} = fullfile(jab_root, sprintf('%s_percent_%d_repeat_%d.jab', filename,  floor(frames_ratio_opts(i)*100), iter));
%             jab_file_list{i}{j} = fullfile(jab_root, sprintf('%s_percent_%d_repeat_%d.jab', filename,  floor(frames_ratio_opts(j)*100), iter));
%         end
%     end

    maxgaps = struct('L', 1, 'WE', 3, 'HB', 1);
    minbouts = struct('L', 3, 'WE', 6, 'HB', 3);
    % minbouts = struct('L', 1, 'WE', 1, 'HB', 1); % Test min bout constraint
    for i=1:length(flymat_id_list)
        if i ~= 2 
            continue; 
        end
        % Predict using the classifiers
        JAABADetect(all_jaaba_folder_list, 'jabfiles', jab_file_list{i});

        thres = threshold_list(i); 
        override_jab_list = struct('behav_sel', {1;2;3}, 'jab_name', jab_file_list{i}');
        if contains(flymat_id_list{i}{1}, 'woRelFeat')
            is_NoRel = true; 
        else
            is_NoRel = false; 
        end

        % Run orgData
        OrgData020816_XuboCopy('D:\xubo\NewTrainingFiles-OriginalCopy', ...
            'TrainingSamples_genotype3.xlsx', thres, is_NoRel, maxgaps, minbouts, flymat_id_list{i}, output_dir);
            
        if length(flymat_id_list{i}) == 1
            flymat_id = flymat_id_list{i}{1};
            behav_sel = [1,2,3];
            
            if use_existing_match
                bout_match_str = strcat('bout_matches_', flymat_id, '.mat');
                frame_match_str = strcat('frame_matches_', flymat_id, '.mat');
            else
                bout_match_str = '';
                frame_match_str = '';
            end

            % Generate bout matches mat and bout-based bar and violin plots
            analyze_human_jaaba_annot_corr_v3(flymat_id, behav_sel, override_jab_list, is_NoRel, bout_match_str, plot_bar_bout_wise, plot_violin_bout_wise, false, {'png', 'eps'});
            % Generate bout-based false positive, false negative and true positive count 
            compile_fp_fn_numbers(behav_sel, flymat_id, false);

            % For wing extension, generate corresponding frame-based plots
            analyze_human_jaaba_annot_corr_frame_wise(flymat_id, [1, 2, 3], override_jab_list, is_NoRel, frame_match_str, plot_bar_frame_wise, plot_violin_frame_wise, {'png', 'eps'});
            compile_fp_fn_numbers([1,2,3], flymat_id, true);
        else    
            for j=1:length(flymat_id_list{i})
                flymat_id = flymat_id_list{i}{j};
               
                if use_existing_match
                    bout_match_str = strcat('bout_matches_', flymat_id, '.mat');
                    frame_match_str = strcat('frame_matches_', flymat_id, '.mat');
                else
                    bout_match_str = '';
                    frame_match_str = '';
                end

                analyze_human_jaaba_annot_corr_v3(flymat_id, j, override_jab_list, is_NoRel, bout_match_str, plot_bar_bout_wise, plot_violin_bout_wise, false, {'png', 'eps'});
                compile_fp_fn_numbers(j, flymat_id, false);
                if j == 2 
                    % For wing extension, generate corresponding frame-based plots
                    analyze_human_jaaba_annot_corr_frame_wise(flymat_id, j, override_jab_list, is_NoRel, frame_match_str, plot_bar_frame_wise, plot_violin_frame_wise, {'png', 'eps'});
                    compile_fp_fn_numbers(j, flymat_id, true);
                end
            end
        end
    end
end