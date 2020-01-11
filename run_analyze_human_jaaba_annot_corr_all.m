load('common-params-annot-analysis.mat', 'all_jaaba_folder_list', 'jab_list');
% Initialize
plot_bar_bout_wise = false; 
plot_violin_bout_wise = false;
plot_bar_frame_wise = true; 
plot_violin_frame_wise = true; 

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

% ===== Top rules removed, re-trained classifiers =====
flymat_id_list = {{'ALL_All_thres0.1_minusTop3'}, {'ALL_All_thres0.1_minusTop5'}, ...
    {'ALL_All_thres0.1_minusTop10'}};
threshold_list = [0.1, 0.1, 0.1];

jab_file_list = repmat(jab_list, 1, 3); 
jab_file_list = reshape(jab_file_list, [], 3)';
jab_file_list = mat2cell(jab_file_list, [1, 1, 1], 3);
minus_rules = [3, 5, 10];
for i=1:length(jab_file_list)
    for j=1:length(jab_list)
        [~, filename, ~] = fileparts(jab_file_list{i}{j});
        jab_file_list{i}{j} = sprintf('%s_minusTop%dNew.jab', filename, minus_rules(i));
    end
end

maxgaps = struct('L', 1, 'WE', 3, 'HB', 1);
minbouts = struct('L', 3, 'WE', 6, 'HB', 3);
% minbouts = struct('L', 1, 'WE', 1, 'HB', 1); % Test min bout constraint
for i=1:length(flymat_id_list)
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
%     OrgData020816_XuboCopy('D:\xubo\NewTrainingFiles-OriginalCopy', ...
%         'TrainingSamples_genotype3.xlsx', thres, is_NoRel, maxgaps, minbouts, flymat_id_list{i});
    
    if length(flymat_id_list{i}) == 1
        flymat_id = flymat_id_list{i}{1};
        behav_sel = [1,2,3];
        
%         % Generate bout matches mat and bout-based bar and violin plots 
%         bout_match_str = strcat('bout_matches_', flymat_id, '.mat');
%         analyze_human_jaaba_annot_corr_v3(flymat_id, behav_sel, override_jab_list, is_NoRel, '', plot_bar_bout_wise, plot_violin_bout_wise, false, {'png', 'eps'});
%         % Generate bout-based false positive, false negative and true positive count 
%         compile_fp_fn_numbers(behav_sel, flymat_id, false);
        
        % For wing extension, generate corresponding frame-based plots
        frame_match_str = strcat('frame_matches_', flymat_id, '.mat');
        if plot_bar_frame_wise
            analyze_human_jaaba_annot_corr_frame_wise(flymat_id, [1, 2, 3], override_jab_list, is_NoRel, '', plot_bar_frame_wise, plot_violin_frame_wise, {'png', 'eps'});
            compile_fp_fn_numbers([1,2,3], flymat_id, true);
        end
    else    
        for j=1:length(flymat_id_list{i})
            flymat_id = flymat_id_list{i}{j};
            
%             bout_match_str = strcat('bout_matches_', flymat_id, '.mat');
%             analyze_human_jaaba_annot_corr_v3(flymat_id, j, override_jab_list, is_NoRel, '', plot_bar_bout_wise, plot_violin_bout_wise, false, {'png', 'eps'});
%             compile_fp_fn_numbers(j, flymat_id, false);
            if j == 2 && plot_bar_frame_wise
                % For wing extension, generate corresponding frame-based plots
                frame_match_str = strcat('frame_matches_', flymat_id, '.mat');
                analyze_human_jaaba_annot_corr_frame_wise(flymat_id, j, override_jab_list, is_NoRel, '', plot_bar_frame_wise, plot_violin_frame_wise, {'png', 'eps'});
                compile_fp_fn_numbers(j, flymat_id, true);
            end
        end
    end
end