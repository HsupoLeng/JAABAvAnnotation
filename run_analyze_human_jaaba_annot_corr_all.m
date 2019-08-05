load('common-params-annot-analysis.mat', 'all_jaaba_folder_list', 'jab_list');
% Initialize
plot_bar_bout_wise = false; 
plot_violin_bout_wise = false;
plot_bar_frame_wise = false; 
plot_violin_frame_wise = false; 
flymat_id_list = {{'ALL_All_thres0.1'}, {'ALL_All_thres0.2'}, {'ALL_All_thres0'}, ...
    {'ALL_All_woRelFeat'}, ...
    {'L2017-0219_MvsMonly2_thres0.1', 'WE2017-0112b_MvsFonly_thres0.1', ...
    'HB2018-0905_FvsFonly_thres0.1'}, ...
    {'L2017-0219_MvsMandMvsF_thres0.1', 'WE2017-0112b_MvsFandFvsF_thres0.1', ...
    'HB2018-0905_FvsFandMvsF_thres0.1'}, ...
    {'L2017-0219_MvsMandFvsF_thres0.1', 'WE2017-0112b_MvsFandMvsM_thres0.1', ...
    'HB2018-0905_FvsFandMvsM_thres0.1', }...
    {'L2017-0219_MvsMandFvsFandMvsF_thres0.1', 'WE2017-0112b_MvsMandMvsFandFvsF_thres0.1', ...
    'HB2018-0905_FvsFandMvsFandMvsM_thres0.1'}
    };
threshold_list = [0.1, 0.2, 0, 0.1, 0.1, 0.1, 0.1, 0.1];

jab_file_list = {jab_list, jab_list, jab_list, ...
    {'LungeNewRetrain2017-0219-woRelFeat_v2.jab', 'WingExtNew2017-0112b_woRelFeat_v2.jab', ...
    'HeadbuttNewRetrain2018-0905-woRelFeat.jab'}, ...
    {'LungeNewRetrain2017-0219_MvsMonly.jab', 'WingExtNew2017-0112b_MvsFonly.jab', ...
    'HeadbuttNewRetrain2018-0905_FvsFonly.jab'}, ...
    {'LungeNewRetrain2017-0219_MvsMandMvsF2.jab', 'WingExtNew2017-0112b_MvsFandFvsF.jab', ...
    'HeadbuttNewRetrain2018-0905_FvsFandMvsF.jab'}, ...
    {'LungeNewRetrain2017-0219_MvsMandFvsF.jab', 'WingExtNew2017-0112b_MvsFandMvsM.jab', ...
    'HeadbuttNewRetrain2018-0905_FvsFandMvsM.jab'}, ...
    {'LungeNewRetrain2017-0219_MvsMandFvsFandMvsF.jab', 'WingExtNew2017-0112b_MvsMandMvsFandFvsF.jab', ...
    'HeadbuttNewRetrain2018-0905_FvsFandMvsMandMvsF.jab'}};
for i=1:length(flymat_id_list)
    % Predict using the classifiers
    % JAABADetect(all_jaaba_folder_list, 'jabfiles', jab_file_list{i});
     
    thres = threshold_list(i); 
    override_jab_list = struct('behav_sel', {1;2;3}, 'jab_name', jab_file_list{i}');
    if i == 4
        is_NoRel = true; 
    else
        is_NoRel = false; 
    end
    
    if length(flymat_id_list{i}) == 1
        flymat_id = flymat_id_list{i}{1};
        behav_sel = [1,2,3];
        
        % Generate bout matches mat and bout-based bar and violin plots 
%         bout_match_str = strcat('bout_matches_', flymat_id, '.mat');
%         analyze_human_jaaba_annot_corr_v3(flymat_id, behav_sel, override_jab_list, is_NoRel, '', plot_bar_bout_wise, plot_violin_bout_wise, false, 'png');
%         analyze_human_jaaba_annot_corr_v3(flymat_id, behav_sel, override_jab_list, is_NoRel, bout_match_str, plot_bar_bout_wise, plot_violin_bout_wise, false, 'eps');
% 
%         % Generate bout-based false positive, false negative and true positive count 
%         compile_fp_fn_numbers(behav_sel, flymat_id, false);
        
        % For wing extension, generate corresponding frame-based plots
        frame_match_str = strcat('frame_matches_', flymat_id, '.mat');
        %analyze_human_jaaba_annot_corr_frame_wise(flymat_id, [1, 2, 3], override_jab_list, is_NoRel, '', plot_bar_frame_wise, plot_violin_frame_wise, 'png');
        %analyze_human_jaaba_annot_corr_frame_wise(flymat_id, 2, override_jab_list, is_NoRel, frame_match_str, plot_bar_frame_wise, plot_violin_frame_wise, 'eps');

        compile_fp_fn_numbers([1,2,3], flymat_id, true);
    else    
        for j=1:length(flymat_id_list{i})
            flymat_id = flymat_id_list{i}{j};
            
%             bout_match_str = strcat('bout_matches_', flymat_id, '.mat');
%             analyze_human_jaaba_annot_corr_v3(flymat_id, j, override_jab_list, is_NoRel, '', plot_bar_bout_wise, plot_violin_bout_wise, false, 'png');
%             analyze_human_jaaba_annot_corr_v3(flymat_id, j, override_jab_list, is_NoRel, bout_match_str, plot_bar_bout_wise, plot_violin_bout_wise, false, 'eps');
% 
%             compile_fp_fn_numbers(j, flymat_id, false);
            %if j == 2
                % For wing extension, generate corresponding frame-based plots
                frame_match_str = strcat('frame_matches_', flymat_id, '.mat');
%                 analyze_human_jaaba_annot_corr_frame_wise(flymat_id, j, override_jab_list, is_NoRel, '', plot_bar_frame_wise, plot_violin_frame_wise, 'png');
%                 analyze_human_jaaba_annot_corr_frame_wise(flymat_id, j, override_jab_list, is_NoRel, frame_match_str, plot_bar_frame_wise, plot_violin_frame_wise, 'eps');
                compile_fp_fn_numbers(j, flymat_id, true);
            %end
        end
    end
end