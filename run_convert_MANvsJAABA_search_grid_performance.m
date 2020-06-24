% Search grid of hyperparameters for precision and recall on the fully
% trained classifiers
close all;
load('common-params-annot-analysis.mat', 'all_jaaba_folder_list', 'jab_list');

file_id_list = {{'L2017-0219_All', 'WE2017-0112b_All', 'HB2018-0905_All'}, ...
    {'L2017-0219_All_woRelFeat', 'WE2017-0112b_All_woRelFeat', 'HB2018-0905_All_woRelFeat'}, ...
    {'L2017-0219_MvsMonly2', 'WE2017-0112b_MvsFonly', ...
    'HB2018-0905_FvsFonly'}, ...
    {'L2017-0219_MvsMandMvsF', 'WE2017-0112b_MvsFandFvsF', ...
    'HB2018-0905_FvsFandMvsF'}, ...
    {'L2017-0219_MvsMandFvsF', 'WE2017-0112b_MvsFandMvsM', ...
    'HB2018-0905_FvsFandMvsM', }...
    {'L2017-0219_MvsMandFvsFandMvsF', 'WE2017-0112b_MvsMandMvsFandFvsF', ...
    'HB2018-0905_FvsFandMvsFandMvsM'}
    }; % Identifiers for FLYMATs
threshold_list = -0.1:0.2:0.3; % Confidence thresholds
jab_file_list = {jab_list, ...
    {'LungeNewRetrain2017-0219-woRelFeat_v2.jab', 'WingExtNew2017-0112b_woRelFeat_v2.jab', ...
    'HeadbuttNewRetrain2018-0905-woRelFeat.jab'}, ...
    {'LungeNewRetrain2017-0219_MvsMonly.jab', 'WingExtNew2017-0112b_MvsFonly.jab', ...
    'HeadbuttNewRetrain2018-0905_FvsFonly.jab'}, ...
    {'LungeNewRetrain2017-0219_MvsMandMvsF2.jab', 'WingExtNew2017-0112b_MvsFandFvsF.jab', ...
    'HeadbuttNewRetrain2018-0905_FvsFandMvsF.jab'}, ...
    {'LungeNewRetrain2017-0219_MvsMandFvsF.jab', 'WingExtNew2017-0112b_MvsFandMvsM.jab', ...
    'HeadbuttNewRetrain2018-0905_FvsFandMvsM.jab'}, ...
    {'LungeNewRetrain2017-0219_MvsMandFvsFandMvsF.jab', 'WingExtNew2017-0112b_MvsMandMvsFandFvsF.jab', ...
    'HeadbuttNewRetrain2018-0905_FvsFandMvsMandMvsF.jab'}}; % JAABA classifiers, each corresponding to one FLYMAT

behav_sel = [1,2,3]; % Behavior labels, [lunge, wingExtension, headbutt]
is_frame_wise = false; % Flag to use frame-based performance metrics (precision and recall) or bout-based

for i=1:length(jab_file_list)
    if i ~= 1 
        continue; 
    end
    % Predict using the classifiers
    JAABADetect(all_jaaba_folder_list, 'jabfiles', jab_file_list{i});
    
    override_jab_list = struct('behav_sel', {1;2;3}, 'jab_name', jab_file_list{i}');
    % For my particular case, the second set of classifiers are w/o relative features, 
    % this requires changes to some field names, which is taken care of by
    % using this flag. Change as needed. 
    if i == 2
        is_NoRel = true; 
    else
        is_NoRel = false; 
    end
    threshold = threshold_list;
    
    for j=1:length(behav_sel) % For each of the three behaviours
        file_id = file_id_list{i}{j}; 
        
        % Each type of behavior have different post-processing parameters.
        % We vary them around the original values that we've been using.
        % See OrgData function. 
        if behav_sel(j) == 1 % lunge
            maxgap = 1:3; % Original 1
            minbehav = 1:2:5; % Original 3
            ourchoice_idx = sub2ind([length(threshold_list), length(maxgap), length(minbehav)], ...
                find(abs(threshold_list-0.1)<0.0001), find(maxgap==1), find(minbehav==3));
        elseif behav_sel(j) == 2 % wing extension
            maxgap = 1:2:5; % Original 3
            minbehav= 4:2:8; % Original 6
            ourchoice_idx = sub2ind([length(threshold_list), length(maxgap), length(minbehav)], ...
                find(abs(threshold_list-0.1)<0.0001), find(maxgap==3), find(minbehav==6));
        elseif behav_sel(j) == 3 % heatbutt
            maxgap = 1:3; % Original 1
            minbehav = 1:2:5; % Original 3
            ourchoice_idx = sub2ind([length(threshold_list), length(maxgap), length(minbehav)], ...
                find(abs(threshold_list-0.1)<0.0001), find(maxgap==1), find(minbehav==3));
        end
        
        % If performance metrics are frame-based, first compile MANvsJAABA,
        % then use MANvsJAABAduration_xubo to compute precision and recall,
        % and visualize the results; 
        % If performance metrics are bout-based, call
        % search_grid_performance_bout_wise to compute precision and
        % recall, which under-the-hood calls OrgData to generate temporary FLYMAT
        % for each set of parameter setting, and also calls
        % analyze_human_jaaba_annot_corr_v3 to generate bout_matches. 
        if is_frame_wise
            convert_to_MANvsJAABA(file_id, behav_sel(j), override_jab_list, is_NoRel);
            MANvsJAABAduration_xubo(sprintf('ManualvsJAABA-%s.mat', file_id), threshold, maxgap, minbehav, ourchoice_idx);
        else
            search_grid_performance_bout_wise(file_id, behav_sel(j), threshold, is_NoRel, maxgap, minbehav, override_jab_list, ourchoice_idx);
        end
    end
end