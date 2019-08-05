annot_file = 'Copy of ManualvsJAABA-MasterComparison-SelectedTrainingMovie-v3.xlsx';
behaviour_list = {'Lunges_New', 'WingExt_New', 'Headbutt_New'};
movie_name_list = {'033015_NPF3-CsChrimsonattP40MG-2', ...
    '041815_1_m_g_Otd-FLPoChrimsonmvenusattP2_f2a20', ...
    '042015_12_m_g_Otd-FLPoChrimsonmvenusattP2_f10a10', ...
    '042415_4_m_g_Otd-FLPoChrimsonTdTomattP2_CsHeis_f2a20', ...
    '042815_assay1', ...
    '042815_assay4', ...
    '050815_assay9', ...
    '082615_CSMH_SF', ...
    '100815_4', ...
    '2016_02_15_CsMH_M_SH2', ...
    '2016_02_15_CsMH_M_SH3'
};

% Make valid field names from behaviour names and movie names
behav_name_pairs = struct('fieldName', matlab.lang.makeValidName(behaviour_list), 'oriName', behaviour_list);
movie_name_pairs = struct('fieldName', matlab.lang.makeValidName(movie_name_list), 'oriName', movie_name_list);

movieLevel_struct_args = [{movie_name_pairs(:).fieldName}; cell(1, size(movie_name_pairs,2))];
movieLevel_struct = struct(movieLevel_struct_args{:});
behavLevel_struct_args = [{behav_name_pairs(:).fieldName}; repmat({movieLevel_struct}, 1, size(behav_name_pairs,2))];
all_behav_annot = struct(behavLevel_struct_args{:});

% Create structure conforming with the format of flymatAll
init_fields = {'movie', 'fly', 'L_t0s', 'L_t1s', 'L_scores', 'L_combined_score', ...
    'HB_t0s', 'HB_t1s', 'HB_scores', 'HB_combined_score', 'WE_t0s', 'WE_t1s', ...
    'WE_scores', 'WE_combined_score'};
init_flymat_args = [init_fields; cell(1, length(init_fields))];
flymatHumanAnnot = struct(init_flymat_args{:});


for i=1:length(behaviour_list)
    [~,~,raw_xls] = xlsread(annot_file, behav_name_pairs(i).oriName);
    raw_xls = cellfun(@(x) num2str(x), raw_xls, 'UniformOutput', false); % Convert all fields to str 
    raw_xls = raw_xls(sum(cellfun(@(c) ~strcmp(c, 'NaN'), raw_xls), 2)>0, :); % Remove rows where entries are all NaN
    
    for j=1:length(movie_name_list)
        movie_start_idx = find(contains(raw_xls, movie_name_pairs(j).oriName));
        if j<length(movie_name_list)
            movie_end_idx = find(contains(raw_xls, movie_name_pairs(j+1).oriName));
        else
            movie_end_idx = numel(raw_xls);
        end
        [movie_start_row, movie_start_col] = ind2sub(size(raw_xls), movie_start_idx);
        [~, movie_end_col] = ind2sub(size(raw_xls), movie_end_idx);
        movie_fields = raw_xls(movie_start_row:size(raw_xls,1), movie_start_col:movie_end_col-1);
        
        fly_start_indices = find(contains(movie_fields, 'Fly'));
        fly_end_indices = find(contains(movie_fields, 'Human cobmbined score'));
        num_of_flies = {movie_fields{fly_start_indices}};
        fly_numbers = regexp(num_of_flies, '[0-9]+', 'match');
        flyLevel_struct = struct('fly_number', fly_numbers, ...
            't0s', cell(1, length(fly_numbers)), ...
            't1s', cell(1, length(fly_numbers)), ...
            'scores', cell(1, length(fly_numbers)), ...
            'human_combined_score', cell(1, length(fly_numbers)));
        
        for k=1:length(num_of_flies)
            [fly_start_row, fly_start_col] = ind2sub(size(movie_fields), fly_start_indices(k));
            [~, fly_end_col] = ind2sub(size(movie_fields), fly_end_indices(k));
            fly_fields = movie_fields(fly_start_row:size(movie_fields,1), fly_start_col:fly_end_col);
            valid_fly_fields = fly_fields(sum(cellfun(@(c) ~strcmp(c, 'NaN'), fly_fields), 2)>0, :);
            annot_start_indices = find(contains(valid_fly_fields, 'Start'));
            annot_fields = valid_fly_fields(annot_start_indices(1):end, :);
            true_positive_fields = annot_fields(cellfun(@(c) str2double(c), annot_fields(:, end)) > 0, :);
            cols_to_keep = boolean(ones(1, size(true_positive_fields, 2)));
            cols_to_keep(4:4:end) = 0;
            true_positive_fields = cellfun(@(c) str2double(c), true_positive_fields(:, cols_to_keep));
            for m=1:size(true_positive_fields, 1)
                flyLevel_struct(k).t0s = true_positive_fields(:, 1:3:end-1);
                flyLevel_struct(k).t1s = true_positive_fields(:, 2:3:end-1);
                flyLevel_struct(k).scores = true_positive_fields(:,3:3:end-1);
                flyLevel_struct(k).human_combined_score = true_positive_fields(:,end);
            end   
        end
        all_behav_annot = setfield(all_behav_annot, ...
            behav_name_pairs(i).fieldName, ...
            movie_name_pairs(j).fieldName, flyLevel_struct);
    end  
end

total_num_of_flies = 0;
%for i=1:length(movie_name_pairs)
%flymat_human_annot = struct([{'movie', 'fly', 'L_t0s', 'L_t1s', 'HB_t0s', 'HB_t1s', 'WE_t0s', 'WE_t1s'};
%    repmat(cell(1, ))])
save('all_behaviour_annotations.mat', 'all_behav_annot', 'behav_name_pairs', ...
    'movie_name_pairs');