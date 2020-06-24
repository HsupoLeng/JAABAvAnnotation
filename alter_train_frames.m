function [x, keep_frames_ratio] = alter_train_frames(jab_file_path, keep_frames_ratio, save_jab)
    jab = load(jab_file_path, '-mat');
    [jab_path, jab_id, ~] = fileparts(jab_file_path);
    behav_name = regexp(jab_id, '([a-zA-Z]+)(\d+)', 'tokens');
    behav_name = erase(behav_name{1}{1}, 'Retrain');
    
    % Compile positive and negative training bouts
    jaaba_all_train_bouts = struct('behav', nan, 'non_behav', nan);
    field_names = fieldnames(jaaba_all_train_bouts);
    for i=1:length(jab.x.labels)
        jaaba_all_train_bouts(i).behav = cell(size(jab.x.labels(i).flies));
        jaaba_all_train_bouts(i).non_behav = cell(size(jab.x.labels(i).flies));
        for k=1:length(jab.x.labels(i).flies)
            t0s = jab.x.labels(i).t0s{k};
            t1s = jab.x.labels(i).t1s{k};
            names = jab.x.labels(i).names{k};

            positive_label_mask = strcmp(names, behav_name);
            
            positive_ts = [t0s(positive_label_mask)', t1s(positive_label_mask)'];
            try 
                positive_ts = array2table(positive_ts, 'VariableNames', {'t0', 't1'});
                jaaba_all_train_bouts(i).behav{k} = table2struct(positive_ts);
            catch 
                jaaba_all_train_bouts(i).behav{k} = struct('t0', nan, 't1', nan);
                jaaba_all_train_bouts(i).behav{k}(1) = [];
                jaaba_all_train_bouts(i).behav{k} = jaaba_all_train_bouts(i).behav{k}';
            end

            negative_ts = [t0s(~positive_label_mask)', t1s(~positive_label_mask)'];
            try 
                negative_ts = array2table(negative_ts, 'VariableNames', {'t0', 't1'});
                jaaba_all_train_bouts(i).non_behav{k} = table2struct(negative_ts);
            catch 
                jaaba_all_train_bouts(i).non_behav{k} = struct('t0', nan, 't1', nan);
                jaaba_all_train_bouts(i).non_behav{k}(1) = [];
                jaaba_all_train_bouts(i).non_behav{k} = jaaba_all_train_bouts(i).non_behav{k}';
            end
        end
    end
    
    all_positive_train_bouts = vertcat(jaaba_all_train_bouts(:).behav);
    all_positive_train_bouts = vertcat(all_positive_train_bouts{:});
    all_negative_train_bouts = vertcat(jaaba_all_train_bouts(:).non_behav);
    all_negative_train_bouts = vertcat(all_negative_train_bouts{:});
    num_positive_train_frames_all = sum(arrayfun(@(s) s.t1 - s.t0, all_positive_train_bouts));
    num_negative_train_frames_all = sum(arrayfun(@(s) s.t1 - s.t0, all_negative_train_bouts));
    
    fprintf('Behavior %s - Number of positive train frames: %d, negative train frames: %d\n', ...
        behav_name, num_positive_train_frames_all, num_negative_train_frames_all);

    if isnan(keep_frames_ratio)
        if strcmp(behav_name, 'LungeNew')
            movie_mask = contains(jab.x.expDirNames, ...
                {'033015_NPF3', '041815_1_m_g_Otd', '042015_12_m_g_Otd', ...
                '042415_4_m_g_Otd', '042815_assay1', '2016_02_15_CsMH_M_SH2', ...
                '2016_02_15_CsMH_M_SH3', '2016-0824_TRH-BL3839', '2016-0824_15A01-AD_71G01-DBD', ...
                '160928_trial1'});
        elseif strcmp(behav_name, 'WingExtNew')
            movie_mask = contains(jab.x.expDirNames, {'042415_Jing', '042815_assay4'});
        elseif strcmp(behav_name, 'HeadbuttNew')
            movie_mask = contains(jab.x.expDirNames, {'050815_assay9', '082615_CSMH_SF', '20160826_TRH-KayserDark', '2017-1125_6'});
        end
        basal_positive_train_bouts = vertcat(jaaba_all_train_bouts(movie_mask).behav);
        basal_positive_train_bouts = vertcat(basal_positive_train_bouts{:});
        basal_negative_train_bouts = vertcat(jaaba_all_train_bouts(movie_mask).non_behav);
        basal_negative_train_bouts = vertcat(basal_negative_train_bouts{:});
        num_positive_train_frames_basal = sum(arrayfun(@(s) s.t1 - s.t0, basal_positive_train_bouts));
        num_negative_train_frames_basal = sum(arrayfun(@(s) s.t1 - s.t0, basal_negative_train_bouts));
        keep_frames_ratio = sum(num_positive_train_frames_basal + num_negative_train_frames_basal)/...
            sum(num_positive_train_frames_all + num_negative_train_frames_all);
    end
    
    % Convert bouts to frames, and keep the specified ratio of frames. The
    % frames kept are chosen randomly
    jaaba_new_train_bouts = struct('behav', nan, 'non_behav', nan);
    for i=1:size(jaaba_all_train_bouts, 2)
        for j=1:length(field_names) 
            jaaba_new_train_bouts(i).(field_names{j}) = cell(size(jaaba_all_train_bouts(i).(field_names{j})));
            for k=1:length(jaaba_new_train_bouts(i).(field_names{j}))
            % Convert bouts into frames
                ori_bouts = jaaba_all_train_bouts(i).(field_names{j}){k};
                ori_frames = arrayfun(@(s) (s.t0:(s.t1-1))', ori_bouts, 'UniformOutput', false);
                ori_frames = vertcat(ori_frames{:});
                if ~isempty(ori_frames)
                    % Randomly select frames to keep
                    keep_frames_idx = randperm(length(ori_frames), ceil(keep_frames_ratio*length(ori_frames)));
                    new_frames = ori_frames(sort(keep_frames_idx));

                    new_changepoint_mask = diff(new_frames) > 1;
                    new_boundaries = [new_frames([true; new_changepoint_mask]), new_frames([new_changepoint_mask; true])+1];
                    new_boundaries = array2table(new_boundaries, 'VariableNames', {'t0', 't1'});
                    jaaba_new_train_bouts(i).(field_names{j}){k} = table2struct(new_boundaries);
                else
                    jaaba_new_train_bouts(i).(field_names{j}){k} = struct('t0', nan, 't1', nan);
                    jaaba_new_train_bouts(i).(field_names{j}){k}(1) = [];
                    jaaba_new_train_bouts(i).(field_names{j}){k} = jaaba_new_train_bouts(i).(field_names{j}){k}';
                end
            end
        end
    end
    
    new_positive_train_bouts = vertcat(jaaba_new_train_bouts(:).behav);
    new_positive_train_bouts = vertcat(new_positive_train_bouts{:});
    new_negative_train_bouts = vertcat(jaaba_new_train_bouts(:).non_behav);
    new_negative_train_bouts = vertcat(new_negative_train_bouts{:});
    num_positive_train_frames_new = sum(arrayfun(@(s) s.t1 - s.t0, new_positive_train_bouts));
    num_negative_train_frames_new = sum(arrayfun(@(s) s.t1 - s.t0, new_negative_train_bouts));
    
    fprintf('Behavior %s - Number of positive train frames: %d, negative train frames: %d\n', ...
        behav_name, num_positive_train_frames_new, num_negative_train_frames_new);
    fprintf('Requested training set percentage %d%%; actual percentage: positive frames %.2f%%, negative frames %.2f%%\n', ...
        floor(keep_frames_ratio*100), num_positive_train_frames_new*100/num_positive_train_frames_all, ...
        num_negative_train_frames_new*100/num_negative_train_frames_all)
    
    new_jab = jab;
    for i=1:length(new_jab.x.labels)
        for k=1:length(new_jab.x.labels(i).flies)
            new_jab.x.labels(i).t0s{k} = [[jaaba_new_train_bouts(i).behav{k}(:).t0], [jaaba_new_train_bouts(i).non_behav{k}(:).t0]];
            new_jab.x.labels(i).t1s{k} = [[jaaba_new_train_bouts(i).behav{k}(:).t1], [jaaba_new_train_bouts(i).non_behav{k}(:).t1]];
            new_jab.x.labels(i).names{k} = [repmat({behav_name}, 1, length(jaaba_new_train_bouts(i).behav{k})), ...
                repmat({'None'}, 1, length(jaaba_new_train_bouts(i).non_behav{k}))];
            
            new_jab.x.labels(i).imp_t0s{k} = sort(new_jab.x.labels(i).t0s{k});
            new_jab.x.labels(i).imp_t1s{k} = sort(new_jab.x.labels(i).t1s{k});
            if ~isempty(new_jab.x.labels(i).t0s{k})
                now_val = now; 
                new_jab.x.labels(i).timestamp{k} = repmat(now_val, 1, length(new_jab.x.labels(i).t0s{k}));
                new_jab.x.labels(i).timelinetimestamp{k} = struct(behav_name, now_val);
            end
        end
    end
    
    x = new_jab.x;
    if save_jab
        save(fullfile(jab_path, '\rand_train_subset', sprintf('%s_percent_%d.jab', jab_id, floor(keep_frames_ratio*100))), 'x');
    end
end