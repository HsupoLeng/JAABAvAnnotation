function [num_positive_train_frames, num_negative_train_frames] = report_training_set_size(jab_file_path)
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
    num_positive_train_frames = sum(arrayfun(@(s) s.t1 - s.t0, all_positive_train_bouts));
    num_negative_train_frames = sum(arrayfun(@(s) s.t1 - s.t0, all_negative_train_bouts));

end