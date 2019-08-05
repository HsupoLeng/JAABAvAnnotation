function jaaba_score_all_train_frames = find_train_frame_jaaba_scores(...
    behav_name, jab_file_name)
    jab = load(jab_file_name, '-mat');
    jaaba_score_all_train_frames = struct('behav', [], 'none_behav', []);
    for i=1:length(jab.x.labels)
        % Skip the one experiment which cannot be classified
        if contains(jab.x.expDirNames{i}, '100215_2')
            continue;
        end
        exp_dir = strsplit(jab.x.expDirNames{i}, '\');
        rel_exp_dir = fullfile(exp_dir{6:end});
        try
            load(fullfile(rel_exp_dir, ...
                strcat('scores_', erase(behav_name, {'_','s'}),'.mat')));
        catch ME
            warning(ME.message);
            continue
        end
        for j=1:length(jab.x.labels(i).flies)
            t0s = jab.x.labels(i).t0s{j};
            t1s = jab.x.labels(i).t1s{j};
            fly_idx = jab.x.labels(i).flies(j);
            behavs = jab.x.labels(i).names{j};
            for k=1:length(t0s)
                scores_seg = allScores.scores{fly_idx}(t0s(k):t1s(k)-1)./...
                    allScores.scoreNorm;
                if isequal(behavs{k}, 'None')
                    behav_opt = 'none_behav';                    
                else
                    behav_opt = 'behav';
                end
                jaaba_score_all_train_frames.(behav_opt) = ...
                    [jaaba_score_all_train_frames.(behav_opt), scores_seg];
            end
        end
    end
end