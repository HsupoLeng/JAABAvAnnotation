function jaaba_score_stats_all_train_bouts = find_train_bout_jaaba_score_stats(...
    behav_name, jab_file_name)
    jab = load(jab_file_name, '-mat');
    jaaba_score_stats_all_train_bouts = struct('behav', nan, 'none_behav', nan);
    jaaba_score_stats_all_train_bouts.behav = struct('avg', cell(1), 'min', cell(1), 'max', cell(1));
    jaaba_score_stats_all_train_bouts.none_behav = struct('avg', cell(1), 'min', cell(1), 'max', cell(1));
    train_bouts_idx = struct('behav', uint16(1), 'none_behav', uint16(1));
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
                jaaba_score_stats_all_train_bouts.(behav_opt).avg{train_bouts_idx.(behav_opt)} = ...
                    mean(scores_seg);
                jaaba_score_stats_all_train_bouts.(behav_opt).min{train_bouts_idx.(behav_opt)} = ...
                    min(scores_seg);
                jaaba_score_stats_all_train_bouts.(behav_opt).max{train_bouts_idx.(behav_opt)} = ...
                    max(scores_seg);
                train_bouts_idx.(behav_opt) = train_bouts_idx.(behav_opt) + 1;
            end
        end
    end
end