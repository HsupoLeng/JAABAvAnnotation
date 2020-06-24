% Alter training frames randomly and re-train for down-sampled classifiers 
load('common-params-annot-analysis.mat', 'behav_list', 'jab_list');
jab_root = 'D:\xubo\classifiers';
movie_root = 'D:\xubo\NewTrainingFiles-OriginalCopy\movies';
% frames_ratio_opts = [0.05, 0.25, 0.5, 0.75];
frames_ratio_opts = [nan, nan, nan]; % frame ratio to keep determined as equal to the base genotype pairing(m-m for lunge; m-f for wing extension; f-f for headbutt)
num_repeat = 10;

for iter=7:num_repeat 
    for j=1:length(frames_ratio_opts)
        for i=1:length(jab_list)
            if i ~= j
                continue;
            end
            
            [~, jab_id, ~] = fileparts(jab_list{i});
            [x, keep_frames_ratio] = alter_train_frames(fullfile(jab_root, jab_list{i}), frames_ratio_opts(j), false);
            if isnan(frames_ratio_opts(i))
                frames_ratio_opts(i) = keep_frames_ratio; 
            end
            
            for k=1:length(x.expDirNames)
                movie_path_parts = strsplit(x.expDirNames{k}, '\');
                x.expDirNames{k} = fullfile(movie_root, movie_path_parts{end-2:end}); 
            end

            ori_training_timestamp = x.classifierStuff.timeStamp;
            new_jab_path = fullfile(jab_root, '\rand_train_subset', sprintf('%s_percent_%d_repeat_%d.jab', jab_id, floor(frames_ratio_opts(j)*100), iter));
            save(new_jab_path, 'x');

            hObj = myStartJAABA(new_jab_path);

            % Train on modified training set
            simulate_keypress = struct('Modifier', 'control', 'Key', 't', 'Character', 't');
            hObj.KeyPressFcn(hObj, simulate_keypress);

            % Save project
            simulate_keypress = struct('Modifier', 'control', 'Key', 's', 'Character', 't');
            hObj.KeyPressFcn(hObj, simulate_keypress);

            % Verify that training is actually done
            load(new_jab_path, '-mat');
            new_training_timestamp = x.classifierStuff.timeStamp;
            assert(~isequal(ori_training_timestamp, new_training_timestamp), 'No training performed!\n');        
        end
    end
end