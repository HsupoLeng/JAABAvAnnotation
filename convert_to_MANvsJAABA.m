% Generate MANvsJAABA, a MATLAB .mat combining human annotation and JAABA
% scores. 
function convert_to_MANvsJAABA(file_id, behav_sel, override_jab_list, is_NoRel)
    load('common-params-annot-analysis.mat', ...
             'annot_file', 'behav_list', 'behav_shorthands', 'jab_list', 'groundtruth_movie_list');
    load('FLYMAT_HumanAnnotation_v3.mat', 'flymatHumanAnnot');
    
    if is_NoRel
        behav_list = {'LungeNewNoRel', 'WingExtNoRel', 'HeadbuttNew'};
    end
    if ~isempty(override_jab_list)
        for i=1:length(override_jab_list)
            jab_list{override_jab_list(i).behav_sel} = override_jab_list(i).jab_name;
        end
    end
    
    movie_name_pairs = cellfun(@(c) strsplit(c, '\'), groundtruth_movie_list, 'UniformOutput', false);
    
    for i=1:behav_sel
        t0s_field = strcat(behav_shorthands{i}, '_t0s');
        t1s_field = strcat(behav_shorthands{i}, '_t1s');
        scores_field = strcat(behav_shorthands{i}, '_scores');
        
        manvjaaba = struct('movie', []);
        for j=1:length(groundtruth_movie_list)
            manvsjaaba_onemovie = struct('name', '', 'flies', [], 'observer1', {}, 'observer2', {}, ...
                'frames', nan, 'allScores', []);
            manvsjaaba_onemovie(1).name = movie_name_pairs{j}{2};
            manvsjaaba_onemovie(1).flies = [flymatHumanAnnot(strcmp({flymatHumanAnnot.movie}, movie_name_pairs{j}{1})).fly];
            manvsjaaba_onemovie(1).observer1 = cell(1, length(manvsjaaba_onemovie.flies));
            for k=1:length(manvsjaaba_onemovie.flies)
                fly_idx = find_fly_in_flymat(flymatHumanAnnot, movie_name_pairs{j}{1}, manvsjaaba_onemovie.flies(k), false);
                manvsjaaba_onemovie.observer1{k} = [flymatHumanAnnot(fly_idx).(t0s_field)(:,1), flymatHumanAnnot(fly_idx).(t1s_field)(:,1), ...
                    flymatHumanAnnot(fly_idx).(scores_field)(:,1)]';
                manvsjaaba_onemovie.observer1{k}(:, isnan(manvsjaaba_onemovie.observer1{k}(3, :))) = [];
                manvsjaaba_onemovie.observer2{k} = [flymatHumanAnnot(fly_idx).(t0s_field)(:,2), flymatHumanAnnot(fly_idx).(t1s_field)(:,2), ...
                    flymatHumanAnnot(fly_idx).(scores_field)(:,2)]';
                manvsjaaba_onemovie.observer2{k}(:, isnan(manvsjaaba_onemovie.observer2{k}(3, :))) = [];
            end
            load(fullfile(groundtruth_movie_list{j}, ...
                        strcat(movie_name_pairs{j}{2},'_JAABA'), ...
                        strcat('scores_', erase(behav_list{i}, {'_','s'}),'.mat')), 'allScores');     
            num_frames = unique(cellfun(@length, allScores.scores));
            if length(num_frames) == 1
                manvsjaaba_onemovie.frames = num_frames; 
            else
                fprintf('Some flies have different number of predicted frames in movie %s, behavior %s.\n', movie_name_pairs{j}{1}, behav_list{i});
                fprintf('Padding with -Inf at the end to make length equal...\n');
                max_num_frames = max(num_frames);
                for k=1:length(allScores.scores)
                    if length(allScores.scores{k}) ~= max_num_frames
                        allScores.scores{k} = padarray(allScores.scores{k}, [0, max_num_frames-length(allScores.scores{k})], -inf, 'post');
                    end
                end
                manvsjaaba_onemovie.frames = max_num_frames; 
            end
            
            manvsjaaba_onemovie.allScores = allScores;
            manvjaaba(j).movie = manvsjaaba_onemovie;
        end
 
        % save the result
        save(sprintf('ManualvsJAABA-%s.mat', file_id), 'manvjaaba');
    end
end