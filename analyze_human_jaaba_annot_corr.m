load('common-params-annot-analysis.mat');
load('FLYMAT_HumanAnnotation.mat');
load('FLYMAT_NewTrainingFiles-OriginalCopy.mat');

movie_name_pairs = cellfun(@(c) strsplit(c, '\'), movie_name_list, 'UniformOutput', false);
init_bout_matches_all_args = [behav_list; cell(1, length(behav_list))];
bout_matches_all = struct(init_bout_matches_all_args{:});

% Compile all the data we need 
for i=1:length(behav_list)
    % Create dynamic field name for this current behaviour
    t0s = strcat(behav_shorthands{i}, '_t0s');
    t1s = strcat(behav_shorthands{i}, '_t1s');
    scores = strcat(behav_shorthands{i}, '_scores');
    combined_score = strcat(behav_shorthands{i}, '_combined_score');
    t4s = strcat(behav_shorthands{i}, '_t4s');
    t5s = strcat(behav_shorthands{i}, '_t5s');
    
    % Create structure for bout matches in one type of behaviour
    % How we deal with no match situations:
    % 1) For a JAABA false negative, we average JAABA score in
    %    [annot_intsct_start, annot_intsct_end]. 
    % 2) For a JAABA false positive, set the annot_score as 0
    bout_matches_fields = {'movie', 'fly', ...
        'annot_union_start', 'annot_union_end', ...
        'annot_intsct_start', 'annot_intsct_end', 'annot_score', ...
        'jaaba_bout_start', 'jaaba_bout_end', 'jaaba_score_avg', ...
        'jaaba_score_min', 'jaaba_score_max', ...
        'jaaba_score_avg_normed', 'jaaba_score_min_normed', 'jaaba_score_max_normed', ...
        'multi_match', 'virtual_jaaba_match', 'invalid_match'};
    init_bout_matches_args = [bout_matches_fields; cell(1,length(bout_matches_fields))];
    bout_matches = struct(init_bout_matches_args{:});
    bout_matches(1) = [];
    
    % Extract annotated bouts
    % Currently, we take the union of two or more annotations to calculate
    % the start and end of the annotated bout. This is rather liberal, and
    % is aimed to make matching with JAABA bouts eaiser. 
    % When one annotator did not label a bout, while another did, NaN
    % values will appear, but they do not affect the min/max operation
    bout_idx = uint16(1);
    nan_init_cell_arr = cell(1,length(bout_matches_fields));
    nan_init_cell_arr(:) = {nan};
    init_bout_matches_args = [bout_matches_fields; nan_init_cell_arr];
    bout_match = struct(init_bout_matches_args{:});
    for j=1:length(flymatHumanAnnot)
        for k=1:length(flymatHumanAnnot(j).(combined_score))
            bout_match.movie = flymatHumanAnnot(j).movie;
            bout_match.fly = flymatHumanAnnot(j).fly;
            bout_match.annot_union_start = min(flymatHumanAnnot(j).(t0s)(k,:));
            bout_match.annot_union_end = max(flymatHumanAnnot(j).(t1s)(k,:));
            bout_match.annot_intsct_start = max(flymatHumanAnnot(j).(t0s)(k,:));
            bout_match.annot_intsct_end = min(flymatHumanAnnot(j).(t1s)(k,:));
            bout_match.annot_score = flymatHumanAnnot(j).(combined_score)(k);
            bout_matches(bout_idx) = bout_match;
            bout_idx = bout_idx + 1;
        end
    end
    
    % Iterate over annotated flies, find them in the JAABA flymat, match
    % JAABA bouts one-by-one with annotated bouts. If no match is found,
    % add the JAABA bout as a new entry in bout_matches
    for j=1:length(flymatHumanAnnot)
        fly_idx = find_fly_in_flymat(flymatAll, flymatHumanAnnot(j).movie, ...
            flymatHumanAnnot(j).fly, true);
        candidates = find_fly_in_flymat(bout_matches, ...
                flymatHumanAnnot(j).movie, flymatHumanAnnot(j).fly, false);
        for k=1:length(flymatAll(fly_idx).(t4s))
            jaaba_bout_start = flymatAll(fly_idx).(t4s)(k);
            jaaba_bout_end = flymatAll(fly_idx).(t5s)(k);     
            matched = 0;
            for m=1:length(candidates)
                % Declare a match if there is any overlap between annotated
                % bout and JAABA bout
                if (jaaba_bout_start <= bout_matches(candidates(m)).annot_union_start ...
                        && jaaba_bout_end > bout_matches(candidates(m)).annot_union_start) ...
                    || ...
                    (jaaba_bout_start > bout_matches(candidates(m)).annot_union_start ...
                        && jaaba_bout_start < bout_matches(candidates(m)).annot_union_end)
                    bout_matches(candidates(m)).jaaba_bout_start = jaaba_bout_start; 
                    bout_matches(candidates(m)).jaaba_bout_end = jaaba_bout_end; 
                    bout_matches(candidates(m)).virtual_jaaba_match = false;
                    bout_matches(candidates(m)).invalid_match = false;
                    if ~matched
                        matched = candidates(m);
                        bout_matches(candidates(m)).multi_match = false;
                    else
                        bout_matches(matched).multi_match = true;
                        bout_matches(candidates(m)).multi_match = true;
                    end
                end
            end
            % If no match is found for a JAABA bout, add it to bout_matches
            if ~matched
                bout_match = struct(init_bout_matches_args{:});
                bout_match.movie = flymatHumanAnnot(j).movie;
                bout_match.fly = flymatHumanAnnot(j).fly;
                bout_match.annot_union_start = nan;
                bout_match.annot_union_end = nan;
                bout_match.annot_score = 0;
                bout_match.jaaba_bout_start = jaaba_bout_start;
                bout_match.jaaba_bout_end = jaaba_bout_end;
                bout_match.virtual_jaaba_match = false;
                bout_match.invalid_match = false;
                bout_matches(bout_idx) = bout_match;
                bout_idx = bout_idx + 1;
            end
        end
    end
    
    % If no match is found for an existing human annotation, use human
    % intersection of human annotation (consensus) as a virtual JAABA bout 
    for j=1:length(bout_matches)
        if isnan(bout_matches(j).jaaba_bout_start)
            bout_matches(j).jaaba_bout_start = bout_matches(j).annot_intsct_start;
            bout_matches(j).jaaba_bout_end = bout_matches(j).annot_intsct_end;
            bout_matches(j).virtual_jaaba_match = true;
        end
    end
    
    % Sort the bouts by movie field
    bout_matches_cell = struct2cell(bout_matches);
    bout_matches_cell = permute(squeeze(bout_matches_cell), [2,1]);
    bout_matches_cell = sortrows(bout_matches_cell, 1:4);
    bout_matches_cell = permute(bout_matches_cell, [2,1]);
    bout_matches = cell2struct(bout_matches_cell, bout_matches_fields, 1)';
    
    % Iterate over the movies, load scores.mat, 
    % calculate the avg, min and max score for each JAABA bout entry 
    % in bout_matches
    for j=1:length(movie_name_list)
        load(fullfile(movie_name_list{j}, ...
            strcat(movie_name_pairs{j}{2},'_JAABA'), ...
            strcat('scores_', erase(behav_list{i}, {'_','s'}),'.mat')));
        bouts = uint16(find_fly_in_flymat(bout_matches, movie_name_pairs{j}{1}, [], false));
        for k=1:length(bouts)
            jaaba_bout_start = bout_matches(bouts(k)).jaaba_bout_start;
            jaaba_bout_end = bout_matches(bouts(k)).jaaba_bout_end;
            if ~isnan(jaaba_bout_start) && ~isnan(jaaba_bout_end)
                scores_slice = allScores.scores{bout_matches(bouts(k)).fly}(jaaba_bout_start:jaaba_bout_end-1);
                bout_matches(bouts(k)).jaaba_score_avg = mean(scores_slice);
                bout_matches(bouts(k)).jaaba_score_min = min(scores_slice);
                bout_matches(bouts(k)).jaaba_score_max = max(scores_slice);
                bout_matches(bouts(k)).jaaba_score_avg_normed = bout_matches(bouts(k)).jaaba_score_avg/allScores.scoreNorm;
                bout_matches(bouts(k)).jaaba_score_min_normed = bout_matches(bouts(k)).jaaba_score_min/allScores.scoreNorm;
                bout_matches(bouts(k)).jaaba_score_max_normed = bout_matches(bouts(k)).jaaba_score_max/allScores.scoreNorm;
            end
        end
    end
    
    bout_matches_all.(behav_list{i}) = bout_matches;
end

% Remove entries which failed the sanity check & adjacent ones affected
invalid_indxs = [1, 464; 1, 499; 1, 535; 1, 536; 1, 587; 1, 764; 2, 1293; 3, 505];
affected_indxs = [1, 461; 1, 500; 1, 537; 1, 591; 1, 766; 3, 504];
remove_indxs = [invalid_indxs; affected_indxs];
bout_masks = cell(1,3);
for i=1:length(bout_masks)
    bout_masks{i} = uint16(setdiff(1:length(bout_matches_all.(behav_list{i})), ...
        remove_indxs(remove_indxs(:,1)==i, 2)));
end

% Plot human annotation score vs. JAABA score for valid bout matches, for
% each behaviour type
stats_to_plot = {'avg', 'min', 'max'};
for i=1:length(behav_list)
    jaaba_score_stats_all_train_bouts = find_train_bout_jaaba_score_stats(...
        behav_list{i}, jab_list{i});
    annot_score = uint8([bout_matches_all.(behav_list{i})(bout_masks{i}).annot_score]);
    false_negat_idxs = [bout_matches_all.(behav_list{i})(bout_masks{i}).virtual_jaaba_match];
    % For annotation scores, the following two names just inherit
    % the template of corresponding jaaba score. Themselves are by no
    % means virtual or true. 
    virt_annot_score = annot_score(false_negat_idxs); 
    true_annot_score = annot_score(~false_negat_idxs);
    for k=1:length(stats_to_plot)
        field_to_plot = strcat('jaaba_score_', stats_to_plot{k}, '_normed');
        jaaba_score = [bout_matches_all.(behav_list{i})(bout_masks{i}).(field_to_plot)];
        virt_jaaba_score = jaaba_score(false_negat_idxs);
        true_jaaba_score = jaaba_score(~false_negat_idxs);

        figure();
        hold on;
        % =====
        % Scatter plot
        % x-axis offset to make distinguish between 2 catagories
        % =====

        scatter(double(true_annot_score)-0.1, true_jaaba_score, [], 'b');
        scatter(double(virt_annot_score)+0.1, virt_jaaba_score, [], 'm');
        legend('true or false positive', 'false negative', 'Location', 'southeast');
        
        % =====
        % Box plot
        % =====
        %{
        boxplot(true_jaaba_score, true_annot_score, 'Positions', double(0:6)-0.15, 'PlotStyle', 'compact', 'Widths', 0.2, 'Colors', 'b');
        boxplot(virt_jaaba_score, virt_annot_score, 'Positions', double(1:6)+0.15, 'PlotStyle', 'compact', 'Widths', 0.2, 'Colors', 'm');
        set(gca, 'xtick', 0:6);
        set(gca, 'xticklabel', 0:6);
        boxs = findobj(gca, 'Tag', 'Box');
        legend(boxs([8, 1]),  'true or false positive', 'false negative', ...
        'Location', 'southeast');
        %}
        % ===== Common plot settings =====
        line1 = line([-0.5, 6.5], repmat(mean(cell2mat(score_benchmrk_train.behav.(stats_to_plot{k}))), [1, 2]), ...
            'Color', 'k', 'LineStyle', '--', 'DisplayName', 'avg. over all positive{\it training} bouts');
        line2 = line([-0.5, 6.5], repmat(mean(cell2mat(score_benchmrk_train.none_behav.(stats_to_plot{k}))), [1, 2]), ...
            'Color', 'k', 'LineStyle', '-.', 'DisplayName', 'avg. over all negative{\it training} bouts');
        title({strcat(stats_to_plot{k}, '. of JAABA score vs. human combined score for test bouts'), ...
            strcat('behaviour:', erase(behav_list{i},'_'))});
        xlabel('Human annotation combined score \in [0,6]');
        ylabel(strcat(stats_to_plot{k}, '. of JAABA score in a bout'));
        axis([-0.5, 6.5, min([true_jaaba_score, virt_jaaba_score])-1.2, ...
            max([true_jaaba_score, virt_jaaba_score])+0.3]);
        hold off;
        saveas(gcf, strcat(strcat('jaaba_score_', stats_to_plot{k}, '_normed'), ...
        '_vs_annot_combined-',erase(behav_list{i},'_'), '-scatter.png'));
        % ===== Plot options end =====

    end
    
    % =====
    % Relative density plot
    % =====
    %{
    figure();
    false_negat_accu = hist(virt_annot_score, double(unique(virt_annot_score)));
    positive_accu = hist(true_annot_score, double(unique(true_annot_score)));
    collect_accu = [positive_accu(2:end)', false_negat_accu'];
    false_negat_ratio = false_negat_accu ./(false_negat_accu+positive_accu(2:end));
    yyaxis left
    bar_plt = bar(1:6, collect_accu, 'stacked', 'FaceAlpha', 0.5);
    yyaxis right
    l_plt = plot(1:6, false_negat_ratio);   
    l_plt.LineWidth = 1;
    title(strcat('JAABA bout classification count-', ...
        erase(behav_list{i},'_')));
    xlabel('Human annotation combined score \in [1,6]');
    yyaxis left
    ylabel('JAABA bout classification count');
    yyaxis right
    ylabel('false negative ratio');
    curr_lim = ylim;
    ylim([curr_lim(1), curr_lim(2)*1.2])
    set(bar_plt, {'FaceColor'}, {'b'; 'm'});
    legend({'true or false positives', 'false negatives'}, 'Location', 'northwest');
    saveas(gcf, strcat('jaaba_bout_classification_count-',erase(behav_list{i},'_'), '-bar.png'));
    %}
end

save('bout_matches_all.mat', 'bout_matches_all');