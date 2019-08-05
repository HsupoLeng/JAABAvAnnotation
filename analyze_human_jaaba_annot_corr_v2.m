function analyze_human_jaaba_annot_corr_v2(flymat_file, behav_sel, override_jab_list, thres, plot_bar, plot_violin, plot_box, img_format)
    load('common-params-annot-analysis.mat', ...
         'annot_file', 'behav_list', 'behav_shorthands', 'jab_list', 'movie_name_list');
    load('FLYMAT_HumanAnnotation_new.mat', 'flymatHumanAnnot');
    load(flymat_file, 'flymatAll');

    [~, flymat_file_name, ~] = fileparts(flymat_file);
    flymat_file_name_components = strsplit(flymat_file_name, '_');
    flymat_id = join(flymat_file_name_components(end-1:end), '_'); 
    flymat_id = flymat_id{:};
    movie_name_pairs = cellfun(@(c) strsplit(c, '\'), movie_name_list, 'UniformOutput', false);
    %behav_list = {'LungeNewNoRel', 'WingExtNoRel', 'HeadbuttNew'};
    if ~isempty(behav_sel)
        for i=1:length(behav_list)
            if ismember(i, behav_sel)
                jab_list{i} = override_jab_list{find(behav_sel==i, 1)};
            else
                jab_list{i} = '';
            end
        end
    end
    
    init_bout_matches_all_args = [behav_list; cell(1, length(behav_list))];
    bout_matches_all = struct(init_bout_matches_all_args{:});
    
    if contains(img_format, 'eps')
        img_format = 'epsc';
    else
        contains(img_format, 'png')
        img_format = 'png';
    end

    % Compile all the data we need 
    for i=1:length(behav_list)
        if ~ismember(i, behav_sel)
            continue;
        end
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
                            % bout_matches(matched).multi_match = true;
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

        % If no match is found for an existing human annotation, use
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
                strcat('scores_', erase(behav_list{i}, {'_','s'}),'.mat')), 'allScores');
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
    % to-do: identify and remove the entry affected by the last invalid entry. 
    invalid_indxs = sanity_check_humanAnnot(behav_list, bout_matches_all);
    affected_indxs = [0, 0];
    for i=1:size(invalid_indxs, 1)
        mismatch_range = [-2, -1, 1, 2];
        for j=1:length(mismatch_range)
            if bout_matches_all.(behav_list{invalid_indxs(i, 1)})(uint16(invalid_indxs(i, 2)+mismatch_range(j))).multi_match == 1 ...
                    && ~any(invalid_indxs(:, 2) == invalid_indxs(i, 2)+mismatch_range(j)) ...
                    && ~any(affected_indxs(:, 2) == invalid_indxs(i, 2)+mismatch_range(j)) 
                affected_indxs = [affected_indxs; [invalid_indxs(i, 1), invalid_indxs(i, 2)+mismatch_range(j)]];
            end
        end
    end
    affected_indxs(1,:) = []; 
    remove_indxs = [invalid_indxs; affected_indxs];
    bout_masks = cell(1,length(behav_list));
    for i=1:length(bout_masks)
        bout_masks{i} = uint16(setdiff(1:length(bout_matches_all.(behav_list{i})), ...
            remove_indxs(remove_indxs(:,1)==i, 2)));
    end

    % Plot human annotation score vs. JAABA score for valid bout matches, for
    % each behaviour type
    stats_to_plot = {'avg', 'min', 'max'};
    human_annot_scores = 0:6;
    xtick_label_cell_arr = strtrim(cellstr(num2str(human_annot_scores')));
    xtick_label_cell_arr{end+1} = 'training';
    recall_prec_args = [behav_list; cell(1, length(behav_list))];
    recalls = struct(recall_prec_args{:});
    precisions = struct(recall_prec_args{:});
    for i=1:length(behav_list)
        if ~ismember(i, behav_sel)
            continue;
        end
        score_benchmrk_train = find_train_bout_jaaba_score_stats(...
            behav_list{i}, jab_list{i});
        annot_score = uint8([bout_matches_all.(behav_list{i})(bout_masks{i}).annot_score]);
        false_negat_idxs = [bout_matches_all.(behav_list{i})(bout_masks{i}).virtual_jaaba_match];
        multi_match_idxs = [bout_matches_all.(behav_list{i})(bout_masks{i}).multi_match];
        multi_match_idxs(isnan(multi_match_idxs)) = 0;
        
        % Calculate recalls without counting bout matches with
        % multi_match=1
        virt_annot_score_no_mm = annot_score(bitand(false_negat_idxs, ~multi_match_idxs)); 
        true_annot_score_no_mm = annot_score(bitand(~false_negat_idxs, ~multi_match_idxs));
        false_negat_accum_no_mm = hist(virt_annot_score_no_mm, human_annot_scores);
        positive_accum_no_mm = hist(true_annot_score_no_mm, human_annot_scores);
        recalls.(behav_list{i}) = sum(positive_accum_no_mm(2:end))/(sum(false_negat_accum_no_mm(2:end))+sum(positive_accum_no_mm(2:end)));
        collect_accum_no_mm = [positive_accum_no_mm(2:end)', false_negat_accum_no_mm(2:end)'];
        recall_rates = collect_accum_no_mm(:,1)./sum(collect_accum_no_mm, 2);
           
        % When plotting and calculating precisions, we count bout matches with
        % multi_match=1
        virt_annot_score = annot_score(false_negat_idxs); 
        true_annot_score = annot_score(~false_negat_idxs);
        false_negat_accum = hist(virt_annot_score, human_annot_scores);
        positive_accum = hist(true_annot_score, human_annot_scores);
        collect_accum = [positive_accum(2:end)', false_negat_accum(2:end)'];
        test_violin_norm = max([false_negat_accum, positive_accum]);
        positive_train_size = length(score_benchmrk_train.behav.avg);
        negat_train_size = length(score_benchmrk_train.none_behav.avg);
        train_violin_norm = max([positive_train_size, negat_train_size]);
        precisions.(behav_list{i}) = sum(positive_accum(2:end))/sum(positive_accum);
        for k=1:length(stats_to_plot)
            field_to_plot = strcat('jaaba_score_', stats_to_plot{k}, '_normed');
            jaaba_score = [bout_matches_all.(behav_list{i})(bout_masks{i}).(field_to_plot)];
            virt_jaaba_score = jaaba_score(false_negat_idxs);
            true_jaaba_score = jaaba_score(~false_negat_idxs);
            positive_train_jaaba_score = cell2mat(score_benchmrk_train.behav.(stats_to_plot{k}));
            negative_train_jaaba_score = cell2mat(score_benchmrk_train.none_behav.(stats_to_plot{k}));
            % =====
            % Box plot
            % =====
            if plot_box
                figure();
                hold on;

                boxplot(true_jaaba_score, true_annot_score, 'Positions', double(human_annot_scores)-0.15, 'PlotStyle', 'compact', 'Widths', 0.2, 'Colors', 'b');
                boxplot(virt_jaaba_score, virt_annot_score, 'Positions', double(human_annot_scores(2:end))+0.15, 'PlotStyle', 'compact', 'Widths', 0.2, 'Colors', 'm');
                boxplot(positive_train_jaaba_score, 'Positions', 7-0.15, 'PlotStyle', 'compact', 'Widths', 0.2, 'Colors', [0, 0.5, 0.5]);
                boxplot(negative_train_jaaba_score, 'Positions', 7+0.15, 'PlotStyle', 'compact', 'Widths', 0.2, 'Colors', [0.75, 0.25, 0]);
                boxs = findobj(gca, 'Tag', 'Box');
                legend(boxs([10, 3, 2, 1]),  {'true or false positive', 'false negative', 'positive training', 'negative training'},...
                'NumColumns', 2, 'Location', 'southwest');
                offset_hack = 0;
                
                all_plotted_scores = [true_jaaba_score, virt_jaaba_score, positive_train_jaaba_score, negative_train_jaaba_score];
                axis([-0.5+offset_hack, 7.5+offset_hack, min(all_plotted_scores)-1.2, ...
                    max(all_plotted_scores)+0.3]);
                line([-0.5, 6.5]+offset_hack, repmat(mean(cell2mat(score_benchmrk_train.behav.(stats_to_plot{k}))), [1, 2]), ...
                    'Color', 'k', 'LineStyle', '--', 'DisplayName', 'avg. over all positive{\it training} bouts');
                line([-0.5, 6.5]+offset_hack, repmat(mean(cell2mat(score_benchmrk_train.none_behav.(stats_to_plot{k}))), [1, 2]), ...
                    'Color', 'k', 'LineStyle', '-.', 'DisplayName', 'avg. over all negative{\it training} bouts');
                ylim_tmp = ylim; 
                train_bkg = area([6.5+offset_hack, 7.5+offset_hack], [ylim_tmp(2),ylim_tmp(2)], ylim_tmp(1), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                train_bkg.Annotation.LegendInformation.IconDisplayStyle = 'off';
                title({strcat(stats_to_plot{k}, '. of JAABA score vs. human combined score for test bouts'), ...
                    strcat('behaviour:', erase(behav_list{i},'_'))});
                xlabel('Human annotation combined score \in [0,6]');
                ylabel(strcat(stats_to_plot{k}, '. of JAABA score in a bout'));
                set(gca, 'xtick', [human_annot_scores, 7]++offset_hack);
                set(gca, 'xticklabel', xtick_label_cell_arr);
                hold off;
                set(gcf,'renderer','Painters');
                saveas(gcf, strcat(strcat('jaaba_score_', stats_to_plot{k}, '_normed'), ...
                    '_vs_annot_combined-',erase(behav_list{i},'_'), '-box-thres', num2str(thres), '.', erase(img_format, 'c')), img_format);
                
            end 
            % =====
            % Violin plot
            % =====
            if plot_violin
                figure();
                hold on;
                v1 = violinplot(true_jaaba_score, double(true_annot_score)-0.3, 'ViolinColor', [0, 0, 1], 'ShowData', false, 'ViolinNorm', test_violin_norm);
                v2 = violinplot([-10, virt_jaaba_score], double([0, virt_annot_score]+0.3), 'ViolinColor', [1, 0, 1], 'ShowData', false, 'ViolinNorm', test_violin_norm);
                v3 = violinplot([repmat(-10, 1, 7), positive_train_jaaba_score], ...
                    [0:6, repmat(7, 1, length(positive_train_jaaba_score))], 'ViolinColor', [0, 0.5, 0.5], 'ShowData', false, 'ViolinNorm', train_violin_norm);
                v4 = violinplot([repmat(-10, 1, 7), negative_train_jaaba_score], ...
                    [0:6, repmat(7, 1, length(negative_train_jaaba_score))], 'ViolinColor', [0.75, 0.25, 0], 'ShowData', false, 'ViolinNorm', train_violin_norm);
                set(gca, 'xtick', [human_annot_scores, 7]+1);
                set(gca, 'xticklabel', xtick_label_cell_arr);
                legend([v1(end).ViolinPlot, v2(2).ViolinPlot, v3(end).ViolinPlot, v4(end).ViolinPlot], {'true or false positive', 'false negative', 'positive training', 'negative training'},...
                'NumColumns', 2, 'Location', 'southwest');
                offset_hack = 1;

                all_plotted_scores = [true_jaaba_score, virt_jaaba_score, positive_train_jaaba_score, negative_train_jaaba_score];
                axis([-0.5+offset_hack, 7.5+offset_hack, min(all_plotted_scores)-1.2, ...
                    max(all_plotted_scores)+0.3]);
                line([-0.5, 6.5]+offset_hack, repmat(mean(cell2mat(score_benchmrk_train.behav.(stats_to_plot{k}))), [1, 2]), ...
                    'Color', 'k', 'LineStyle', '--', 'DisplayName', 'avg. over all positive{\it training} bouts');
                line([-0.5, 6.5]+offset_hack, repmat(mean(cell2mat(score_benchmrk_train.none_behav.(stats_to_plot{k}))), [1, 2]), ...
                    'Color', 'k', 'LineStyle', '-.', 'DisplayName', 'avg. over all negative{\it training} bouts');
                ylim_tmp = ylim; 
                train_bkg = area([6.5+offset_hack, 7.5+offset_hack], [ylim_tmp(2),ylim_tmp(2)], ylim_tmp(1), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                train_bkg.Annotation.LegendInformation.IconDisplayStyle = 'off';
                title({strcat(stats_to_plot{k}, '. of JAABA score vs. human combined score for test bouts'), ...
                    strcat('behaviour:', erase(behav_list{i},'_'))});
                xlabel('Human annotation combined score \in [0,6]');
                ylabel(strcat(stats_to_plot{k}, '. of JAABA score in a bout'));
                set(gca, 'xtick', [human_annot_scores, 7]++offset_hack);
                set(gca, 'xticklabel', xtick_label_cell_arr);
                hold off;
                set(gcf,'renderer','Painters');
                saveas(gcf, strcat(strcat('jaaba_score_', stats_to_plot{k}, '_normed'), ...
                '_vs_annot_combined-',erase(behav_list{i},'_'), '-violin-thres', num2str(thres), '.', erase(img_format, 'c')), img_format);
            end
        end 

        % =====
        % Bar plot
        % =====
        if plot_bar
            figure();
            hold on;
            yyaxis left
            bar_plt = bar(1:6, collect_accum, 'stacked', 'FaceAlpha', 0.5);
            false_positive_bar = bar(0, positive_accum(1), 'FaceColor', 'b', 'FaceAlpha', 0.5);
            yyaxis right
            plot(1:6, recall_rates*100, 'LineWidth', 1.5);   
            yyaxis left
            sum_collect_accum_no_mm = sum(collect_accum_no_mm, 2);
            % Put annotations of the binned recall rate
            for j=0:1
                binned_recall = uint8(sum(collect_accum_no_mm(1+j*3:3+j*3, 1))/sum(sum_collect_accum_no_mm(1+j*3:3+j*3))*100);
                curr_lim = ylim;
                bin_label_yloc = max(sum_collect_accum_no_mm(1+j*3:3+j*3))+20;
                axes_colors = get(gcf, 'defaultAxesColorOrder');
                drawbrace([1+j*3, bin_label_yloc], [3+j*3, bin_label_yloc], 7, 'Color', axes_colors(2,:));
                text(2+j*3, bin_label_yloc+40*(curr_lim(2)/500), strcat(num2str(binned_recall), '%'), 'Color', axes_colors(2,:));
            end
            title({'JAABA bout classification count', ...
                strcat('behaviour:', erase(behav_list{i},'_')), ...
                sprintf('recall %d%%, precision %d%%', uint8(recalls.(behav_list{i})*100), uint8(precisions.(behav_list{i})*100))});
            xlabel('Human annotation combined score \in [1,6]');
            yyaxis left
            ylabel('JAABA bout classification count');
            curr_lim = ylim;
            ylim([curr_lim(1), curr_lim(2)*1.2]);
            yyaxis right
            ylabel('recall in %');
            ylim([0, 100]);
            curr_lim = xlim;
            xlim([curr_lim(1)+0.3, curr_lim(2)-0.3])
            set(bar_plt, {'FaceColor'}, {'b'; 'm'});
            legend({'true or false positives', 'false negatives'}, 'Location', 'northwest');
            hold off;
            set(gcf,'renderer','Painters');
            saveas(gcf, strcat('jaaba_bout_classification_count-',erase(behav_list{i},'_'), '-bar-thres', num2str(thres), '.', erase(img_format, 'c')), img_format);
        end
    end

    save(strcat('bout_matches_', flymat_id, '-thres', num2str(thres), '.mat'), 'bout_matches_all');
end