function search_grid_performance_bout_wise(file_id, behav_sel, threshold, is_NoRel, maxgap, minbout, override_jab_list, ourchoice_idx)
    load('common-params-annot-analysis.mat', 'behav_list', 'behav_shorthands');
    human_annot_scores = 0:6;
    if is_NoRel
        behav_list = {'LungeNewNoRel', 'WingExtNoRel', 'HeadbuttNew'};
    end
    
    MVJ = zeros(length(threshold)*length(maxgap)*length(minbout), 7);
    % For each combination of confidence threshold, min bout length and max
    % gap between bouts to merge, generate FLYMAT, generate bout_matches
    % between human annotated bouts and JAABA bouts, finally, compute
    % precision and recall based on bouts
    for i=1:length(threshold)
        for j=1:length(maxgap)
            for k=1:length(minbout)
                maxgap_allbehav = struct();
                minbout_allbehav = struct();
                for l=1:length(behav_shorthands)
                    maxgap_allbehav.(behav_shorthands{l}) = 0;
                    minbout_allbehav.(behav_shorthands{l}) = 0;
                end
                maxgap_allbehav.(behav_shorthands{behav_sel}) = maxgap(j);
                minbout_allbehav.(behav_shorthands{behav_sel}) = minbout(k);
                % The following FLYMAT is only correct for the selected
                % behavior. 
                OrgData020816_XuboCopy('D:\xubo\NewTrainingFiles-OriginalCopy', ...
                    'TrainingSamples_genotype3.xlsx', threshold(i), is_NoRel, maxgap_allbehav, minbout_allbehav);
                analyze_human_jaaba_annot_corr_v3('temp', behav_sel, override_jab_list, is_NoRel, '', false, false, false, '');
                
                load(strcat('bout_matches_temp.mat'), 'bout_matches_all', 'bout_masks');
                
                % The following section uses the conventional way of
                % defining precision and recall, i.e. 
                % precision = TP/(TP + FP); recall = (TP/(TP + FN)). 
                % It does not distinguish between high-confidence bouts and
                % low-confidence bouts. 
                % One twist: for recall, since one JAABA bout might be
                % matched to multiple human annotated bouts, we do not
                % count more than one bout match for each JAABA bout. 
                % Another twist: all human annotated bouts that overlaps
                % are counted as one single annotated bout with score being
                % the maximum combined confidence score
                annot_score = cellfun(@(scores) max(scores), {bout_matches_all.(behav_list{behav_sel})(bout_masks{behav_sel}).annot_score});
                false_negat_idxs = [bout_matches_all.(behav_list{behav_sel})(bout_masks{behav_sel}).virtual_jaaba_match];
                multi_match_idxs = [bout_matches_all.(behav_list{behav_sel})(bout_masks{behav_sel}).multi_match];
                multi_match_idxs(isnan(multi_match_idxs)) = 0;
                virt_annot_score_no_mm = annot_score(bitand(false_negat_idxs, ~multi_match_idxs)); 
                true_annot_score_no_mm = annot_score(bitand(~false_negat_idxs, ~multi_match_idxs));
                false_negat_accum_no_mm = hist(virt_annot_score_no_mm, human_annot_scores);
                positive_accum_no_mm = hist(true_annot_score_no_mm, human_annot_scores);
                recall = sum(positive_accum_no_mm(2:end))/(sum(false_negat_accum_no_mm(2:end))+sum(positive_accum_no_mm(2:end)));

                true_annot_score = annot_score(~false_negat_idxs);
                positive_accum = hist(true_annot_score, human_annot_scores);
                precision = sum(positive_accum(2:end))/sum(positive_accum);
                
                config_idx = sub2ind([length(minbout), length(maxgap), length(threshold)], k, j, i);
                MVJ(config_idx, [1,2,3]) = [threshold(i), maxgap(j), minbout(k)];
                MVJ(config_idx, 6) = precision; 
                MVJ(config_idx, 7) = recall; 
            end
        end
    end
    
    % Plot the results
    for i=1:2
        if i == 1
            continue; % We will be using the conventional definition of precision and accuracy
        end
        
        figure();
        % Plot shapes indicating three parameters on three different axes
        conf_ax = axes;
        max_ax = axes;
        min_ax = axes;
        % Generate color values corresponding to the parameter values.
        % Larger the parameter, brighter the color. 
        % Each parameter takes one color channel (r, g, or b)
        confcolor_mat = [zeros(size(MVJ, 1), 1), (MVJ(:, 1)/(0.2*2)) + 0.25, zeros(size(MVJ, 1), 1)];
        maxcolor_mat = [(MVJ(:, 2)/max(MVJ(:,2))), zeros(size(MVJ, 1), 1), zeros(size(MVJ, 1), 1)];
        mincolor_mat = [zeros(size(MVJ, 1), 1), zeros(size(MVJ, 1), 1), (MVJ(:,3)/max(MVJ(:, 3)))];
        
        hold on;
        scatter(conf_ax, MVJ(:, 4+2*(i-1)), MVJ(:, 5+2*(i-1)), 500, confcolor_mat, 's', 'LineWidth', 1.5);
        scatter(max_ax, MVJ(:, 4+2*(i-1)), MVJ(:, 5+2*(i-1)), 200, maxcolor_mat, 'd', 'LineWidth', 1.5);
        scatter(min_ax, MVJ(:, 4+2*(i-1)), MVJ(:, 5+2*(i-1)), 75, mincolor_mat, 'o', 'LineWidth', 1.5);
        linkaxes([conf_ax, max_ax, min_ax]);
        max_ax.Visible = 'off';
        max_ax.XTick = [];
        max_ax.YTick = [];
        min_ax.Visible = 'off';
        min_ax.XTick = [];
        min_ax.YTick = [];
        % One colormap for each parameter
        colormap(conf_ax, unique(confcolor_mat, 'rows'));
        colormap(max_ax, unique(maxcolor_mat, 'rows'));
        colormap(min_ax, unique(mincolor_mat, 'rows'));
        set([conf_ax, max_ax, min_ax], 'Position', [0.10, 0.11, 0.685, 0.815]);
        colorbar(conf_ax, 'Position', [0.79, 0.11, 0.0675, 0.815], 'Ticks', []);
        colorbar(max_ax, 'Position', [0.86, 0.11, 0.0675, 0.815], 'Ticks', []);
        colorbar(min_ax, 'Position', [0.93, 0.11, 0.0675, 0.815], 'Ticks', []);
        % plot the value you chose in the end as a black tick mark
        scatter(min_ax, MVJ(ourchoice_idx, 4+2*(i-1)), MVJ(ourchoice_idx, 5+2*(i-1)), 75,'+','MarkerEdgeColor',[0 0 0],'LineWidth',1.5)
        hold off
        xlabel(conf_ax, 'Precision');
        ylabel(conf_ax, 'Recall');

        if i==1
            saveas(gcf, strcat('ManualvsJAABA_grid_performance_original-', file_id, '.png'), 'png');
            saveas(double(gcf), strcat('ManualvsJAABA_grid_performance_original-', file_id, '.eps'), 'eps');
        else
            saveas(gcf, strcat('ManualvsJAABA_grid_performance_bout_new-', file_id, '.png'), 'png');
            saveas(double(gcf), strcat('ManualvsJAABA_grid_performance_bout_new-', file_id, '.eps'), 'eps');
        end
    end
end