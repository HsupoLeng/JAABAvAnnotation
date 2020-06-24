% Analyze correlation between JAABA confidence score and human annotation
% combined confidence score (by frame)
% Output the violin plot
function score_post_analysis_struct = analyze_human_jaaba_annot_corr_frame_wise(bout_match_sel_str, behav_sel, override_jab_list, is_noRel, frame_matches_mat, plot_bar, plot_violin, img_formats)
    load('common-params-annot-analysis.mat', ...
        'annot_file', 'behav_list', 'behav_shorthands', 'jab_list', 'groundtruth_movie_list');
    load(strcat('bout_matches_', bout_match_sel_str, '.mat'), 'bout_matches_all');
    [~, annot_file_name, ~] = fileparts(annot_file);
    annot_file_name_elements = strsplit(annot_file_name, '-');
    annot_file_id = annot_file_name_elements{end};
    load(sprintf('FLYMAT_HumanAnnotation_%s.mat', annot_file_id), 'flymatHumanAnnot');
     
    movie_name_pairs = cellfun(@(c) strsplit(c, '\'), groundtruth_movie_list, 'UniformOutput', false);
    movie_name_pairs = [movie_name_pairs{:}];
    
    frame_matches_all = struct;

    if is_noRel
        behav_list = {'LungeNewNoRel', 'WingExtNoRel', 'HeadbuttNew'};
    end
    
    if ~isempty(override_jab_list)
        for i=1:length(override_jab_list)
            jab_list{override_jab_list(i).behav_sel} = override_jab_list(i).jab_name;
        end
    end
    
    if ~strcmp(frame_matches_mat, '')
        load(frame_matches_mat, 'frame_matches_all');
    else
        for i=1:length(behav_list)
            if ~ismember(i, behav_sel)
                continue; 
            end
            t0s_field = strcat(behav_shorthands{i}, '_t0s');
            t1s_field = strcat(behav_shorthands{i}, '_t1s');
            scores_field = strcat(behav_shorthands{i}, '_scores');

            bout_matches = bout_matches_all.(behav_list{i});

            frame_mat_fields = {'movie', 'fly', 'frame_num', ...
                'jaaba_score', 'annot_score', 'annot_score_pair'...
                'is_false_negative'};
            init_frame_mat_args = [frame_mat_fields; cell(1,length(frame_mat_fields))];
            frame_matches = struct(init_frame_mat_args{:});
            frame_matches(1) = [];

            frame_count = 0;
            curr_fly_idx = 0;

            for j=1:length(bout_matches)
                bout_match = bout_matches(j);
                fly_idx = find_fly_in_flymat(flymatHumanAnnot, bout_match.movie, bout_match.fly, false);
                if fly_idx ~= curr_fly_idx
                    [~, per_frame_annotation] = convert_bout_annot_into_frame_annot(...
                        flymatHumanAnnot(fly_idx).(t0s_field), flymatHumanAnnot(fly_idx).(t1s_field), ...
                        flymatHumanAnnot(fly_idx).(scores_field));
                    
                    movie_idx = find(contains(groundtruth_movie_list, bout_match.movie));
                    load(fullfile(groundtruth_movie_list{movie_idx}, ...
                        strcat(movie_name_pairs{2*movie_idx},'_JAABA'), ...
                        strcat('scores_', erase(behav_list{i}, {'_','s'}),'.mat')), 'allScores');
                    curr_fly_idx = fly_idx; 
                end
                scores = allScores.scores{bout_match.fly}./allScores.scoreNorm;

                if bout_match.virtual_jaaba_match % false negative
                    annot_union = bout_match.annot_union_start:(bout_match.annot_union_end-1); 
                    annot_intsct_all = arrayfun(@(s,e) s:(e-1), bout_match.annot_intsct_start,bout_match.annot_intsct_end, 'UniformOutput', false);
                    annot_intsct = horzcat(annot_intsct_all{:});
                    stray_fp_frames = [];
                    stray_fn_frames = annot_union; 
                elseif bout_match.annot_score == 0 % false positive
                    annot_union = bout_match.jaaba_bout_start:(bout_match.jaaba_bout_end-1); 
                    annot_intsct = [];
                    stray_fp_frames = annot_union; 
                    stray_fn_frames = [];
                else
                    annot_union = bout_match.annot_union_start:(bout_match.annot_union_end-1);
                    annot_intsct_all = arrayfun(@(s,e) s:(e-1), bout_match.annot_intsct_start, bout_match.annot_intsct_end, 'UniformOutput', false);
                    annot_intsct = horzcat(annot_intsct_all{:});
                    jaaba_all = arrayfun(@(s,e) s:(e-1), bout_match.jaaba_bout_start,bout_match.jaaba_bout_end, 'UniformOutput', false);
                    jaaba_all = horzcat(jaaba_all{:});
                    stray_fp_frames = setdiff(jaaba_all, annot_union);
                    stray_fn_frames = setdiff(annot_union, jaaba_all);                
                end

                frame_slice = unique([annot_union, stray_fp_frames]);
                for k=1:length(frame_slice)
                    frame_match = struct(init_frame_mat_args{:});
                    frame_match.movie = bout_match.movie;
                    frame_match.fly = bout_match.fly; 
                    frame_match.frame_num = frame_slice(k);
                    frame_match.is_false_negative = ismember(frame_slice(k), stray_fn_frames); 
                    frame_match.jaaba_score = scores(frame_slice(k));
                    if ~ismember(frame_slice(k),annot_intsct)
                        if ismember(frame_slice(k), stray_fp_frames)
                            frame_match.annot_score = 0;
                            frame_match.annot_score_pair = [0, 0];
                        else
                            frame_match.annot_score = sum(per_frame_annotation(:, frame_slice(k)));
                            frame_match.annot_score_pair = per_frame_annotation(:, frame_slice(k));
                        end
                    else
                        if length(bout_match.annot_score) ~= length(annot_intsct_all)
                            fprintf('Length of combined annotation score does not match with number of intersections for this bout match\n');
                        end
                        frame_match.annot_score = bout_match.annot_score(cellfun(@(intsct) ismember(frame_slice(k), intsct), annot_intsct_all));
                        frame_match.annot_score_pair = bout_match.annot_score_pair(cellfun(@(intsct) ismember(frame_slice(k), intsct), annot_intsct_all), :);
                    end

                    % If an entry for the frame already exists (due to multi-match), 
                    % and the new incoming entry is neither a false positive nor
                    % a false negative, 
                    % we overwrite that entry record
                    search_range = min(length(frame_matches)-1, 1000); 
                    existing_entry_idx = find(bitand(bitand(strcmp({frame_matches(end-search_range:end).movie}, frame_match.movie), ...
                        [frame_matches(end-search_range:end).fly] == frame_match.fly), [frame_matches(end-search_range:end).frame_num] == frame_match.frame_num));
                    if existing_entry_idx
                        if ~any([frame_match.is_false_negative, ~frame_match.annot_score]) % Takes care of JAABA multi-match
                            existing_entry_idx = existing_entry_idx + length(frame_matches) - search_range - 1; 
                            frame_matches(existing_entry_idx) = frame_match;
                        end
                    else
                        frame_matches(frame_count+1) = frame_match; 
                        frame_count = frame_count + 1;
                    end
                end
                
                if ~mod(j, 200)
                    fprintf('Now at bout match No. %d\n', j);
                end
            end

            frame_matches_all.(behav_list{i}) = frame_matches;     
        end
        save(strcat('frame_matches_', bout_match_sel_str, '.mat'), 'frame_matches_all');
    end
    
    human_annot_scores = 0:6;
    xtick_label_cell_arr = strtrim(cellstr(num2str(human_annot_scores')));
    xtick_label_cell_arr{end+1} = 'training';
    score_post_analysis_struct = struct('annot_score', [], 'jaaba_score', [], 'false_negat_idxs', []);
    score_post_analysis_struct(1) = [];
    jaaba_score_quantile_struct_all = struct();
    for i=1:length(behav_list)
        if ~ismember(i, behav_sel)
            continue; 
        end
        
        frame_matches = frame_matches_all.(behav_list{i});
        % Produce the figures
        score_benchmrk_train = find_train_frame_jaaba_scores(...
            behav_list{i}, jab_list{i});
        jaaba_score = [frame_matches(:).jaaba_score];
        annot_score = uint8(full([frame_matches(:).annot_score]));
        false_negat_mask = [frame_matches(:).is_false_negative];
        positive_train_jaaba_score = score_benchmrk_train.behav;
        negative_train_jaaba_score = score_benchmrk_train.none_behav;
        virt_jaaba_score = jaaba_score(false_negat_mask);
        true_jaaba_score = jaaba_score(~false_negat_mask);
        virt_annot_score = annot_score(false_negat_mask); 
        true_annot_score = annot_score(~false_negat_mask);
        
        score_post_analysis_struct(find(behav_sel == i)).annot_score = annot_score; 
        score_post_analysis_struct(find(behav_sel == i)).jaaba_score = jaaba_score; 
        score_post_analysis_struct(find(behav_sel == i)).false_negat_idxs = find(false_negat_mask);
        
        false_negat_accum = hist(virt_annot_score, human_annot_scores);
        positive_accum = hist(true_annot_score, human_annot_scores);
        collect_accum = [positive_accum(2:end)', false_negat_accum(2:end)'];
        precisions.(behav_list{i}) = sum(positive_accum(2:end))/sum(positive_accum);
        recalls.(behav_list{i}) = sum(positive_accum(2:end))/(sum(false_negat_accum(2:end))+sum(positive_accum(2:end)));
        recall_rates = collect_accum(:,1)./sum(collect_accum, 2);
        
        test_violin_norm = max([false_negat_accum, positive_accum]);
        positive_train_size = length(score_benchmrk_train.behav);
        negat_train_size = length(score_benchmrk_train.none_behav);
        train_violin_norm = max([positive_train_size, negat_train_size]);
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
            
            jaaba_score_quantile_struct = struct();
            for m=human_annot_scores
                jaaba_score_quantile_struct(m+1).data_type = 'test';
                jaaba_score_quantile_struct(m+1).annot_score = m;
                scores_one_cat = true_jaaba_score(true_annot_score == m);
                if ~isempty(scores_one_cat)
                    jaaba_score_quantile_struct(m+1).jaaba_posit_quantiles = quantile(scores_one_cat, [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]);
                else
                    jaaba_score_quantile_struct(m+1).jaaba_posit_quantiles = nan(5, 1);
                end
                scores_one_cat = virt_jaaba_score(virt_annot_score == m);
                if ~isempty(scores_one_cat)
                    jaaba_score_quantile_struct(m+1).jaaba_negat_quantiles = quantile(scores_one_cat, [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]);
                else
                    jaaba_score_quantile_struct(m+1).jaaba_negat_quantiles = nan(5, 1);
                end
                scores_one_cat = horzcat(true_jaaba_score(true_annot_score == m), virt_jaaba_score(virt_annot_score == m));
                if ~isempty(scores_one_cat)
                    jaaba_score_quantile_struct(m+1).jaaba_combined_quantiles = quantile(scores_one_cat, [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]);
                else
                    jaaba_score_quantile_struct(m+1).jaaba_combined_quantiles = nan(5, 1);
                end

            end
            jaaba_score_quantile_struct(8).data_type = 'train';
            jaaba_score_quantile_struct(8).annot_score = nan;
            jaaba_score_quantile_struct(8).jaaba_posit_quantiles = quantile(positive_train_jaaba_score, [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]);
            jaaba_score_quantile_struct(8).jaaba_negat_quantiles = quantile(negative_train_jaaba_score, [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]);
            jaaba_score_quantile_struct(8).jaaba_combined_quantiles = quantile(horzcat(positive_train_jaaba_score, negative_train_jaaba_score), [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]);
            jaaba_score_quantile_struct_all.(behav_list{i}) = jaaba_score_quantile_struct; 


            all_plotted_scores = [true_jaaba_score, virt_jaaba_score, positive_train_jaaba_score, negative_train_jaaba_score];
            axis([-0.5+offset_hack, 7.5+offset_hack, min(all_plotted_scores)-1.2, ...
                max(all_plotted_scores)+0.3]);
            line([-0.5, 6.5]+offset_hack, repmat(mean(score_benchmrk_train.behav), [1, 2]), ...
                'Color', 'k', 'LineStyle', '--', 'DisplayName', 'avg. over all positive{\it training} bouts');
            line([-0.5, 6.5]+offset_hack, repmat(mean(score_benchmrk_train.none_behav), [1, 2]), ...
                'Color', 'k', 'LineStyle', '-.', 'DisplayName', 'avg. over all negative{\it training} bouts');
            ylim_tmp = ylim; 
            train_bkg = area([6.5+offset_hack, 7.5+offset_hack], [ylim_tmp(2),ylim_tmp(2)], ylim_tmp(1), 'FaceColor', [0.3, 0.3, 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            train_bkg.Annotation.LegendInformation.IconDisplayStyle = 'off';
            title({strcat('JAABA score vs. human combined score for test frames'), ...
                strcat('behaviour:', erase(behav_list{i},'_'))});
            xlabel('Human annotation combined score \in [0,6]');
            ylabel(strcat('JAABA score (for a frame)'));
            set(gca, 'xtick', [human_annot_scores, 7]++offset_hack);
            set(gca, 'xticklabel', xtick_label_cell_arr);
            hold off;
            set(gcf,'renderer','Painters');
            for j=1:length(img_formats)
                if contains(img_formats{j}, 'eps')
                    img_format = 'epsc';
                else
                    img_format = 'png';
                end
                saveas(gcf, strcat(strcat('jaaba_score_perframe_normed'), ...
                    '_vs_annot_perframe-',erase(behav_list{i},'_'), '-violin-', bout_match_sel_str, '.', erase(img_format, 'c')), img_format);
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
            sum_collect_accum = sum(collect_accum, 2);
            % Put annotations of the binned recall rate
            for j=0:1
                binned_recall = uint8(sum(collect_accum(1+j*3:3+j*3, 1))/sum(sum_collect_accum(1+j*3:3+j*3))*100);
                curr_lim = ylim;
                bin_label_yloc = max(sum_collect_accum(1+j*3:3+j*3))+20;
                axes_colors = get(gcf, 'defaultAxesColorOrder');
                drawbrace([1+j*3, bin_label_yloc], [3+j*3, bin_label_yloc], 7, 'Color', axes_colors(2,:));
                text(2+j*3, bin_label_yloc+40*(curr_lim(2)/500), strcat(num2str(binned_recall), '%'), 'Color', axes_colors(2,:));
            end
            title({'JAABA frame classification count', ...
                strcat('behaviour:', erase(behav_list{i},'_')), ...
                sprintf('recall %d%%, precision %d%%', uint8(recalls.(behav_list{i})*100), uint8(precisions.(behav_list{i})*100))});
            xlabel('Human annotation combined score \in [1,6]');
            yyaxis left
            ylabel('JAABA frame classification count');
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
            for j=1:length(img_formats)
                if contains(img_formats{j}, 'eps')
                    img_format = 'epsc';
                else
                    img_format = 'png';
                end
                saveas(gcf, strcat('jaaba_frame_classification_count-',erase(behav_list{i},'_'), '-bar-', bout_match_sel_str, '.', erase(img_format, 'c')), img_format);
            end
        end
    end
    save(sprintf('jaaba_score_quantile_struct_all-%s_frame_wise.mat', bout_match_sel_str), 'jaaba_score_quantile_struct_all');
end