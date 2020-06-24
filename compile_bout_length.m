% Compile bout length histogram for human annotation and JAABA prediction
close all;
load('bout_matches_ALL_all_thres0.1_nominlength.mat', 'bout_matches_all');
load('common-params-annot-analysis.mat');
[~, annot_file_name, ~] = fileparts(annot_file);
annot_file_name_elements = strsplit(annot_file_name, '-');
annot_file_id = annot_file_name_elements{end};
load(sprintf('FLYMAT_HumanAnnotation_%s.mat', annot_file_id), 'flymatHumanAnnot');


human_annot_scores = 0:6; 
wingext_bin_edges = [0,6,30,60:60:300,600, inf]-0.5;
% wingext_bin_edges = (1:16)-0.5;
legend_wingext = arrayfun(@(bin_start, bin_end) sprintf('%d~%d', bin_start, bin_end), wingext_bin_edges(1:end-1), wingext_bin_edges(2:end), 'UniformOutput', false);
annot_intsct_bout_length_binned_histo_cell = cell(length(behav_list), length(human_annot_scores));
annot_union_bout_length_binned_histo_cell = cell(length(behav_list), length(human_annot_scores));
jaaba_bout_length_binned_histo_cell = cell(length(behav_list), length(human_annot_scores));

for i=1:length(behav_list)
    behav = behav_list{i};
    bout_matches = bout_matches_all.(behav);
    
    % Pre-process the data: bin bout lengths according to their scores
    annot_intsct_bout_length_binned = cell(size(human_annot_scores));
    annot_union_bout_length_binned = cell(size(human_annot_scores));
    jaaba_bout_length_binned = cell(size(human_annot_scores));
    
    % Remove matches that have negative annotation intersection lengths
    % and report the number of such occurrences
    bout_mask_length_problem = arrayfun(@(s) min([s.annot_intsct_end - s.annot_intsct_start]), ...
        bout_matches) < 0;
    bout_matches(bout_mask_length_problem) = [];
    if sum(bout_mask_length_problem)
        fprintf('Removing %d bout matches with negative length of bout intersection\n', sum(bout_mask_length_problem));
    end
    for j=1:length(human_annot_scores)
        annot_score = human_annot_scores(j);
        bout_mask = arrayfun(@(s) max(s.annot_score), bout_matches) == annot_score;
        annot_intsct_bout_length_cells = arrayfun(@(s) s.annot_intsct_end - s.annot_intsct_start, ...
            bout_matches(bout_mask), 'UniformOutput', false);
        annot_intsct_bout_length_binned{j} = horzcat(annot_intsct_bout_length_cells{:});
        annot_union_bout_length_cells = arrayfun(@(s) s.annot_union_end - s.annot_union_start, ...
            bout_matches(bout_mask), 'UniformOutput', false);
        annot_union_bout_length_binned{j} = horzcat(annot_union_bout_length_cells{:});
        bout_mask_for_jaaba_score = bitand(bout_mask, ~[bout_matches(:).virtual_jaaba_match]); % Remove false negatives
        jaaba_bout_length_cells = arrayfun(@(s) s.jaaba_bout_end - s.jaaba_bout_start, ...
            bout_matches(bout_mask_for_jaaba_score), 'UniformOutput', false);
        jaaba_bout_length_binned{j} = horzcat(jaaba_bout_length_cells{:});
    end
    
    % Lump everthing together and plot a histogram
    figure();
    if i == 2
        h_obj = histogram(horzcat(annot_intsct_bout_length_binned{:}), 'BinEdges', wingext_bin_edges);
%         xticks(h_obj.BinEdges);
        xticks(h_obj.BinEdges+0.5);
        xticklabels(cellstr(num2str(h_obj.BinEdges'+0.5)));
        xlim([0, 600]);
    else
        h_obj = histogram(horzcat(annot_intsct_bout_length_binned{:}), 'BinMethod', 'integers');
        xticks(h_obj.BinEdges+0.5);
        xticklabels(cellstr(num2str(h_obj.BinEdges'+0.5)));
    end
    title({'Histogram of human annotated bout lengths (intersection)', ...
        sprintf('behaviour: %s', strrep(behav, '_', ''))});
%     xlim([0, inf]);
    set(gcf,'renderer','Painters');
%     saveas(gcf, sprintf('human_annotated_bout_length_intersection_histogram-%s.png', behav), 'png');
%     saveas(gcf, sprintf('human_annotated_bout_length_intersection_histogram-%s.eps', behav), 'epsc');
    
    figure();
    if i == 2
        h_obj = histogram(horzcat(annot_union_bout_length_binned{:}), 'BinEdges', wingext_bin_edges);
%         xticks(h_obj.BinEdges);
        xticks(h_obj.BinEdges+0.5);
        xticklabels(cellstr(num2str(h_obj.BinEdges'+0.5)));
        xlim([0, 600]);
    else
        h_obj = histogram(horzcat(annot_union_bout_length_binned{:}), 'BinMethod', 'integers');
        xticks(h_obj.BinEdges+0.5);
        xticklabels(cellstr(num2str(h_obj.BinEdges'+0.5)));
    end
    title({'Histogram of human annotated bout lengths (union)', ...
        sprintf('behaviour: %s', strrep(behav, '_', ''))});
%     xlim([0, inf]);
    set(gcf,'renderer','Painters');
%     saveas(gcf, sprintf('human_annotated_bout_length_union_histogram-%s.png', behav), 'png');
%     saveas(double(gcf), sprintf('human_annotated_bout_length_union_histogram-%s.eps', behav), 'epsc');

    figure();
    if i == 2
        h_obj = histogram(horzcat(jaaba_bout_length_binned{:}), 'BinEdges', wingext_bin_edges);
        xticks(h_obj.BinEdges+0.5);
        xticklabels(cellstr(num2str(h_obj.BinEdges'+0.5)));
        xlim([0, 600]);
    else
        h_obj = histogram(horzcat(jaaba_bout_length_binned{:}), 'BinMethod', 'integers');
        xticks(h_obj.BinEdges+0.5);
        xticklabels(cellstr(num2str(h_obj.BinEdges'+0.5)));
    end
    title({'Histogram of JAABA bout lengths', ...
        sprintf('behaviour: %s', strrep(behav, '_', ''))});
%     xlim([0, inf]);
    if i ~= 2 
        xticks(h_obj.BinEdges+0.5);
        xticklabels(cellstr(num2str(h_obj.BinEdges'+0.5)));
    end
    set(gcf,'renderer','Painters');
%     saveas(gcf, sprintf('jaaba_bout_length_histogram-%s.png', behav), 'png');
%     saveas(double(gcf), sprintf('jaaba_bout_length_histogram-%s.eps', behav), 'epsc');
    
    % Draw area plots for binned bout lengths, for human annotation
    % intersection, union, and JAABA annotaion, respectively
    % Human annotation intersection
    figure();
    if i == 2
        histcount_cell = cellfun(@(lengths) histcounts(lengths, 'BinEdges', wingext_bin_edges), annot_intsct_bout_length_binned, 'UniformOutput', false);
        histcount_cell_augmented = cellfun(@(h) [wingext_bin_edges(1:end-1)', wingext_bin_edges(2:end)', h'], histcount_cell, 'UniformOutput', false);
    else
        min_length = min(horzcat(annot_intsct_bout_length_binned{:}));
        % max_length = max(horzcat(annot_intsct_bout_length_binned{:}));
        if i == 1
            max_length = 20;
        else
            max_length = 15;
        end
        histcount_cell = cellfun(@(lengths) histcounts(lengths, 'BinMethod', 'integers', 'BinLimits', [min_length, max_length]), annot_intsct_bout_length_binned, 'UniformOutput', false);
        histcount_cell_augmented = cellfun(@(h) [(min_length:max_length)', h'], histcount_cell, 'UniformOutput', false);
    end
    annot_intsct_bout_length_binned_histo_cell(i, :) = histcount_cell_augmented; 
    
    area(vertcat(histcount_cell{:}), 'FaceColor', 'flat');
    xticklabels(human_annot_scores);
    if i == 2
        legend_wingext = arrayfun(@(bin_start, bin_end) sprintf('%d~%d', bin_start, bin_end), wingext_bin_edges(1:end-1)+0.5, wingext_bin_edges(2:end)+0.5, 'UniformOutput', false);
        legend(legend_wingext, 'Location', 'northwest');
    else
        legend(cellstr(num2str((min_length:max_length)')), 'Location', 'northwest');
    end
    ylim(2.5*ylim);
    xlabel('Human annotation combined score (in range [0, 6])');
    ylabel('Number of bouts');
    title({'Human annotated bout length (intersection) distribution', 'by human annotation combined score', ...
        sprintf('behaviour: %s', strrep(behav, '_', ''))});
    set(gcf,'renderer','Painters');
%     saveas(gcf, sprintf('human_annotated_bout_length_intersection_disttribution-%s.png', behav), 'png');
%     saveas(double(gcf), sprintf('human_annotated_bout_length_intersection_disttribution-%s.eps', behav), 'epsc');
    
    % Human annotation union
    figure();
    if i == 2
        histcount_cell = cellfun(@(lengths) histcounts(lengths, 'BinEdges', wingext_bin_edges), annot_union_bout_length_binned, 'UniformOutput', false);
        histcount_cell_augmented = cellfun(@(h) [wingext_bin_edges(1:end-1)', wingext_bin_edges(2:end)', h'], histcount_cell, 'UniformOutput', false);
    else
        min_length = min(horzcat(annot_union_bout_length_binned{:}));
        % max_length = max(horzcat(annot_intsct_bout_length_binned{:}));
        if i == 1
            max_length = 20;
        else
            max_length = 15;
        end
        histcount_cell = cellfun(@(lengths) histcounts(lengths, 'BinMethod', 'integers', 'BinLimits', [min_length, max_length]), annot_union_bout_length_binned, 'UniformOutput', false);
        histcount_cell_augmented = cellfun(@(h) [(min_length:max_length)', h'], histcount_cell, 'UniformOutput', false);
    end
    annot_union_bout_length_binned_histo_cell(i, :) = histcount_cell_augmented; 
    
    area(vertcat(histcount_cell{:}), 'FaceColor', 'flat');
    xticklabels(human_annot_scores);
    if i == 2
        legend_wingext = arrayfun(@(bin_start, bin_end) sprintf('%d~%d', bin_start, bin_end), wingext_bin_edges(1:end-1)+0.5, wingext_bin_edges(2:end)+0.5, 'UniformOutput', false);
        legend(legend_wingext, 'Location', 'northwest');
    else
        legend(cellstr(num2str((min_length:max_length)')), 'Location', 'northwest');
    end
    ylim(2.5*ylim);
    xlabel('Human annotation combined score (in range [0, 6])');
    ylabel('Number of bouts');
    title({'Human annotated bout length (union) distribution', 'by human annotation combined score', ...
        sprintf('behaviour: %s', strrep(behav, '_', ''))});
    set(gcf,'renderer','Painters');
%     saveas(gcf, sprintf('human_annotated_bout_length_union_disttribution-%s.png', behav), 'png');
%     saveas(double(gcf), sprintf('human_annotated_bout_length_union_disttribution-%s.eps', behav), 'epsc');
    
    % Jaaba annotation
    figure();
    if i == 2
        histcount_cell = cellfun(@(lengths) histcounts(lengths, 'BinEdges', wingext_bin_edges), jaaba_bout_length_binned, 'UniformOutput', false);
        histcount_cell_augmented = cellfun(@(h) [wingext_bin_edges(1:end-1)', wingext_bin_edges(2:end)', h'], histcount_cell, 'UniformOutput', false);
    else
        min_length = min(horzcat(jaaba_bout_length_binned{:}));
        % max_length = max(horzcat(annot_intsct_bout_length_binned{:}));
        if i == 1
            max_length = 20;
        else
            max_length = 15;
        end
        histcount_cell = cellfun(@(lengths) histcounts(lengths, 'BinMethod', 'integers', 'BinLimits', [min_length, max_length]), jaaba_bout_length_binned, 'UniformOutput', false);
        histcount_cell_augmented = cellfun(@(h) [(min_length:max_length)', h'], histcount_cell, 'UniformOutput', false);
    end
    jaaba_bout_length_binned_histo_cell(i, :) = histcount_cell_augmented; 
    
    area(vertcat(histcount_cell{:}), 'FaceColor', 'flat');
    xticklabels(human_annot_scores);
    if i == 2
        legend_wingext = arrayfun(@(bin_start, bin_end) sprintf('%d~%d', bin_start, bin_end), wingext_bin_edges(1:end-1)+0.5, wingext_bin_edges(2:end)+0.5, 'UniformOutput', false);
        legend(legend_wingext, 'Location', 'northwest');
    else
        legend(cellstr(num2str((min_length:max_length)')), 'Location', 'northwest');
    end
    ylim(2.5*ylim);
    xlabel('Human annotation combined score (in range [0, 6])');
    ylabel('Number of bouts');
    title({'JAABA bout length distribution', 'by human annotation combined score', ...
        sprintf('behaviour: %s', strrep(behav, '_', ''))});
    set(gcf,'renderer','Painters');
%     saveas(gcf, sprintf('jaaba_annotated_bout_length_disttribution-%s.png', behav), 'png');
%     saveas(double(gcf), sprintf('jaaba_annotated_bout_length_disttribution-%s.eps', behav), 'epsc');
end
