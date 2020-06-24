% Compile and draw figures for distribution of human annotation confidence
% (by bout)
load('common-params-annot-analysis.mat');
load('bout_matches_ALL_All_thres0.1.mat');

struct_init_args = [behav_list; squeeze(mat2cell(zeros(4,4,3),4,4,[1,1,1]))'];
human_annot_dist_collective = struct(struct_init_args{:});

annotators_by_behav = cell(length(behav_list), 1);
for i=1:length(behav_shorthands)
    annot_score_all = cellfun(@(scores) max(scores), {bout_matches_all.(behav_list{i}).annot_score});
    annotators_by_behav{i} = unique([bout_matches_all.(behav_list{i})(annot_score_all > 0).annotators]);
end
human_annot_dist_by_annotator = struct;
for i=1:length(behav_shorthands)
    human_annot_dist_by_annotator.(behav_list{i}) = struct; 
    for j=1:length(annotators_by_behav{i})
        human_annot_dist_by_annotator.(behav_list{i}).(annotators_by_behav{i}{j}) = zeros(4, 4);
    end
end

[score_a, score_b] = ind2sub([4, 4], 1:16);
score_a = score_a - 1;
score_b = score_b - 1;
score_pair = [score_a; score_b]';
score_pair_cell = num2cell(score_pair, 2);

for i=1:length(behav_shorthands)
    for k=1:length(bout_matches_all.(behav_list{i}))
        annot_score_all = bout_matches_all.(behav_list{i})(k).annot_score;
        if ~(annot_score_all)
            continue;
        end
        [first_annotator, second_annotator] = bout_matches_all.(behav_list{i})(k).annotators{:};
        scores = bout_matches_all.(behav_list{i})(k).annot_score_pair;
        scores(isnan(scores)) = 0;

        for l=1:size(scores, 1)
            score_pair = scores(l, :) + 1;
            human_annot_dist_collective.(behav_list{i})(score_pair(1), score_pair(2)) = ...
                human_annot_dist_collective.(behav_list{i})(score_pair(1), score_pair(2)) + 1;

            human_annot_dist_by_annotator.(behav_list{i}).(first_annotator)(score_pair(2), score_pair(1)) = ...
                human_annot_dist_by_annotator.(behav_list{i}).(first_annotator)(score_pair(2), score_pair(1)) + 1;
            human_annot_dist_by_annotator.(behav_list{i}).(second_annotator)(score_pair(1), score_pair(2)) = ...
                human_annot_dist_by_annotator.(behav_list{i}).(second_annotator)(score_pair(1), score_pair(2)) + 1;
        end
    end
end
    
for i=1:length(behav_shorthands)
    human_annot_dist_collective_percent.(behav_list{i}) = ...
                    human_annot_dist_collective.(behav_list{i})./sum(sum(human_annot_dist_collective.(behav_list{i})));   
end

for i=1:length(behav_shorthands)
    for j=1:length(annotators_by_behav{i})
        human_annot_dist_by_annotator_percent.(behav_list{i}).(annotators_by_behav{i}{j}) = ...
                        human_annot_dist_by_annotator.(behav_list{i}).(annotators_by_behav{i}{j})/...
                        sum(sum(human_annot_dist_by_annotator.(behav_list{i}).(annotators_by_behav{i}{j})));   
    end
end

text_label_x_offset = -0.2;

% ----- Generate plots for the collective analysis ----
% Simply plot the array as is
for i=1:length(behav_shorthands)
    figure(i);
    hold on;
    imagesc(human_annot_dist_collective_percent.(behav_list{i}));
    set(gca,'YDir','normal');
    colormap parula;
    cb = colorbar;
    cb.TickLabels = strcat(cellstr(num2str(round(cb.Ticks(:)*100, 1))), '%');
    dont_care_gray = [0.4, 0.4, 0.4];
    dont_care_line_color = 'k';
    patch([0.5, 0.5, 1.5, 1.5], [1.5, 0.5, 0.5, 1.5], dont_care_gray, 'EdgeColor', 'none');
    plot(0.5:1.5, 0.5:1.5, dont_care_line_color, 0.5:1.5, 1.5:-1:0.5, dont_care_line_color);
    hold off;
    text_xloc = repelem(1:4, 4)+text_label_x_offset;
    text_yloc = repmat(1:4, 1, 4);
    text_contents = strcat(cellstr(num2str(round(human_annot_dist_collective_percent.(behav_list{i})(1:end)'*100, 1))), '%');
    text_contents = cellfun(@(x,y) {x;strcat('(', num2str(y), ')')}, text_contents, num2cell(human_annot_dist_collective.(behav_list{i})(:), 2), 'UniformOutput', false);
    text(text_xloc(2:end), text_yloc(2:end), text_contents(2:end));
    title({'Human annotation distribution', strcat('behaviour: ', erase(behav_list{i}, '_'))});
    xticks(1:4);
    xticklabels(0:3);
    yticks(1:4);
    yticklabels(0:3);
    set(gcf,'renderer','Painters');
    saveas(gcf, strcat('human_annotation_dist-order_preserved-', erase(behav_list{i}, '_'), '.eps'), 'epsc');
    saveas(gcf, strcat('human_annotation_dist-order_preserved-', erase(behav_list{i}, '_'), '.png'));
end

% Ignore the annotator axis 
% by adding the lower triangle of the plot to cooresponding entry in the
% upper triangle
for i=1:length(behav_shorthands)
    trigu_scores = triu(human_annot_dist_collective.(behav_list{i})) + ...
        triu(human_annot_dist_collective.(behav_list{i})', 1);
    trigu_scores_percent = triu(human_annot_dist_collective_percent.(behav_list{i})) + ...
        triu(human_annot_dist_collective_percent.(behav_list{i})', 1);
    figure(i+3);
    hold on;
    imagesc(trigu_scores_percent);
    dont_care_gray = [0.4, 0.4, 0.4];
    dont_care_line_color = 'k';
    patch([0.5, 0.5, 1.5, 1.5], [4.5, 0.5, 0.5, 4.5], dont_care_gray, 'EdgeColor', 'none');
    patch([1.5, 1.5, 2.5, 2.5], [4.5, 2.5, 2.5, 4.5], dont_care_gray, 'EdgeColor', 'none');
    patch([2.5, 2.5, 3.5, 3.5], [4.5, 3.5, 3.5, 4.5], dont_care_gray, 'EdgeColor', 'none');
    plot(0.5:1.5, 0.5:1.5, dont_care_line_color, 0.5:1.5, 1.5:-1:0.5, dont_care_line_color);
    plot(0.5:1:3.5, 1.5:1:4.5, dont_care_line_color, 0.5:1:2.5, 2.5:1:4.5, dont_care_line_color, 0.5:1.5, 3.5:4.5, dont_care_line_color);
    plot(0.5:1.5, 2.5:-1:1.5, dont_care_line_color, 0.5:1.5, 3.5:-1:2.5, dont_care_line_color, 0.5:2.5, 4.5:-1:2.5, dont_care_line_color, 1.5:2.5, 4.5:-1:3.5, dont_care_line_color, 2.5:3.5, 4.5:-1:3.5, dont_care_line_color);
    plot(0.5:1.5, [1.5, 1.5],dont_care_line_color, 0.5:2.5, [2.5, 2.5, 2.5], dont_care_line_color, 0.5:3.5, [3.5, 3.5, 3.5, 3.5], dont_care_line_color);
    plot([1.5, 1.5, 1.5, 1.5, 1.5], 0.5:4.5, dont_care_line_color, [2.5, 2.5, 2.5], 2.5:4.5, dont_care_line_color, [3.5, 3.5], 3.5:4.5, dont_care_line_color);
    hold off;
    colormap parula;
    cb = colorbar;
    cb.TickLabels = strcat(cellstr(num2str(round(cb.Ticks(:)*100, 1))), '%');
    text_contents = strcat(cellstr(num2str(round(trigu_scores_percent(trigu_scores_percent~=0)*100, 1))), '%');
    text_contents = cellfun(@(x,y) {x;strcat('(', num2str(y), ')')}, text_contents, num2cell(trigu_scores(trigu_scores~=0), 2), 'UniformOutput', false);
    text([2,2,3,3,3,4,4,4,4]+text_label_x_offset, [1,2,1,2,3,1,2,3,4], text_contents);
    title({'Human annotation distribution (annotator info. ignored)', strcat('behaviour: ', erase(behav_list{i}, '_'))});
    xticks(1:4);
    xticklabels(0:3);
    yticks(1:4);
    yticklabels(0:3);
    set(gcf,'renderer','Painters');
    saveas(gcf, strcat('human_annotation_dist-order_ignored-', erase(behav_list{i}, '_'), '.eps'), 'epsc');
    saveas(gcf, strcat('human_annotation_dist-order_ignored-', erase(behav_list{i}, '_'), '.png'));
end

% Save the data into more reader-friendly format
for i=1:length(behav_shorthands)
    human_annot_dist_collective.(behav_list{i}) = struct('score_pair', score_pair_cell, 'count', num2cell(human_annot_dist_collective.(behav_list{i})(:)));
    human_annot_dist_by_annotator.(behav_list{i}) = structfun(@(arr) struct('score_pair_otherVSself', score_pair_cell, 'count', num2cell(arr(:))), human_annot_dist_by_annotator.(behav_list{i}), 'UniformOutput', false);
end

save('human_annot_dist.mat', 'human_annot_dist_collective', 'human_annot_dist_by_annotator');
