load('common-params-annot-analysis.mat');
[~, annot_file_name, ~] = fileparts(annot_file);
annot_file_name_elements = strsplit(annot_file_name, '-');
annot_file_id = annot_file_name_elements{end};
load(sprintf('FLYMAT_HumanAnnotation_%s.mat', annot_file_id), 'flymatHumanAnnot');
load('frame_matches_ALL_All_thres0.1.mat', 'frame_matches_all');
construct_from_frame_matches = true;

struct_init_args = [behav_list; squeeze(mat2cell(zeros(4,4,3),4,4,[1,1,1]))'];
human_annot_dist_by_frame = struct(struct_init_args{:});

[score_a, score_b] = ind2sub([4, 4], 1:16);
score_a = score_a - 1;
score_b = score_b - 1;
score_pair = [score_a; score_b]';
score_pair_cell = num2cell(score_pair, 2);

for i=1:length(behav_list)
    if construct_from_frame_matches % Construct from frame matches, which is built on bout matches
        frame_matches = frame_matches_all.(behav_list{i});
        for j=1:length(frame_matches)
            score_pair = frame_matches(j).annot_score_pair;
            score_pair(isnan(score_pair)) = 0;
            if any(score_pair)
                score_pair = score_pair + 1;
                human_annot_dist_by_frame.(behav_list{i})(score_pair(1), score_pair(2)) = ...
                    human_annot_dist_by_frame.(behav_list{i})(score_pair(1), score_pair(2)) + 1;
            else % For human annotation, we do not count FP frames
                continue; 
            end
        end
    else % Construct from original FlymatHumanAnnot
        for j=1:length(flymatHumanAnnot)
            t0s_field = strcat(behav_shorthands{i}, '_t0s');
            t1s_field = strcat(behav_shorthands{i}, '_t1s');
            scores_field = strcat(behav_shorthands{i}, '_scores');
            [~, sparse_score_mat] = convert_bout_annot_into_frame_annot(flymatHumanAnnot(j).(t0s_field), ...
                flymatHumanAnnot(j).(t1s_field), flymatHumanAnnot(j).(scores_field));
            [~, labelled_frame] = find(sparse_score_mat);
            labelled_frame = unique(labelled_frame);
            if ~isempty(labelled_frame)
                for k=1:length(labelled_frame)
                    score_pair = full(sparse_score_mat(:, labelled_frame(k))) + 1;
                    human_annot_dist_by_frame.(behav_list{i})(score_pair(1), score_pair(2)) = ...
                        human_annot_dist_by_frame.(behav_list{i})(score_pair(1), score_pair(2)) + 1;
                end
            end
        end
    end
end

for i=1:length(behav_shorthands)
    human_annot_dist_by_frame_percent.(behav_list{i}) = ...
                    human_annot_dist_by_frame.(behav_list{i})./sum(sum(human_annot_dist_by_frame.(behav_list{i})));   
end

text_label_x_offset = -0.2;
for i=1:length(behav_shorthands)
    trigu_scores = triu(human_annot_dist_by_frame.(behav_list{i})) + ...
        triu(human_annot_dist_by_frame.(behav_list{i})', 1);
    trigu_scores_percent = triu(human_annot_dist_by_frame_percent.(behav_list{i})) + ...
        triu(human_annot_dist_by_frame_percent.(behav_list{i})', 1);
    figure(i);
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
    title({'Frame-wise human annotation distribution (annotator info. ignored)', strcat('behaviour: ', erase(behav_list{i}, '_'))});
    xticks(1:4);
    xticklabels(0:3);
    yticks(1:4);
    yticklabels(0:3);
    set(gcf,'renderer','Painters');
    saveas(gcf, strcat('human_annotation_dist_frame_wise-order_ignored-', erase(behav_list{i}, '_'), '.eps'), 'epsc');
    saveas(gcf, strcat('human_annotation_dist_frame_wise-order_ignored-', erase(behav_list{i}, '_'), '.png'));
end

% Save the data into more reader-friendly format
for i=1:length(behav_shorthands)
    human_annot_dist_by_frame.(behav_list{i}) = struct('score_pair', score_pair_cell, 'count', num2cell(human_annot_dist_by_frame.(behav_list{i})(:)));
end

save('human_annot_dist_by_frame.mat', 'human_annot_dist_by_frame');

