load('common-params-annot-analysis.mat');
load('FLYMAT_HumanAnnotation_v3.mat');

struct_init_args = [behav_list; squeeze(mat2cell(zeros(4,4,3),4,4,[1,1,1]))'];
human_annot_dist_by_frame = struct(struct_init_args{:});


for i=1:length(behav_list)
    count_total = 0;
    for j=1:length(flymatHumanAnnot)
        t0s_field = strcat(behav_shorthands{i}, '_t0s');
        t1s_field = strcat(behav_shorthands{i}, '_t1s');
        scores_field = strcat(behav_shorthands{i}, '_scores');
        [count, sparse_score_mat] = convert_bout_annot_into_frame_annot(flymatHumanAnnot(j).(t0s_field), ...
            flymatHumanAnnot(j).(t1s_field), flymatHumanAnnot(j).(scores_field));
        count_total = [count_total, count];
        labelled_inds = find(sparse_score_mat);
        if ~isempty(labelled_inds)
            labelled_by_annotator_1 = labelled_inds(logical(mod(labelled_inds, 2)));
            labelled_only_by_annotator_2 = setdiff(labelled_inds, [labelled_by_annotator_1, labelled_by_annotator_1+1]);
            for k=1:length(labelled_by_annotator_1)
                score_pair = full([sparse_score_mat(labelled_by_annotator_1(k)), ...
                    sparse_score_mat(labelled_by_annotator_1(k)+1)]) + 1;
                human_annot_dist_by_frame.(behav_list{i})(score_pair(1), score_pair(2)) = ...
                    human_annot_dist_by_frame.(behav_list{i})(score_pair(1), score_pair(2)) + 1;
            end
            for k=1:length(labelled_only_by_annotator_2)
                score_pair = full([0, sparse_score_mat(labelled_only_by_annotator_2(k))]) + 1;
                human_annot_dist_by_frame.(behav_list{i})(score_pair(1), score_pair(2)) = ...
                    human_annot_dist_by_frame.(behav_list{i})(score_pair(1), score_pair(2)) + 1;
            end
        end
    end
    %h = histogram(count_total);
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
end


