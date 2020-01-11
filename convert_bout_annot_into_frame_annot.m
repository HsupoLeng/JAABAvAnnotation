function [count, frame_score_mat_sparse] = convert_bout_annot_into_frame_annot(t0s, t1s, scores)
    % This program constructs per frame annotation from FLYMAT-like input. 
    % Since FLYMAT can have inaccurate entries, use of this program is not
    % recommended. 
    % Use convert_bout_match_into_frame_match instead.
    count = 0;
    max_timeloc = max(t1s(:));
    if isempty(t0s)
        frame_score_mat_sparse = sparse([]);
        return;
    end
    num_of_nonzero_frames = sum(max(t1s, [], 2) - min(t0s, [], 2) + 1);
    frame_score_mat_sparse = spalloc(2, max_timeloc, num_of_nonzero_frames);
    scores(isnan(scores)) = 0;
    
    for i=1:size(t0s, 1)
        if t0s(i, 1) < t0s(i, 2)
            pre_segment = t0s(i,1):t0s(i,2)-1;
            frame_score_mat_sparse(:, pre_segment) = ...
                repmat([scores(i, 1); 0], 1, length(pre_segment));
            
            if(scores(i,1)) <= 3 
                if length(pre_segment) > 500
                    fprintf('Something is wrong! Length of pre-segment:%d\n', length(pre_segment));
                end
                count = [count, length(pre_segment)];
            end
        elseif t0s(i, 1) > t0s(i, 2)
            pre_segment = t0s(i,2):t0s(i,1)-1;
            frame_score_mat_sparse(:, pre_segment) = ...
                repmat([0; scores(i, 2)], 1, length(pre_segment));
            
            if(scores(i,2)) <= 3
                 if length(pre_segment) > 500
                    fprintf('Something is wrong! Length of pre-segment:%d\n', length(pre_segment));
                end
                count = [count, length(pre_segment)];
            end
        end
        overlap_segment = max(t0s(i,:)):(min(t1s(i,:)) - 1);
        for j=1:length(overlap_segment)
            frame_score_mat_sparse(:, overlap_segment(j)) = scores(i, :);
        end
        if ~isempty(overlap_segment)
            if(scores(i,2) == 3 && scores(i,1) == 0) || (scores(i,1) == 3 && scores(i,2) == 0)
                if length(overlap_segment) > 800
                    fprintf('Something is wrong! Length of overlap-segment:%d\n', length(overlap_segment));
                end
                count = [count, length(overlap_segment)];        
            end
        end
        if t1s(i, 1) > t1s(i, 2)
            post_segment = t1s(i,2):(t1s(i,1)-1);
            frame_score_mat_sparse(:, post_segment) = ...
                repmat([scores(i, 1); 0], 1, length(post_segment));
            
            if(scores(i,1)) <= 3
                if length(post_segment) > 500
                    fprintf('Something is wrong! Length of post-segment:%d\n', length(post_segment));
                end
                count = [count, length(post_segment)];
            end
        elseif t1s(i, 1) < t1s(i, 2)
            post_segment = t1s(i,1):(t1s(i,2)-1);
            frame_score_mat_sparse(:, post_segment) = ...
                repmat([0; scores(i, 2)], 1, length(post_segment));
            
            if(scores(i,2)) <= 3 
                 if length(post_segment) > 500
                    fprintf('Something is wrong! Length of post-segment:%d\n', length(post_segment));
                end
                count = [count, length(post_segment)];
            end
        end
    end
end