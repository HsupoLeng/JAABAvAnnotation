function intsct_violt = sanity_check_humanAnnot(behav_list, bout_matches_all)
    if ~exist('behav_list', 'var') || isempty(behav_list)
        load('common-params-annot-analysis.mat');
    end
    if ~exist('bout_matches_all', 'var') || isempty(bout_matches_all)
        load('bout_matches_all.mat');
    end

    intsct_violt = [];
    for i=1:length(behav_list)
        for j=1:length(bout_matches_all.(behav_list{i}))
            entry = bout_matches_all.(behav_list{i})(j);
            if ~isnan(entry.annot_union_start)
                if entry.annot_intsct_start > entry.annot_intsct_end
                    intsct_violt = [intsct_violt, [i,j]];
                end
            end
        end
    end
    intsct_violt = reshape(intsct_violt, 2, [])';
end
