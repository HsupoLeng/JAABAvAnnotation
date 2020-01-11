bout_fn_fp_list_ori = dir('D:\xubo\code\annot-analysis\true_false_count\bout_fn_fp_all');
bout_fn_fp_list_ori = {bout_fn_fp_list_ori(:).name};
bout_fn_fp_list_new = dir('D:\xubo\code\annot-analysis\true_false_count_new\bout_fn_fp_all');
bout_fn_fp_list_new = {bout_fn_fp_list_new(:).name};

bout_mismatch_list = {};
for i=3:length(bout_fn_fp_list_ori)
    s_ori = load(fullfile('D:\xubo\code\annot-analysis\true_false_count\bout_fn_fp_all', bout_fn_fp_list_ori{i}));
    s_new = load(fullfile('D:\xubo\code\annot-analysis\true_false_count_new\bout_fn_fp_all', bout_fn_fp_list_new{i}));
    if ~isequal(s_ori, s_new)
        bout_mismatch_list = {bout_mismatch_list{:}, bout_fn_fp_list_ori{i}};
    end
end