load('common-params-annot-analysis.mat');
all_paths_in_dir = dir('D:\xubo\code\annot-analysis\bout_frame_matches');
all_files = all_paths_in_dir(~[all_paths_in_dir.isdir]);
bout_frame_matches_file_mask = regexp({all_files(:).name}, '.+_matches_.+\.mat');
bout_frame_matches_file_mask = ~cellfun(@isempty, bout_frame_matches_file_mask);
all_bout_frame_matches = all_files(bout_frame_matches_file_mask);

behav_sel_strs = ['all', behav_shorthands];
behav_sel_options = {[1,2,3], 1,2,3};

for i=1:length(all_bout_frame_matches)
    match_file = all_bout_frame_matches(i).name; 
    [~, match_str, ~] = fileparts(match_file);
    match_str_elems = strsplit(match_str, '_');
    behav_sel_str = regexp(match_str_elems{3}, '[a-zA-Z]*', 'match', 'once');
    behav_sel = behav_sel_options{strcmpi(behav_sel_strs, behav_sel_str)};
    if strcmp(match_str_elems{1}, 'frame')
        is_frame_matches = true;
    else
        is_frame_matches = false; 
    end
    compile_fp_fn_numbers(behav_sel, strjoin(match_str_elems(3:end), '_'), is_frame_matches);
end