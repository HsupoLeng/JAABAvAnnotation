function fly_idx = find_fly_in_flymat(flymat, movie_name, fly_number, is_cell_arr)
    if is_cell_arr
        movie_arr = [flymat.movie];
    else
        movie_arr = {flymat.movie};
    end
    
    if fly_number
        match_fly_num = [flymat.fly] == fly_number;
    else
        match_fly_num = ones(size(movie_arr));
    end
    
    fly_idx = find(bitand(strcmp(movie_arr, movie_name), ...
        match_fly_num));
end