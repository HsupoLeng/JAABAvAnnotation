function B = smoothing(allScores,frameshift,max_gap,min_bout);
% takes scored behaviors and gets rid of behaviors below min_bout

% INPUTS: max_gap = maximum isi for which two bouts will be smooshed together
% min_bout = minimum number of frames to call something a behavior

% OUTPUT: B - allScores structure with added fieldnames and values that 
% reflect recalculated behavior data based on frameshift, max_gap, and min_bout 


% For behaviors that didn't need to go through thresholding (don't have
% thresholded behavior data)
if ~isfield(allScores,'startth')   
    allScores.startth = allScores.t2s;
    allScores.endth = allScores.t3s;
end;

% Obtain the number of flies in a given movie
flies_n = size(allScores.t0s, 2);   
frames = length(allScores.postprocessed{1});
for p = 1:flies_n;
    bouts_n = size(allScores.startth{p}, 2); 
    
    % Apply frame shift specified in info file
    startth = allScores.startth{p} + frameshift;
    endtth = allScores.endth{p} + frameshift;
    
    %when 
    if ~(length(startth) == length(endtth))
        disp('start and end frames do not have equal elements')
        length(startth)
        length(endtth)
        if length(startth) > length(endtth)
            startth = startth(1:length(endtth));
        else
            endtth = endtth(1:length(startth));
        end;
    end;
    
    %Merging bouts with short gap (2 frame seems reasonable for WingExts).
    isi = startth(2:end)- endtth(1:(end-1)); %Gap between neighboring bouts
    startth((find(isi <= max_gap) + 1))=[]; 
    endtth(find(isi <= max_gap))=[]; 
    %only keeping behaviors that last longer than or equal to min_bout frames 
    boutlength = endtth - startth;
    %get rid of bouts that are too short
    start = startth(find(boutlength >= min_bout));
    %get rid of t4s that may be have been shifted past the end of
    %the movie
    start = start(find(start < frames)); 
    finish = endtth(find(boutlength >= min_bout));
     %get rid of t5s that may be have been shifted past the end of
    %the movie
    finish = finish(find(finish <= frames));
    if length(start) > length(finish)
        finish=[finish,frames];
    end;
    allScores.t4s{p} = start;
    allScores.t5s{p} = finish;
    bouts = length(start); 
    binary = zeros(1,frames);
    for b = 1:bouts
    if start(b) > 0;
    binary(start(b):finish(b))=1;
    end;
    end;
    allScores.binary{p} = binary;
end;

B=allScores;
    
    
% %         







