%start with manvjaaba Structure.  length(manvjaaba)is the number of movies
%being annotated.  Within manvjaaba.movie you will find the name of the
%movie, the IDs of the flies used.  You will find the start frame, end
%frame, and confidence value of a bout in observer1 and observer2.
%Predicted is the start and end frame generated by the JAABA classifier.  

%this code will combine observer confidence values for each frame and 
%then compare it to predictions.  This will be concatenated into a
%2x(frames*flies) matrix for all movies were (1,:) is confidence sum 
%and (2,:) is the binary behavior decision of the classifier

% MODIFIED from MANvsJAABAduration, with following changes
% 1. compute conventional precision and recall
% 2. automatically draw colormap legend
% 3. automatically label the chosen set of parameters
function MVJ=MANvsJAABAduration_xubo(filename,thresholds,maxgaps,minbehavs, ourchoice_idx)
    load(filename);
    cnt=0;
    for th=1:length(thresholds)
        for gaps=1:length(maxgaps)
            for beh=1:length(minbehavs)
                    allobsframes=[]; allpreframes=[];
                    threshold=thresholds(th);
                    maxgap=maxgaps(gaps); 
                    minbehav=minbehavs(beh);
                    cnt=cnt+1;
                for n=1:length(manvjaaba)
                %access movie 
                frs=manvjaaba(n).movie.frames; %need to know the number of frames to look for
                flies=manvjaaba(n).movie.flies;
                for fl=1:length(flies)
                    allS=manvjaaba(n).movie.allScores;
                    A = changeJAABAconfidenceValCutoff(allS,threshold); %threshold
                    B = smoothBehavBouts(A,0,maxgap,minbehav); %close gaps, eliminate short behaviors
                    obs1=manvjaaba(n).movie.observer1{fl};
                    obs2=manvjaaba(n).movie.observer2{fl};
                    tmpobs=zeros(2,frs); %prealocate frames array
                    tmppre=zeros(1,frs);
                    for o=1:size(obs1,2) %transform observer 1 start and end frames into 1 X length(frames) array of confidence values
                        row=obs1(:,o); %pull out each bout
    %                     if row(2)-row(1)-1 >= minbehav %get rid of short behaviors
                        tmpobs(1,row(1):(row(2)-1))=row(3); %populate the zeros array with confidence values
    %                     end;
                    end;
                    for t=1:size(obs2,2)
                        row=obs2(:,t);
    %                     if row(2)-row(1)-1 >= minbehav 
                        tmpobs(2,row(1):(row(2)-1))=row(3);
    %                     end;
                        tmpobs=tmpobs(:,1:frs);
                    end
                    pre=[];
                    if ~isempty([B.startsm{flies(fl)}])
                    pre(1,:)=[B.startsm{flies(fl)}];
                    pre(2,:)=[B.endsm{flies(fl)}];
                    end;
                    for p=1:size(pre,2)
                        row=pre(:,p);
                        tmppre(row(1):(row(2)-1))=1;
                        tmppre=tmppre(:,1:frs);
                    end;   
                    tmpsum=sum(tmpobs); %add confidence values of two observers for each frame
                    allobsframes=[allobsframes,tmpsum]; %make one big array of frame values for all movies
                    allpreframes=[allpreframes,tmppre];
                    
                end;
    end;

    for n=1:6 %6 is amount of possible confidence values
        tmpfr=allpreframes(find(allobsframes==n));
        yesclass(n)=length(find(tmpfr==1)); %sum of frames classified as a behavior for each confidence value sum
        noclass(n)=length(find(tmpfr==0)); %frames not behavior
    end;
    %find false positives
    FPobs=allobsframes(find(allpreframes==1));
    FP_original=length(find(FPobs<4)); %find all false positives
    FP_new = length(find(FPobs<1));
    HCPdetected=sum(yesclass(4:6));
    HCPmissed=sum(noclass(4:6));
    LCPdetected=sum(yesclass(1:3));
    LCPmissed=sum(noclass(1:3));
    precision_original=HCPdetected/(HCPdetected+FP_original);
    recall_original=HCPdetected/(HCPdetected+HCPmissed);
    %The following two lines compute the conventional precision and recall, i.e. 
    % precision = TP/(TP + FP); recall = (TP/(TP + FN)). 
    % It does not distinguish between high-confidence bouts and
    % low-confidence bouts. 
    precision_new = sum(LCPdetected+HCPdetected)/(sum(LCPdetected+HCPdetected)+FP_new);
    recall_new = sum(LCPdetected+HCPdetected)/(sum(LCPdetected+HCPdetected)+sum(HCPmissed+LCPmissed));
    MVJ(cnt,:)=[threshold maxgap minbehav precision_original recall_original precision_new recall_new];
                end;
            end;
          MVJ
    end;

    %GRAPH VALUES
    for i=1:2
        if i == 1
            continue; % We will be using the conventional definition of precision and accuracy
        end
        figure();
        % Generate color values corresponding to the parameter values.
        % Larger the parameter, brighter the color. 
        % Each parameter takes one color channel (r, g, or b)
        conf_ax = axes;
        max_ax = axes;
        min_ax = axes;
        confcolor_mat = [zeros(size(MVJ, 1), 1), (MVJ(:, 1)/(0.2*2)) + 0.25, zeros(size(MVJ, 1), 1)];
        maxcolor_mat = [(MVJ(:, 2)/max(MVJ(:,2))), zeros(size(MVJ, 1), 1), zeros(size(MVJ, 1), 1)];
        mincolor_mat = [zeros(size(MVJ, 1), 1), zeros(size(MVJ, 1), 1), (MVJ(:,3)/max(MVJ(:, 3)))];
        
        hold on;
        scatter(conf_ax, MVJ(:, 4+2*(i-1)), MVJ(:, 5+2*(i-1)), 500, confcolor_mat, 's', 'LineWidth', 1.5);
        scatter(max_ax, MVJ(:, 4+2*(i-1)), MVJ(:, 5+2*(i-1)), 200, maxcolor_mat, 'd', 'LineWidth', 1.5);
        scatter(min_ax, MVJ(:, 4+2*(i-1)), MVJ(:, 5+2*(i-1)), 75, mincolor_mat, 'o', 'LineWidth', 1.5);
        linkaxes([conf_ax, max_ax, min_ax]);
        max_ax.Visible = 'off';
        max_ax.XTick = [];
        max_ax.YTick = [];
        min_ax.Visible = 'off';
        min_ax.XTick = [];
        min_ax.YTick = [];
        % One colormap for each parameter
        colormap(conf_ax, unique(confcolor_mat, 'rows'));
        colormap(max_ax, unique(maxcolor_mat, 'rows'));
        colormap(min_ax, unique(mincolor_mat, 'rows'));
        set([conf_ax, max_ax, min_ax], 'Position', [0.10, 0.11, 0.685, 0.815]);
        colorbar(conf_ax, 'Position', [0.79, 0.11, 0.0675, 0.815], 'Ticks', []);
        colorbar(max_ax, 'Position', [0.86, 0.11, 0.0675, 0.815], 'Ticks', []);
        colorbar(min_ax, 'Position', [0.93, 0.11, 0.0675, 0.815], 'Ticks', []);
        %plot the value you chose in the end as a black tick mark
        scatter(min_ax, MVJ(ourchoice_idx, 4+2*(i-1)), MVJ(ourchoice_idx, 5+2*(i-1)), 75,'+','MarkerEdgeColor',[0 0 0],'LineWidth',1.5)
        hold off
        xlabel(conf_ax, 'Precision');
        ylabel(conf_ax, 'Recall');

        [~, file_id, ~] = fileparts(filename);
        file_id_elems = strsplit(file_id, '-');
        file_id = strjoin(file_id_elems(2:end), '-');
        if i==1
            saveas(gcf, strcat('ManualvsJAABA_grid_performance_original-', file_id, '.png'), 'png');
            saveas(double(gcf), strcat('ManualvsJAABA_grid_performance_original-', file_id, '.eps'), 'eps');
        else
            saveas(gcf, strcat('ManualvsJAABA_grid_performance_new-', file_id, '.png'), 'png');
            saveas(double(gcf), strcat('ManualvsJAABA_grid_performance_new-', file_id, '.eps'), 'eps');
        end
    end
end
