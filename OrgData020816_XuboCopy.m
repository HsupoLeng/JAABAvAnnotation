function OrgData020816_XuboCopy(expfolder,infofile,thres,isNoRel, maxgaps, minbouts, flymat_id, output_dir)
    if length(strfind(expfolder, '/') > 0)
        parts = strsplit(expfolder,'/');
    elseif length(strfind(expfolder,'\') > 0)
        parts = strsplit(expfolder,'\');
    end;
    experimentname = parts(end);
    src_folder = cd (expfolder); %enter experimental folder
    infos=readtable(infofile); %load information mat file
    movieflies = infos.movie;
    movies = unique(movieflies); %pick out movies represented in infofile
    frameshift=0;
    dates = dir;
    cnt=0; 
    for d=3:length(dates);
        datefolder = dates(d).name;
        if isdir(datefolder) 
        cd (datefolder);  %enter date folder
        datesessfold=dir;
        for s=3:length(datesessfold) 
            sessname = datesessfold(s).name;
            ls=length(sessname);
            for i=1:size(movies,1)
                if strcmp(sessname,strtrim(movies(i))) == 1  %if movie folder is represented in infofile
                   index=find(strcmp(movieflies,movies{i})); %find whicth entries of info file correspond to that movie
                    if any(strcmp('frameshift',fieldnames(infos)))
                        frameshift=infos.frameshift(index(1)); %only shift frames if a value is in spreadsheet
                    end;
                   if isdir(sessname)
                   cd(sessname) %enter movie folder
                   dir2 = dir;
                   for w=3:size(dir2)
                   folder3 = dir2(w).name; %no longer the movie folder- cycle through dir2
                   if isdir(folder3);
                   cd (folder3); %enter secondary movie folder
                   dir3 = dir;
                   for j=3:length(dir3);
                    if length(strfind(dir3(j).name,'JAABA')) >= 1;
                    JAABAfolder = dir3(j).name; 
                    cd (JAABAfolder); %enter JAABA folder
                    sprintf('organizing %s',JAABAfolder)
                    dir4 = dir;
                    cnt=cnt+1;
                    for n=3:length(dir4);
                        tempname = dir4(n).name;
                        if (length(strfind(tempname,'scores')) >= 1)
                            if (length(strfind(tempname,'unge')) >= 1) &&(isempty(strfind(tempname,'~')))
                            if ~isNoRel && contains(tempname, 'NoRel')
                                continue;
                            end
                            clear scoresmat; 
                            scoresmat = load(tempname);
                            A = confidence_changer_with_processingMPW(scoresmat.allScores, thres); %(Scores, threshold)
                            B = smoothing(A,frameshift,maxgaps.L,minbouts.L);  %(Scores, frameshift, max gap, min bout)
                            scoresgalore(cnt).lunge = B;
                            elseif (length(strfind(tempname,'butt')) >= 1) &&(isempty(strfind(tempname,'~')))
                            if ~isNoRel && contains(tempname, 'NoRel')
                            continue;
                            end
                            clear scoresmat; 
                            scoresmat = load(tempname);
                            A = confidence_changer_with_processingMPW(scoresmat.allScores,thres);
                            B = smoothing(A,frameshift,maxgaps.HB,minbouts.HB);
                            scoresgalore(cnt).headbutt = B;
                            elseif (length(strfind(tempname,'xt')) >= 1) &&(isempty(strfind(tempname,'~')))
                            if ~isNoRel && contains(tempname, 'NoRel')
                            continue;
                            end
                            clear scoresmat; 
                            scoresmat = load(tempname);
                            A = confidence_changer_with_processingMPW(scoresmat.allScores,thres);
                            B = smoothing(A,frameshift,maxgaps.WE,minbouts.WE);
                            scoresgalore(cnt).wingext = B;
                            elseif (length(strfind(tempname,'harge')) >= 1) &&(isempty(strfind(tempname,'~')))
                            clear scoresmat; 
                            scoresmat = load(tempname);
                            A = confidence_changer_with_processingMPW(scoresmat.allScores,thres);
                            B = smoothing(A,frameshift,maxgaps.CH,minbouts.CH);
                            scoresgalore(cnt).charge = B;
                            elseif (length(strfind(tempname,'hreat')) >= 1) &&(isempty(strfind(tempname,'~')))
                            clear scoresmat; 
                            scoresmat = load(tempname);
                            A = confidence_changer_with_processingMPW(scoresmat.allScores,thres);
                            B = smoothing(A,frameshift,maxgaps.WT,minbouts.WT);
                            scoresgalore(cnt).wingthreat = B;
                            elseif (length(strfind(tempname,'ussling')) >= 1) &&(isempty(strfind(tempname,'~')))
                            clear scoresmat; 
                            scoresmat = load(tempname);
                            A = confidence_changer_with_processingMPW(scoresmat.allScores,thres);
                            B = smoothing(A,frameshift,maxgaps.TU,minbouts.TU);
                            scoresgalore(cnt).tussling = B;
                            elseif (length(strfind(tempname,'olding')) >= 1) &&(isempty(strfind(tempname,'~')))
                            clear scoresmat; 
                            scoresmat = load(tempname);
                            A = confidence_changer_with_processingMPW(scoresmat.allScores,thres);
                            B = smoothing(A,frameshift,maxgaps.HO,minbouts.HO);
                            scoresgalore(cnt).holding = B;
                            end;
                        scoresgalore(cnt).movie = sessname;
                        scoresgalore(cnt).index = index; %adding index of flies in infofile
                        end;
                    end;
                    %you could just add information to "infos" instead of
                    %making a new thing called flymat...
                    cd ..
                    end;
                   end;
                   cd ..
                   end;
                    end;
                   cd ..
                end;
                end;
             end;
           end;
                cd ..
         end;
    end;
    filename = char(strcat('scoresgalore_',experimentname,'.mat'));
    % save(filename,'scoresgalore','-v7.3');
    % save(filename,'scoresgalore');  
    %now organize data by fly
    cnt=0;  
    flies=size(scoresgalore,2);
    for n=1:flies
        index = scoresgalore(n).index;
        sprintf('organizing %s',scoresgalore(n).movie)
    %     if length(index) > length(scoresgalore(n).lunge.scores)
    %         sprintf('more flies reported in excel then there are scores for %s',scoresgalore(n).movie)
    %         index=1:length(scoresgalore(n).lunge.scores);
    %     end;

        for s=1:length(index) % #flies per movie
                ind = index(s);
                clear flymat
                cnt=cnt+1;
                %make it so any label will be its own column.
                vnames = infos.Properties.VariableNames;
                for t=1:length(vnames);
                matname=strcat('flymat.',vnames(t));
                evalc([matname{1},' = infos.',vnames{t},'(ind)']); 
                end;
                %%I could just use infos as the base and add on to it
                if isfield(scoresgalore,'lunge')&& ~isempty(scoresgalore(n).lunge) %check if scores are there
                flymat.L_t0s=scoresgalore(n).lunge.t0s{s};
                flymat.L_t1s=scoresgalore(n).lunge.t1s{s};
                flymat.L_t4s=scoresgalore(n).lunge.t4s{s};
                flymat.L_t5s=scoresgalore(n).lunge.t5s{s};
                flymat.L_binary=scoresgalore(n).lunge.binary{s};
                flymat.L_total=length(scoresgalore(n).lunge.t5s{s});
                flymat.L_dur=sum(flymat.L_binary);
                if any(strcmp('fps',fieldnames(infos)))
                    [mb,mbd] = makeminbins(flymat.L_t4s,flymat.L_binary,infos.fps(1));
                    flymat.L_minbin = mb;
                    flymat.L_minbindur = mbd;
                end;
                end;
                if isfield(scoresgalore,'headbutt') && ~isempty(scoresgalore(n).headbutt) %check if scores are there
                flymat.HB_t0s=scoresgalore(n).headbutt.t0s{s};
                flymat.HB_t1s=scoresgalore(n).headbutt.t1s{s};            
                flymat.HB_t4s=scoresgalore(n).headbutt.t4s{s};
                flymat.HB_t5s=scoresgalore(n).headbutt.t5s{s};
                flymat.HB_binary=scoresgalore(n).headbutt.binary{s};
                flymat.HB_total=length(scoresgalore(n).headbutt.t5s{s});
                flymat.HB_dur=sum(flymat.HB_binary);
                if any(strcmp('fps',fieldnames(infos)))
                    [mb,mbd] = makeminbins(flymat.HB_t4s,flymat.HB_binary,infos.fps(1));
                    flymat.HB_minbin = mb;
                    flymat.HB_minbindur = mbd;
                end;
                end;
                if isfield(scoresgalore,'wingext')&& ~isempty(scoresgalore(n).wingext) %check if scores are there
                flymat.WE_t0s=scoresgalore(n).wingext.t0s{s};
                flymat.WE_t1s=scoresgalore(n).wingext.t1s{s};
                flymat.WE_t4s=scoresgalore(n).wingext.t4s{s};
                flymat.WE_t5s=scoresgalore(n).wingext.t5s{s};
                flymat.WE_binary=scoresgalore(n).wingext.binary{s};
                flymat.WE_total=length(scoresgalore(n).wingext.t5s{s});
                flymat.WE_dur=sum(flymat.WE_binary);
                if any(strcmp('fps',fieldnames(infos)))
                [mb,mbd] = makeminbins(flymat.WE_t4s,flymat.WE_binary,infos.fps(1));
                flymat.WE_minbin = mb;
                flymat.WE_minbindur = mbd;
                end;
                end;
                if isfield(scoresgalore,'charge')&& ~isempty(scoresgalore(n).charge) %check if scores are there
                flymat.C_t0s=scoresgalore(n).charge.t0s{s};
                flymat.C_t1s=scoresgalore(n).charge.t1s{s};
                flymat.C_t4s=scoresgalore(n).charge.t4s{s};
                flymat.C_t5s=scoresgalore(n).charge.t5s{s};
                flymat.C_binary=scoresgalore(n).charge.binary{s};
                flymat.C_total=length(scoresgalore(n).charge.t5s{s});
                flymat.C_dur=sum(flymat.C_binary);
                if any(strcmp('fps',fieldnames(infos)))
                    [mb,mbd] = makeminbins(flymat.C_t4s,flymat.C_binary,infos.fps(1));
                    flymat.C_minbin = mb;
                    flymat.C_minbindur = mbd;
                end;
                end;
                if isfield(scoresgalore,'wingthreat')&& ~isempty(scoresgalore(n).wingthreat) %check if scores are there
                flymat.WT_t0s=scoresgalore(n).wingthreat.t0s{s};
                flymat.WT_t1s=scoresgalore(n).wingthreat.t1s{s};
                flymat.WT_t4s=scoresgalore(n).wingthreat.t4s{s};
                flymat.WT_t5s=scoresgalore(n).wingthreat.t5s{s};
                flymat.WT_binary=scoresgalore(n).wingthreat.binary{s};
                flymat.WT_total=length(scoresgalore(n).wingthreat.t5s{s});
                flymat.WT_dur=sum(flymat.WT_binary);
                if any(strcmp('fps',fieldnames(infos)))
                    [mb,mbd] = makeminbins(flymat.WT_t4s,flymat.WT_binary,infos.fps(1));
                    flymat.WT_minbin = mb;
                    flymat.WT_minbindur = mbd;
                end;
                end;
                if isfield(scoresgalore,'holding')&& ~isempty(scoresgalore(n).holding) %check if scores are there
                flymat.HO_t0s=scoresgalore(n).holding.t0s{s};
                flymat.HO_t1s=scoresgalore(n).holding.t1s{s};
                flymat.HO_t4s=scoresgalore(n).holding.t4s{s};
                flymat.HO_t5s=scoresgalore(n).holding.t5s{s};
                flymat.HO_binary=scoresgalore(n).holding.binary{s};
                flymat.HO_total=length(scoresgalore(n).holding.t5s{s});
                flymat.HO_dur=sum(flymat.HO_binary);
                if any(strcmp('fps',fieldnames(infos)))
                    [mb,mbd] = makeminbins(flymat.HO_t4s,flymat.HO_binary,infos.fps(1));
                    flymat.HO_minbin = mb;
                    flymat.HO_minbindur = mbd;
                end;
                end;
                if isfield(scoresgalore,'tussling')&& ~isempty(scoresgalore(n).tussling) %check if scores are there
                flymat.T_t0s=scoresgalore(n).tussling.t0s{s};
                flymat.T_t1s=scoresgalore(n).tussling.t1s{s};
                flymat.T_t4s=scoresgalore(n).tussling.t4s{s};
                flymat.T_t5s=scoresgalore(n).tussling.t5s{s};
                flymat.T_binary=scoresgalore(n).tussling.binary{s};
                flymat.T_total=length(scoresgalore(n).tussling.t5s{s});
                flymat.T_dur=sum(flymat.T_binary);
                if any(strcmp('fps',fieldnames(infos)))
                    [mb,mbd] = makeminbins(flymat.T_t4s,flymat.T_binary,infos.fps(1));
                    flymat.T_minbin = mb;
                    flymat.T_minbindur = mbd;
                end;
                end;
                flymatAll(cnt)=flymat;  %add all info for that fly to the big experimental structure
    %         end;
        end;
        end;
        
    for i=1:length(flymat_id)
        filename = fullfile(output_dir, strcat('FLYMAT_MNL-KA JAABA training samples_', flymat_id{i}, '.mat'));
        save(filename,'flymatAll','-v7.3');
    end
    cd(src_folder);
end