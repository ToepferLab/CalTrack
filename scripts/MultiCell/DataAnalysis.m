
% DataAnalysis_MultiCell.m

path = cd;
file = 'CalciumData_MultiCell.xlsx';
s_names = sheetnames(fullfile(path,file));

for i=1:length(s_names)
    ca_matrix_original(i).data = readmatrix(fullfile(path,file),'Sheet',i);
end

global AcquisitionFrequency
global cutoff
global PacingFrequency

for i=1:length(ca_matrix_original)
    figure('Name', strcat('video',num2str(i)))
    for j=1:size(ca_matrix_original(i).data,2)
        subplot(5,ceil(size(ca_matrix_original(i).data,2)/5),j)
        t = (1:length(ca_matrix_original(i).data(:,j)))/AcquisitionFrequency;
        plot(t,ca_matrix_original(i).data(:,j)), hold on, grid on

        if M==0
            plot(t(N+1:end),ca_matrix_original(i).data(N+1:end,j))
            ca_matrix(i).data(:,j) = ca_matrix_original(i).data(N+1:end,j);
        else
            plot(t(N+1:M),ca_matrix_original(i).data(N+1:M,j))
            ca_matrix(i).data(:,j) = ca_matrix_original(i).data(N+1:M,j);
        end
    
        title(strcat('cell',num2str(j)))
        yticklabels(' ')
    end
    ca_matrix(i).name = s_names(i);

end

remove =questdlg('Would you like to remove any cell?','Filter undesired cells',...
    'Yes','No ', 'Yes');

if remove == 'Yes'
    prompt = {'\bf \fontsize{12} Please enter the video number for those videos that contain the cells to be removed (from figures title): e.g. [1 2 3] (if more than one cell from the same video, than repeat the video number)',...
    '\bf \fontsize{12} Please enter the cell number for the cells to be removed (subplots title): e.g. [1 2 3] (their order needs to match the videos order and when you want to remove more than 1 cell from the same video, report the cell number from the largest to the smallest'};
    dlgtitle = 'Remove cells';
    dims = [1 120];
    definput = {'','', '', '','','','',''};
    opts.Interpreter = 'tex';
    RC = inputdlg(prompt,dlgtitle,dims,definput, opts);
    rem_vid = str2num(RC{1});
    rem_cell = str2num(RC{2});
end

if remove == 'y'
    for i=1:length(rem_vid)
        ca_matrix(rem_vid(i)).data(:,rem_cell(i)) = [];
    end
end  

j=1;
for i=1:length(ca_matrix)
    if (isempty(ca_matrix(i).data)==1)
    id(j)=i;
    j=j+1;
    end
end
if exist('id') == 1;
    ca_matrix(id) = [];
end

photobl =questdlg('Do you need to apply photo bleach correction?','Photobleach correction toggle',...
    'Yes','No ', 'Yes');

video = [];

if segmentation == 'n' 
    [video, Er, Sk, extra_beat] = BeatSegmentation(ca_matrix);
elseif segmentation == 'y'
    T0 = xlsread('T0.xlsx'); 
    [video, Er, Sk, extra_beat] = BeatSegmentation2(ca_matrix,T0);
end

if isempty(Sk)
    clearvars Sk
end
if isempty(extra_beat)
    clearvars extra_beat
end
if isempty(Er)
    clearvars Er
end


if photobl == 'Yes'
   
    [video] = PhotoBleachCorrection(video);
    disp('Photo Bleach Correction done')
    for i = 1 : length(video)
        mat = [];
        for j = 1: length(video(i).single_traces)
            len(j) = length(video(i).single_traces(j).BaselineCorrectedTraces);
        end
        for j = 1: length(video(i).single_traces)
            mat = [mat , video(i).single_traces(j).BaselineCorrectedTraces(1:min(len))];
        end
        ca_matrix_phbl(i).data = mat;
        ca_matrix_phbl(i).name = video(i).name;
    end
    
    if segmentation == 'n'
        [video_phbl, errors_phbl, skipped_phbl, extra_beat_phbl] = BeatSegmentation(ca_matrix_phbl);
    elseif segmentation == 'y'
        T02 = T0;
        T02(:,2) = ones(size(T02(:,2)));
        [video_phbl, errors_phbl, skipped_phbl, extra_beat_phbl] = BeatSegmentation2(ca_matrix_phbl,T02);
    end

    if isempty(skipped_phbl)
        clearvars skipped_phbl
    elseif exist('skipped')
        skipped = [skipped skipped_phbl];
    else
        skipped = skipped_phbl;
    end
    
    if isempty(extra_beat_phbl)
        clearvars extra_beat_phbl
    elseif exist('extra_beat')
        extra_beat = [extra_beat extra_beat_phbl];
    else
        extra_beat = extra_beat_phbl;
    end
    
    if isempty(errors_phbl)
        clearvars errors_phbl
    elseif exist('errors')
        errors = [errors errors_phbl];
    else
        errors = errors_phbl;
    end
  
    
    for i=1:length(video)
        figure('Name',['Baseline Correction: video',num2str(i)]), clf
        for j = 1: length(video(i).single_traces)
            subplot(5,ceil(length(video(i).single_traces)/5),j)

            selected_trace = [];

            cll = {video(i).single_traces(j).data.beats}.';
            selected_trace = cat(1,cll{:});
            
            t = 1:length(selected_trace);     
            plot(t,selected_trace), hold on, grid on

            tt = 1:length(video(i).single_traces(j).BaselineCorrectedTraces);
            plot(tt,video(i).single_traces(j).BaselineCorrectedTraces)

            title(strcat('cell',num2str(j)))
            yticklabels(' ')
            xticklabels(' ')
        end
    end
    
    video_before_phbl = video;
    video = [];
    video = video_phbl;

    remove =questdlg('Would you like to remove any cell?','Filter undesired cells',...
        'Yes','No ', 'Yes');

    if remove == 'Yes'
        prompt = {'\bf \fontsize{12} Please enter the video number for those videos that contain the cells to be removed (from figures title): e.g. [1 2 3] (if more than one cell from the same video, than repeat the video number)',...
        '\bf \fontsize{12} Please enter the cell number for the cells to be removed (subplots title): e.g. [1 2 3] (their order needs to match the videos order and when you want to remove more than 1 cell from the same video, report the cell number from the largest to the smallest'};
        dlgtitle = 'Remove cells';
        dims = [1 120];
        definput = {'','', '', '','','','',''};
        opts.Interpreter = 'tex';
        RC = inputdlg(prompt,dlgtitle,dims,definput, opts);
        rem_vid = str2num(RC{1});
        rem_cell = str2num(RC{2});
    end

    if remove == 'Yes'
        for i=1:length(rem_vid)
            video(rem_vid(i)).single_traces(rem_cell(i) )= [];
        end
    end

end

if_skipped = exist('Sk');
if_error = exist('Er');
if_extra_beat = exist('extra_beat');

%%
V = [];
vid = [];
non_regular_beat = [];

% if param_single_beat == 'y'
%     
%    k = 1; 
%     for v = 1:length(video)
%         V(v).name  = video(v).name;
%         ii = 1;
%         for i=1:length(video(v).single_traces)
%             j = 1;
%             for jj=1:length(video(v).single_traces(i).data)
%                 yy = video(v).single_traces(i).data(jj).beats;
% 
%                 % fit
%                 [PP,XX] = findpeaks(yy,'SortStr','descend');
%                 half_vv = yy(XX(1):end);
%                 xx = 1:length(half_vv);
%                 CL = 1/PacingFrequency; %s
%                 n_time = 0.2*CL;
%                 n = ceil(n_time*AcquisitionFrequency);
%                 bsl = mean(half_vv(end-n:end));
%                 bsln = ((ones (size(half_vv,1),1) .*bsl))';
%                 [xint, yint] = intersections(xx, half_vv, xx, bsln, ' robust');
%                 maxx = floor(min(xint))+ 3;
%                 if maxx > length (half_vv)
%                     maxx = length (half_vv);
%                 end
%                 half_vv2 = half_vv(1:maxx);
%                 xx2 = 1:length(half_vv2);
%                 half_vv3 = half_vv(maxx+1:end);
%                 xx3 = 1:length(half_vv3);
%                 [p_EAD] = findpeaks(half_vv2,'MinPeakProminence',0.2*(max(half_vv2-min(half_vv2))));
% 
% 
%                 [fitresult2, gof2, xData2, yData2] = fit_tau(xx2, half_vv2');
%                 vid(v).Sbeat(i).beat(jj).tau = ((1/(-fitresult2.b)).*(1/AcquisitionFrequency))*1000; %ms
%                 vid(v).Sbeat(i).beat(jj).a = fitresult2.a;
%                 vid(v).Sbeat(i).beat(jj).c = fitresult2.c;
%                 vid(v).Sbeat(i).beat(jj).rsquare = gof2.rsquare; 
% 
%                 if length(half_vv3) > 3
%                     pks = [];
%                     warning off
%                     if AcquisitionFrequency > 100
%                         [pks,locs,w,p] = findpeaks(movmean(half_vv3,10),'MinPeakHeight',bsl+0.1*bsl);  
%                     else
%                         [pks,locs,w,p] = findpeaks(half_vv3,'MinPeakProminence',0.03*mean(half_vv3));
%                     end
%                 end
% 
%                 if ((isempty(p_EAD) == 0) || (isempty(pks) == 0) 
%                     non_regular_beat(k).video_n = v;
%                     non_regular_beat(k).cell_n = i;
%                     non_regular_beat(k).beat_n = jj;
%                     non_regular_beat(k).cell_trace = video(v).single_traces(i).original;
%                     non_regular_beat(k).beat_trace = yy;
%                     k = k+1;
% 
%                 else
%                     if AcquisitionFrequency>100
%                         Y = movmean(yy,4);
%                     else
%                         Y = yy;
%                     end
% 
%                     [CalciumPeak, CTD, CTD90, CTD50, CTD10, Time_C_peak, T_C_relax, AvgCalciumTrace, DoverD0,...
%                     yAVGCalciumTrace_baseline, Time_to_90a, Time_to_50a, Time_to_10a, Time_to_10Relax,...
%                     Time_to_50Relax, Time_to_90Relax, CalciumMagnitude,CalciumTDPre,CalciumTDPost] = parameters(Y',AcquisitionFrequency);
%                
%                     p = [CalciumPeak, CTD, CTD90, CTD50, CTD10, Time_C_peak, T_C_relax, DoverD0,...
%                     yAVGCalciumTrace_baseline, Time_to_90a, Time_to_50a, Time_to_10a, Time_to_10Relax,...
%                     Time_to_50Relax, Time_to_90Relax, CalciumMagnitude,CalciumTDPre,CalciumTDPost];
% 
%                     if length(p) == 18
%                         V(v).single_beat(ii).beat(j).DoverD0 = DoverD0;
%                         V(v).single_beat(ii).beat(j).Baseline = yAVGCalciumTrace_baseline;
%                         V(v).single_beat(ii).beat(j).CalciumMagnitude = CalciumMagnitude;
%                         V(v).single_beat(ii).beat(j).Mean_Calcium_CD =  CTD * 1000;
%                         V(v).single_beat(ii).beat(j).Mean_Calcium_CD90 = CTD90 * 1000;
%                         V(v).single_beat(ii).beat(j).Mean_Calcium_CD50 = CTD50 * 1000; 
%                         V(v).single_beat(ii).beat(j).Mean_Calcium_CD10Max = CTD10 * 1000;
%                         V(v).single_beat(ii).beat(j).Time_CalciumPeak = Time_C_peak * 1000;
%                         V(v).single_beat(ii).beat(j).Time_CalciumRelax = T_C_relax * 1000;
%                         V(v).single_beat(ii).beat(j).Time_to_10Contr = Time_to_90a * 1000;
%                         V(v).single_beat(ii).beat(j).Time_to_50Contr = Time_to_50a * 1000;
%                         V(v).single_beat(ii).beat(j).Time_to_90Contr = Time_to_10a * 1000;
%                         V(v).single_beat(ii).beat(j).Time_to_10_Relax = Time_to_10Relax * 1000;
%                         V(v).single_beat(ii).beat(j).Time_to_50_Relax = Time_to_50Relax * 1000;
%                         V(v).single_beat(ii).beat(j).Time_to_90_Relax = Time_to_90Relax * 1000;
%                         V(v).single_beat(ii).beat(j).t0 = CalciumTDPre * 1000;
%                         V(v).single_beat(ii).beat(j).tend = CalciumTDPost * 1000;
% 
%                         V(v).single_beat(ii).beat(j).tau = ((1/(-fitresult2.b)).*(1/AcquisitionFrequency))*1000; %ms
%                         V(v).single_beat(ii).beat(j).a = fitresult2.a;
%                         V(v).single_beat(ii).beat(j).c = fitresult2.c;
%                         V(v).single_beat(ii).beat(j).rsquare = gof2.rsquare;
% 
%                         V(v).single_beat(ii).data(j).trace = yy;
%                          
%                         j=j+1;
%                     end
%                 end
%             end
%             V(v).single_beat(ii).BR = video(v).single_traces(i).BR;
%             V(v).single_beat(ii).BBdistance = video(v).single_traces(i).BBdistance;
%             V(v).single_beat(ii).Nperiod = video(v).single_traces(i).Nperiod;
% 
%             ii = ii+1;
%         
%         end
%     end
%     
%     for v = 1:length(V)
%         j = 1;
%         index = [];
%         for i = 1:length(V(v).single_beat)
%             if (isempty(V(v).single_beat(i).beat))
%                 index(j) = i;
%                 j = j+1;
%             end       
%         end
%         V(v).single_beat(index) = [];
%     end
%     
%     STD_v = [];
%     for v = 1:length(V)
%         for i = 1:length(V(v).single_beat)
%             matr = [];
%             matr = [V(v).single_beat(i).beat.Baseline; V(v).single_beat(i).beat.DoverD0;...
%              V(v).single_beat(i).beat.Mean_Calcium_CD;...
%             V(v).single_beat(i).beat.Mean_Calcium_CD90; V(v).single_beat(i).beat.Mean_Calcium_CD50;...
%             V(v).single_beat(i).beat.Mean_Calcium_CD10Max; V(v).single_beat(i).beat.Time_CalciumPeak;...
%             V(v).single_beat(i).beat.Time_CalciumRelax; V(v).single_beat(i).beat.Time_to_10Contr;...
%             V(v).single_beat(i).beat.Time_to_50Contr; V(v).single_beat(i).beat.Time_to_90Contr;...
%             V(v).single_beat(i).beat.Time_to_10_Relax; V(v).single_beat(i).beat.Time_to_50_Relax;...
%             V(v).single_beat(i).beat.Time_to_90_Relax;...
%             V(v).single_beat(i).beat.tau; V(v).single_beat(i).beat.a; V(v).single_beat(i).beat.c; V(v).single_beat(i).beat.rsquare;...
%             V(v).single_beat(i).beat.t0; V(v).single_beat(i).beat.tend].';
% 
%             [N_beat(i),C] = size(matr);
%             if length(V(v).single_beat(i).beat)==1
%                 STD_v(v).std_beat_param(i,:) = zeros(1,20);
%             else
%                 STD_v(v).std_beat_param(i,:) = std(matr);
%             end
% 
%         end
%     end
% 
%     % average of parameters computed on single beat traces
%     for v=1:length(V)
%         for i = 1:length(V(v).single_beat)
% 
%             AvgParam(v).AvgParam(i).Baseline = mean([V(v).single_beat(i).beat.Baseline]);           
% 
%             AvgParam(v).AvgParam(i).Fmax_F0 = mean([V(v).single_beat(i).beat.DoverD0]);
% 
%             AvgParam(v).AvgParam(i).CD = mean([V(v).single_beat(i).beat.Mean_Calcium_CD]);
%             AvgParam(v).AvgParam(i).CD90 = mean([V(v).single_beat(i).beat.Mean_Calcium_CD90]);
%             AvgParam(v).AvgParam(i).CD50 = mean([V(v).single_beat(i).beat.Mean_Calcium_CD50]);
%             AvgParam(v).AvgParam(i).CD10 = mean([V(v).single_beat(i).beat.Mean_Calcium_CD10Max]);
%             AvgParam(v).AvgParam(i).Ton =  mean([V(v).single_beat(i).beat.Time_CalciumPeak]);
%             AvgParam(v).AvgParam(i).Toff =  mean([V(v).single_beat(i).beat.Time_CalciumRelax]);  
% 
%             AvgParam(v).AvgParam(i).T10on = mean([V(v).single_beat(i).beat.Time_to_10Contr]);
%             AvgParam(v).AvgParam(i).T50on = mean([V(v).single_beat(i).beat.Time_to_50Contr]);
%             AvgParam(v).AvgParam(i).T90on = mean([V(v).single_beat(i).beat.Time_to_90Contr]);
%             AvgParam(v).AvgParam(i).T10off = mean([V(v).single_beat(i).beat.Time_to_10_Relax]);
%             AvgParam(v).AvgParam(i).T50off = mean([V(v).single_beat(i).beat.Time_to_50_Relax]);
%             AvgParam(v).AvgParam(i).T90off = mean([V(v).single_beat(i).beat.Time_to_90_Relax]);
% 
%             AvgParam(v).AvgParam(i).tau = mean([V(v).single_beat(i).beat.tau]);
%             AvgParam(v).AvgParam(i).fit_a = mean([V(v).single_beat(i).beat.a]);
%             AvgParam(v).AvgParam(i).fit_c = mean([V(v).single_beat(i).beat.c]);
%             AvgParam(v).AvgParam(i).fit_rsquare = mean([V(v).single_beat(i).beat.rsquare]);
% 
%             AvgParam(v).AvgParam(i).t0 = mean([V(v).single_beat(i).beat.t0]);
%             AvgParam(v).AvgParam(i).tend = mean([V(v).single_beat(i).beat.tend]);
%         
%         end
%         
%     writetable(struct2table(AvgParam(v).AvgParam), fullfile(path,'Calcium_measurements_AvgParam.xlsx'),'Sheet',v );
%     
%     end
%         
% 
%     if exist('non_regular_beat')==1
%         nonRegBeat_all = [];
%         for k = 1:length(non_regular_beat)
%             nonRegBeat = [];
%             v = zeros(2,length(non_regular_beat(k).cell_trace));
%             v(1,:) = non_regular_beat(k).cell_trace;
%             v(2,1:length(non_regular_beat(k).beat_trace)) = non_regular_beat(k).beat_trace;
%             nonRegBeat = table([{'video ',num2str(non_regular_beat(k).video_n)};{'video ',num2str(non_regular_beat(k).video_n)}],[{'cell ',num2str(non_regular_beat(k).cell_n)};{'beat ',num2str(non_regular_beat(k).beat_n)}],v);
%             nonRegBeat_all = [nonRegBeat_all; nonRegBeat];
%         end
% 
%     writetable(nonRegBeat_all, fullfile(path,'Calcium_Traces_NonRegularBeats.xlsx'));
%     end
% 
%     warning off
%     
% 
%     for i = 1:length(STD_v)
%         writetable(table(STD_v(i).std_beat_param),fullfile(path,'Calcium_measurements_STDparam.xlsx'),'Sheet',i);
%     end
% 
% else
    
    V = video;
    [V.single_beat] = V.single_traces;
    V = rmfield(V,'single_traces');
    for i = 1:length(V)
        for j = 1:length(V(i).single_beat)
            [V(i).single_beat(j).data.trace] = V(i).single_beat(j).data.beats;
            V(i).single_beat(j).data = rmfield(V(i).single_beat(j).data,'beats');
        end
    end
% end

%% 
for i = 1:length(V)
    
    lengths = [];
    for j = 1:length(V(i).single_beat)
        for c = 1:length(V(i).single_beat(j).data)
            lengths(c,1) = length(V(i).single_beat(j).data(c).trace);       
        end
        min_length(j) = min(lengths);
    
        matrix = [];
        for c = 1:length(V(i).single_beat(j).data)
            V(i).single_beat(j).data(c).beats_cut  = V(i).single_beat(j).data(c).trace(1:min_length(j));
            matrix = [matrix; V(i).single_beat(j).data(c).beats_cut'];       
        end
        V(i).single_beat(j).mean(1,:) = mean(matrix,1);
    end
end


for i=1:length(V) 
    figure('Name',strcat('video ',num2str(i))),clf
    for j=1:length(V(i).single_beat) 
            subplot(5,ceil(length(V(i).single_beat)/5),j)
            for c = 1:length(V(i).single_beat(j).data) 
                t = (1:length(V(i).single_beat(j).data(c).beats_cut))/AcquisitionFrequency;
                plot(t,V(i).single_beat(j).data(c).beats_cut), hold on
                title(strcat('cell ',num2str(j)))
                xticklabels('')
                yticklabels('')
                grid on
            end
            tt = (1:length(V(i).single_beat(j).mean))/AcquisitionFrequency;
            plot(tt,V(i).single_beat(j).mean,'k')
    end
end


remove =questdlg('Would you like to remove any cell?','Filter undesired cells',...
    'Yes','No ', 'Yes');

if remove == 'Yes'
    prompt = {'\bf \fontsize{12} Please enter the video number for those videos that contain the cells to be removed (from figures title): e.g. [1 2 3] (if more than one cell from the same video, than repeat the video number)',...
    '\bf \fontsize{12} Please enter the cell number for the cells to be removed (subplots title): e.g. [1 2 3] (their order needs to match the videos order and when you want to remove more than 1 cell from the same video, report the cell number from the largest to the smallest'};
    dlgtitle = 'Remove cells';
    dims = [1 120];
    definput = {'','', '', '','','','',''};
    opts.Interpreter = 'tex';
    RC = inputdlg(prompt,dlgtitle,dims,definput, opts);
    rem_vid = str2num(RC{1});
    rem_cell = str2num(RC{2});
end

if remove == 'Yes'
    for i=1:length(rem_vid)
        V(rem_vid(i)).single_beat(rem_cell(i)) = [];
    end
end  

j=1;
id=[];
for i =1:length(V)
    if (isempty(V(i).single_beat)==1)
        id(j)=i;
        j=j+1;
    end
end
if isempty('id') == 0
    V(id) = [];
end


for i = 1:length(V)
    for j = 1:length(V(i).single_beat)
    [P,X] = findpeaks(V(i).single_beat(j).mean,'SortStr','descend');
    half_v = V(i).single_beat(j).mean(X(1):end);
    x = 1:length(half_v);
    CL = 1/PacingFrequency; %s
    n_time = 0.2*CL; % 20% of CL
    n = ceil(n_time*AcquisitionFrequency);    
    bsl = mean(half_v(end-n:end));
    bsln = ((ones (size(half_v,2),1) .*bsl))';
    [xint, yint] = intersections(x, half_v, x, bsln, ' robust');
    maxx = floor(min(xint))+ 3;
    if maxx > length (half_v)
    maxx = length (half_v);
    end
    half_vv = half_v(1:maxx);
    xx = 1:length(half_vv);
    %
    [fitresult, gof, xData, yData] = fit_tau(xx, half_vv);
    V(i).single_beat(j).tau = ((1/(-fitresult.b)).*(1/AcquisitionFrequency))*1000; %ms
    V(i).single_beat(j).a = fitresult.a;
    V(i).single_beat(j).c = fitresult.c;
    V(i).single_beat(j).rsquare = gof.rsquare;    
    end
end


warning('off');
b = 1;
measurements_table_complete = [];
for i=1:length(V)
    measurements_table = [];
    clearvars    Mean_Calcium_CD Mean_Calcium_CD90 Mean_Calcium_CD50 Mean_Calcium_CD10 Time_CalciumPeak...
        Time_CalciumRelax AVGCalcTrace DoverD0_values Baseline ...
        Time_to_10 Time_to_50 Time_to_90 Time_to_10Max_Relax ...
        Time_to_50_Relax Time_to_90_Relax mean_Calcium_magn tau ...
        fit_a fit_c fit_rsquare BR t0 tend Cell_id BBdistance Nperiod CD CD90 CD50 CD10 Ton...
        Toff  Fmax_F0 ...
        T10on T50on T90on T10off ...
        T50off T90off CaAmplitude
    bb = 1;
    for file = 1:length(V(i).single_beat)

        if AcquisitionFrequency>100
            YY = movmean(V(i).single_beat(file).mean,4);
        else
            YY = V(i).single_beat(file).mean;
        end
        
        [CalciumPeak, CTD, CTD90, CTD50, CTD10, Time_C_peak, T_C_relax, AvgCalciumTrace, DoverD0,...
            yAVGCalciumTrace_baseline, Time_to_90a, Time_to_50a, Time_to_10a, Time_to_10Relax,...
            Time_to_50Relax, Time_to_90Relax, CalciumMagnitude,CalciumTDPre,CalciumTDPost] = parameters(YY,AcquisitionFrequency);% parameters(movmean(V(i).single_beat(file).mean,5),AcquisitionFrequency);
   
        try
        CD(file,1) = CTD;
        CD90(file,1) = CTD90;
        CD50(file,1) = CTD50;
        CD10(file,1) = CTD10;
        Ton (file,1) =  Time_C_peak;
        Toff (file,1) =  T_C_relax;   
        AVGCalcTrace(file).Avg = V(i).single_beat(file).mean';%AvgCalciumTrace;   

        Fmax_F0(file,1) = DoverD0;
        Baseline(file,1) = yAVGCalciumTrace_baseline;

        T10on(file,1) = Time_to_90a;
        T50on(file,1) = Time_to_50a;
        T90on(file,1) =Time_to_10a;
        T10off(file,1) = Time_to_10Relax;
        T50off(file,1) = Time_to_50Relax;
        T90off(file,1) = Time_to_90Relax; 

        CaAmplitude(file,1) = CalciumMagnitude;

        tau(file,1) = V(i).single_beat(file).tau;
        fit_a(file,1) = V(i).single_beat(file).a;
        fit_c(file,1) = V(i).single_beat(file).c;
        fit_rsquare(file,1) = V(i).single_beat(file).rsquare;

        t0(file,1) = CalciumTDPre;
        tend(file,1) = CalciumTDPost;

        BR(file,1) = V(i).single_beat(file).BR;
        BBdistance(file,1) = V(i).single_beat(file).BBdistance* 1000;
        Nperiod(file,1) = V(i).single_beat(file).Nperiod;

        catch
            disp('>>>>>> error occurred in computing the parameters. Trace saved in "error_parameters" (when an error occurs param have value 0)')
            Er_param(b).data = V(i).single_beat(file).mean';
            Er_param(b).Video_id = i;
            Er_param(b).Video_name = V(i).name;
            Er_param(b).Cell_id = file;
            Cell_id(bb) = file;
            b = b+1;
            bb = bb + 1; 
            
            CD(file,1) = NaN;
            CD90(file,1) = NaN;
            CD50(file,1) = NaN;
            CD10(file,1) = NaN;
            Ton (file,1) =  NaN;
            Toff (file,1) =  NaN;   
            AVGCalcTrace(file).Avg = NaN;   
            Fmax_F0(file,1) = NaN;
            Baseline(file,1) = NaN;
            T10on(file,1) = NaN;
            T50on(file,1) = NaN;
            T90on(file,1) =NaN;
            T10off(file,1) = NaN;
            T50off(file,1) = NaN;
            T90off(file,1) = NaN; 
            CaAmplitude(file,1) = NaN;
            tau(file,1) = NaN;
            fit_a(file,1) = NaN;
            fit_c(file,1) = NaN;
            fit_rsquare(file,1) = NaN;
            t0(file,1) = NaN;
            tend(file,1) = NaN;
            BR(file,1) = NaN;
            BBdistance(file,1) = NaN;
            Nperiod(file,1) = NaN;
            SNR(file,1) = NaN;
        end 
    end
    
    if_error_param = exist('Er_param');
    

    K = 1:length(V(i).single_beat); 
    Cell_n = K';
    Video = cell(size(Cell_n));
    Video(:) = {V(i).name}; 

    CD =  CD * 1000;
    CD90 = CD90 * 1000;
    CD50 = CD50 * 1000; 
    CD10 = CD10 * 1000;
    Ton = Ton * 1000;
    Toff = Toff * 1000;
    T90on = T90on * 1000;
    T50on = T50on * 1000;
    T10on = T10on * 1000;
    T10off = T10off * 1000;
    T50off = T50off * 1000;
    T90off = T90off * 1000;
    t0 = t0 * 1000;
    tend = tend * 1000;

    measurements_table = table(Video, Cell_n,...
      Baseline, Fmax_F0, CD, CD90, CD50, ...
      CD10, Ton, Toff, ...
      T10on, T50on, T90on, T10off, ...
      T50off, T90off, tau, fit_a, fit_c,fit_rsquare, BR, t0, tend, BBdistance, Nperiod);%, BBtime_mean,Nperiod);

      if exist('Cell_id') == 1
        measurements_table([Cell_id],:) = [];
      end
      writetable(measurements_table, fullfile(path,'Calcium_measurements_forPlot.xlsx'),'Sheet',i );


    matrix = [];
    matrix = V(i).single_beat;
    SNR = [];
    SNRN = [];
    [SNR] = CalculateSNR(matrix, AcquisitionFrequency, cutoff, PacingFrequency);%
    SNRN = cell2mat({SNR.SignalToNoiseRatio})';
    Signal_to_Noise_Ratio = table(SNRN);

    Signal_to_Noise_Ratio.Properties.VariableNames = {'Signal_to_Noise_Ratio'};
    
    measurements_table = horzcat(measurements_table , Signal_to_Noise_Ratio);

    measurements_table_complete = [measurements_table_complete ; measurements_table];


end

writetable(measurements_table_complete, fullfile(path,'Calcium_measurements.xlsx') );

            
for m = 1:length(V)
    
    T = [];
    Celll = [];
    AVGCalcTrace_cut = [];
    for i = 1:length(V(m).single_beat)
        s = V(m).single_beat(i).mean';
        lengths(i) = length(s);
    end
    L = min(lengths);
    LL = max(lengths);
    for i = 1:length(V(m).single_beat)
        AVGCalcTrace_cut(i).Avg = V(m).single_beat(i).mean(1:L)';
        AVGCalcTrace_ext(i).Avg = [V(m).single_beat(i).mean';V(m).single_beat(i).mean(end)*ones(LL-length(V(m).single_beat(i).mean),1)];

        T = [T, AVGCalcTrace_ext(i).Avg];
    end

    if if_error_param == 1
        M = (m == [Er_param.Video_id]);
        if sum(M) >0 
            T(:,[Er_param(M).Cell_id]) = [];
        end
    end
    
    Celll = num2cell (T);
    writetable(cell2table(Celll),fullfile(path,'Calcium_Traces.xlsx'),'Sheet',video(m).name)

end


if if_skipped == 1
    [video_id,i1,i2] = unique([Sk.Video_id]);
    video_names = {Sk([i1]).Video_name}; %;
    for v = 1:length(video_id)
            video_data = Sk([Sk.Video_id]==video_id(v));
        for i=1:length(video_data)
            leng(i) = length(video_data(i).data);
        end
        leng_min = min(leng);
        leng_max = max(leng);
        Y = [];
        for i=1:length(video_data)
            ext = [];
            ext = [video_data(i).data;video_data(i).data(end)*ones(leng_max-length(video_data(i).data),1)];
            Y = [Y ext];
        end
        writematrix(Y,fullfile(path,'Calcium_Traces_noisy.xlsx'),'Sheet',video_names{1,v})
    end
end


video_id = [];
video_data = [];

            
if if_error == 1
    [video_id,i3,i4] = unique([Er.Video_id]);
    video_names = {Er([i3]).Video_name}; %;

    for v = 1:length(video_id)
        video_data = Er([Er.Video_id]==video_id(v));
    for i=1:length(video_data)
        leng(i) = length(video_data(i).data);
    end
    leng_min = min(leng);
    leng_max = max(leng);

    Y = [];
    for i=1:length(video_data)
        ext = [];
        ext = [video_data(i).data;video_data(i).data(end)*ones(leng_max-length(video_data(i).data),1)];
        Y = [Y ext];
    end
    writematrix(Y,fullfile(path,'Calcium_Traces_errors.xlsx'),'Sheet',video_names{1,v})
    end
end


video_id = [];
video_data = [];

        
if if_error_param == 1
    [video_id,j1,j2] = unique([Er_param.Video_id]);
    video_names = {Er_param([j1]).Video_name}; %;

    for v = 1:length(video_id)
        video_data = Er_param([Er_param.Video_id]==video_id(v));
    for i=1:length(video_data)
        leng(i) = length(video_data(i).data);
    end
    leng_min = min(leng);
    leng_max = max(leng);

    Y = [];
    for i=1:length(video_data)
        ext = [];
        ext = [video_data(i).data;video_data(i).data(end)*ones(leng_max-length(video_data(i).data),1)];
        Y = [Y ext];
    end
    writematrix(Y,fullfile(path,'Calcium_Traces_error_param.xlsx'),'Sheet',video_names{1,v})
    end
end


video_id = [];
video_data = [];

            
if if_extra_beat == 1
    [video_id,i5,i6] = unique([extra_beat.Video_id]);
    video_names = {extra_beat([i5]).Video_name}; %;

    for v = 1:length(video_id)
        video_data = extra_beat([extra_beat.Video_id]==video_id(v));
    for i=1:length(video_data)
        leng(i) = length(video_data(i).data);
    end
    leng_min = min(leng);
    leng_max = max(leng);

    Y = [];
    for i=1:length(video_data)
        ext = [];
        ext = [video_data(i).data;video_data(i).data(end)*ones(leng_max-length(video_data(i).data),1)];
        Y = [Y ext];
    end
    writematrix(Y,fullfile(path,'Calcium_Traces_extra_beat.xlsx'),'Sheet',video_names{1,v})
    end
    EB_table = table(s_names([extra_beat.Video_id]),[extra_beat.BBdistance]',[extra_beat.Nperiod]','VariableNames',{'video','BBdistance','Nperiod'});
    writetable(EB_table,fullfile(path,'extra_beat_measurements.xlsx'))
end


%% Kymographs

if exist('Er_param')
    rem_vid = [Er_param.Video_id];
    rem_cell = [Er_param.Cell_id];

    for i=1:length(rem_vid)
        pair(i,1:2)=[rem_vid(i) rem_cell(i)];
    end   
    else
        pair = [0 0];
end
    

mkdir (path,'Kymographs')
kpath = [path,'/Kymographs/'];
for i = 1:length(video)
    for  l = 1:length(video(i).single_traces)
        il=[i l];
        cut = il ==pair;
        cut2 = sum(cut,2);
        cut3 = (cut2 == 2);
        if sum(cut3) == 1

        else
            tr = {video(i).single_traces(l).data.beats}.';
            FT = cat(1,tr{:});
            FullTrace = (FT-min(FT))/(max((FT-min(FT))));
            sizey = round((size(FullTrace,1))/7.5);
            Kymo = ones(sizey,(size(FullTrace,1))) .* FullTrace';

            grayKymo = uint8(floor(Kymo* 255));
            IndexedKymo = ind2rgb(grayKymo, autumn);

            FileName = [kpath,video(i).name{1,1},' Cell', ' ',num2str(l) , '.jpeg'];

            imwrite(IndexedKymo, FileName, 'jpeg');
        end
    end
end



%%
fprintf (1, '>>>   Times are given in msec.  \n')
fprintf (1, '>>>   Parameters saved in excel spreadsheet "Calcium_measurements" in the folder "calcium_analysis".  \n')
fprintf (1, '>>>   Calcium Traces saved in excel spreadsheet "Calcium_Traces" in the folder "calcium_analysis".  \n')
fprintf (1, '>>>   Calcium Traces not analysed because noisy saved in excel spreadsheet "Calcium_Traces_noisy".  \n')
fprintf (1, '>>>   Calcium Traces not analysed because of an error occurred saved in excel spreadsheet "Calcium_Traces_errors".  \n')
fprintf (1, '>>>   Calcium Traces not analysed because of an error occurred when computing parameters saved in excel spreadsheet "Calcium_Traces_error_param".  \n')
fprintf (1, '>>>   Calcium Traces with extra beat(s) saved in excel spreadsheet "Calcium_Traces_extra_beat".  \n')
fprintf (1, '>>>   Kymographs saved in folder "Kymographs".  \n')


%%

function [fitresult, gof,xData, yData] = fit_tau(x, half_v)

[xData, yData] = prepareCurveData( x, half_v );

ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf 0];
opts.Upper = [Inf Inf Inf 0];

[fitresult, gof] = fit( xData, yData, ft, opts );


end

