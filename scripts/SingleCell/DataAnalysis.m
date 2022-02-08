
% DataAnalysis_SingleCell.m

path = cd;
if quantitative_data == 'y'
    calcium_file = 'Calcium_Traces_quantitative.mat';
    load(fullfile(path,calcium_file))
    Calcium_Traces = Calcium_Traces_quantitative;
else
    calcium_file = 'Calcium_Traces.mat';
    load(fullfile(path,calcium_file))
end

global AcquisitionFrequency
global cutoff
global PacingFrequency

figure('Name', 'Selected frames'), clf
for i=1:length(Calcium_Traces)
    subplot(5,ceil(length(Calcium_Traces)/5),i)
    t = (1:length(Calcium_Traces(i).data));
    plot(t,Calcium_Traces(i).data), hold on, grid on
    if M==0
        plot(t(N+1:end),Calcium_Traces(i).data(N+1:end))
        Calcium_Traces_analysis(i).data = Calcium_Traces(i).data(N+1:end);
    else
        plot(t(N+1:M),Calcium_Traces(i).data(N+1:M))
        Calcium_Traces_analysis(i).data = Calcium_Traces(i).data(N+1:M);
    end
    Calcium_Traces_analysis(i).name = Calcium_Traces(i).name;
    title(strcat('cell',num2str(i)))
    yticklabels(' ')
    xticklabels(' ')
end

remove =questdlg('Based on the figure, would you like to remove any cell?','Filter undesired cells',...
    'Yes','No ', 'Yes');

switch remove
    case 'Yes'
    prompt = {'\bf \fontsize{12} Please enter the cell number as shown in the sub-plot title (eg [1 2 3]):',...
    };
    dlgtitle = 'Remove cells';
    dims = [1 88];
    definput = {'','', '', '','','','',''};
    opts.Interpreter = 'tex';
    RC = [];
    RC = inputdlg(prompt,dlgtitle,dims,definput, opts);
    rem_cell = str2num(RC{1});

    Calcium_Traces_analysis(rem_cell)=[];
end

Names = {Calcium_Traces_analysis.name}.';

%%
photobl =questdlg('Do you need to apply photo bleach correction?','Photobleach correction toggle',...
    'Yes','No ', 'Yes');

single_traces = [];

if segmentation == 'n'
    [single_traces, skipped, extra_beat, errors] = BeatSegmentation({Calcium_Traces_analysis.name}.',{Calcium_Traces_analysis.data}.',param_single_beat);
elseif segmentation == 'y'
    T0 = xlsread('T0.xlsx');
    [single_traces, skipped, extra_beat, errors] = BeatSegmentation2({Calcium_Traces_analysis.name}.',{Calcium_Traces_analysis.data}.',T0);
end

if isempty(skipped)
    clearvars skipped
end
if isempty(extra_beat)
    clearvars extra_beat
end
if isempty(errors)
    clearvars errors
end


switch photobl
    case 'Yes'
    
    [single_traces] = PhotoBleachCorrection(single_traces);
    
    if segmentation == 'n'
        [single_traces_phbl, skipped_phbl, extra_beat_phbl, errors_phbl] = BeatSegmentation({single_traces.name}.',{single_traces.BaselineCorrectedTraces}.',param_single_beat);
    elseif segmentation == 'y'
        T02 = T0;
        T02(:,2) = ones(size(T02(:,2)));
        [single_traces_phbl, skipped_phbl, extra_beat_phbl, errors_phbl] = BeatSegmentation2({single_traces.name}.',{single_traces.BaselineCorrectedTraces}.',T02);
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


    figure('Name','Baseline Correction'), clf
    for i=1:length(single_traces)
        subplot(5,ceil(length(single_traces)/5),i)

        selected_trace = [];

        cll = struct2cell(single_traces(i).data);
        sqcll = squeeze(cll);
        selected_trace = cat(1,sqcll{:});
        t = 1:length(selected_trace);     
        plot(t,selected_trace), hold on, grid on

        tt = 1:length(single_traces(i).BaselineCorrectedTraces);
        plot(tt,single_traces(i).BaselineCorrectedTraces)

        title(strcat('cell',num2str(i)))
        yticklabels(' ')
        xticklabels(' ')
    end
    
    single_traces_before_phbl = single_traces;
    single_traces = [];
    single_traces = single_traces_phbl;
    
    remove =questdlg('Based on the figure, would you like to remove any cell?','Filter undesired averaged traces',...
    'Yes','No ', 'Yes');
    switch remove
        case 'Yes'
            prompt = {'\bf \fontsize{12} Please enter the cell number as shown in the sub-plot title (eg [1 2 3]):',...
            };
            dlgtitle = 'Remove cells';
            dims = [1 88];
            definput = {'','', '', '','','','',''};
            opts.Interpreter = 'tex';
            RC = [];
            RC = inputdlg(prompt,dlgtitle,dims,definput, opts);
            rem_cell = str2num(RC{1});

            single_traces(rem_cell)=[];
    end
    names = {single_traces.name}';

end

if_skipped = exist('skipped');
if_extra_beat = exist('extra_beat');
if_error = exist('errors');


%%
tot_cell_number = length(single_traces);

if param_single_beat == 'y'
    
    k = 1;
   %% 
    for i=1:length(single_traces)
        j = 1;
        jjj = 1;
        for jj=1:length(single_traces(i).data)
            
            try
            yy = single_traces(i).data(jj).beats; 
            
            % fit
            [PP,XX] = findpeaks(yy,'SortStr','descend');
            half_vv = yy(XX(1):end);
            xx = 1:length(half_vv);
            CL = 1/PacingFrequency; 
            n_time = 0.2*CL;
            n = ceil(n_time*AcquisitionFrequency);
            bsl = mean(half_vv(end-n:end));
            bsln = ((ones (size(half_vv,1),1) .*bsl))';
            [xint, yint] = intersections(xx, half_vv, xx, bsln, ' robust');
            maxx = floor(min(xint))+ 3;
            if maxx > length (half_vv)
                maxx = length (half_vv);
            end
            half_vv2 = half_vv(1:maxx);
            xx2 = 1:length(half_vv2);
            
            half_vv3 = half_vv(maxx+1:end);
            xx3 = 1:length(half_vv3);
            
            [p_EAD] = findpeaks(half_vv2,'MinPeakProminence',0.2*(max(half_vv2-min(half_vv2))));                       
           
            [fitresult2, gof2, xData2, yData2] = fit_tau(xx2, half_vv2');
            Sbeat(i).beat(jjj).tau = ((1/(-fitresult2.b)).*(1/AcquisitionFrequency))*1000; %ms
            Sbeat(i).beat(jjj).a = fitresult2.a;
            Sbeat(i).beat(jjj).c = fitresult2.c;
            Sbeat(i).beat(jjj).rsquare = gof2.rsquare; 
            
            jjj = jjj+1;

            pks = [];
            % baseline check
            if length(half_vv3) > 3
                
                warning off
                if AcquisitionFrequency > 100
                    [pks,locs,w,p] = findpeaks(movmean(half_vv3,10),'MinPeakHeight',bsl+0.1*bsl);  
                else
                    [pks,locs,w,p] = findpeaks(half_vv3,'MinPeakProminence',0.03*mean(half_vv3));
                end
                
               
            end
            
            if ((isempty(p_EAD) == 0) || (isempty(pks) == 0))
                non_regular_beat(k).cell_n = i;
                non_regular_beat(k).beat_n = jj;
                non_regular_beat(k).cell_trace = single_traces(i).original;
                non_regular_beat(k).beat_trace = yy;
                non_regular_beat(k).EAD = 0;
                non_regular_beat(k).DAD = 0;
                non_regular_beat(k).Nbeats = single_traces(i).Nbeats;
                if (isempty(p_EAD) == 0)
                    non_regular_beat(k).EAD = 1;
                end
                if (isempty(pks) == 0)
                    non_regular_beat(k).DAD = 1;
                end
                k = k+1;

               
            else
                
                if AcquisitionFrequency>100
                    Y = movmean(yy,4);
                else
                    Y = yy;
                end
                
                [CalciumPeak, CTD, CTD90, CTD50, CTD10, Time_C_peak, T_C_relax, AvgCalciumTrace, DoverD0,...
                yAVGCalciumTrace_baseline, Time_to_90a, Time_to_50a, Time_to_10a, Time_to_10Relax,...
                Time_to_50Relax, Time_to_90Relax, CalciumMagnitude,CalciumTDPre,CalciumTDPost] = parameters(Y',AcquisitionFrequency);

                p = [CalciumPeak, CTD, CTD90, CTD50, CTD10, Time_C_peak, T_C_relax, DoverD0,...
                yAVGCalciumTrace_baseline, Time_to_90a, Time_to_50a, Time_to_10a, Time_to_10Relax,...
                Time_to_50Relax, Time_to_90Relax, CalciumMagnitude,CalciumTDPre,CalciumTDPost];

                if length(p) == 18 
                    single_beat(i).beat(j).Fmax_F0 = DoverD0;
                    single_beat(i).beat(j).Baseline = yAVGCalciumTrace_baseline;
                    single_beat(i).beat(j).CaAmplitude = CalciumMagnitude;
                    single_beat(i).beat(j).CD =  CTD * 1000;
                    single_beat(i).beat(j).CD90 = CTD90 * 1000;
                    single_beat(i).beat(j).CD50 = CTD50 * 1000; 
                    single_beat(i).beat(j).CD10 = CTD10 * 1000;
                    single_beat(i).beat(j).Ton = Time_C_peak * 1000;
                    single_beat(i).beat(j).Toff = T_C_relax * 1000;
                    single_beat(i).beat(j).T10on = Time_to_90a * 1000;
                    single_beat(i).beat(j).T50on = Time_to_50a * 1000;
                    single_beat(i).beat(j).T90on = Time_to_10a * 1000;
                    single_beat(i).beat(j).T10off = Time_to_10Relax * 1000;
                    single_beat(i).beat(j).T50off = Time_to_50Relax * 1000;
                    single_beat(i).beat(j).T90off = Time_to_90Relax * 1000;
                    single_beat(i).beat(j).t0 = CalciumTDPre * 1000;
                    single_beat(i).beat(j).tend = CalciumTDPost * 1000;
                    single_beat(i).beat(j).tau = ((1/(-fitresult2.b)).*(1/AcquisitionFrequency))*1000; %ms
                    single_beat(i).beat(j).a = fitresult2.a;
                    single_beat(i).beat(j).c = fitresult2.c;
                    single_beat(i).beat(j).rsquare = gof2.rsquare;
                    single_beat(i).name = single_traces(i).name;
                    single_beat(i).traces(j).data = yy;
                    
                    j=j+1;
                end
                
            end
            
            catch
                
            end
            
        end
        
    end
    %%
    j=1;
    ID = [];
    if exist('single_beat')
        for i =1:length(single_beat)
            if isempty(single_beat(i).beat)
                ID(j) = i;
                j=j+1;
            end
        end
        single_beat(ID) = [];
    end

    
    if exist('single_beat')
        for i = 1:length(single_beat)
            matr = [];
            matr = [[single_beat(i).beat.Baseline].', [single_beat(i).beat.Fmax_F0].'...
             [single_beat(i).beat.CD].',...
            [single_beat(i).beat.CD90].', [single_beat(i).beat.CD50].',...
            [single_beat(i).beat.CD10].', [single_beat(i).beat.Ton].',...
            [single_beat(i).beat.Toff].', [single_beat(i).beat.T10on].',...
            [single_beat(i).beat.T50on].', [single_beat(i).beat.T90on].',...
            [single_beat(i).beat.T10off].', [single_beat(i).beat.T50off].',...
            [single_beat(i).beat.T90off].',...
            [single_beat(i).beat.tau].', [single_beat(i).beat.a].', [single_beat(i).beat.c].', [single_beat(i).beat.rsquare].',...
            [single_beat(i).beat.t0].', [single_beat(i).beat.tend].'];

            [N_beat(i),C] = size(matr);
            std_beat_param(i,:) = std(matr,0,1);
            
            writetable(struct2table(single_beat(i).beat ), fullfile(path,'Calcium_measurements_IndividualBeat.xlsx'),'Sheet',single_beat(i).name);

            for j = 1:length(single_beat(i).traces)
                len(j) = length(single_beat(i).traces(j).data);
            end
            for j = 1:length(single_beat(i).traces)
               single_beat_traces(i).data(:,j) =  single_beat(i).traces(j).data(1:min(len));
            end
            
            writematrix(single_beat_traces(i).data, fullfile(path,'Calcium_Traces_IndividualBeat.xlsx'),'Sheet',single_beat(i).name);

            
        end

        for i=1:length(single_beat)
            
            AvgParam(i).Baseline = mean([single_beat(i).beat.Baseline]);           
            AvgParam(i).Fmax_F0 = mean([single_beat(i).beat.Fmax_F0]);
            AvgParam(i).CD = mean([single_beat(i).beat.CD]);
            AvgParam(i).CD90 = mean([single_beat(i).beat.CD90]);
            AvgParam(i).CD50 = mean([single_beat(i).beat.CD50]);
            AvgParam(i).CD10 = mean([single_beat(i).beat.CD10]);
            AvgParam(i).Ton =  mean([single_beat(i).beat.Ton]);
            AvgParam(i).Toff =  mean([single_beat(i).beat.Toff]);  
            AvgParam(i).T10on = mean([single_beat(i).beat.T10on]);
            AvgParam(i).T50on = mean([single_beat(i).beat.T50on]);
            AvgParam(i).T90on = mean([single_beat(i).beat.T90on]);
            AvgParam(i).T10off = mean([single_beat(i).beat.T10off]);
            AvgParam(i).T50off = mean([single_beat(i).beat.T50off]);
            AvgParam(i).T90off = mean([single_beat(i).beat.T90off]);
            AvgParam(i).tau = mean([single_beat(i).beat.tau]);
            AvgParam(i).fit_a = mean([single_beat(i).beat.a]);
            AvgParam(i).fit_c = mean([single_beat(i).beat.c]);
            AvgParam(i).fit_rsquare = mean([single_beat(i).beat.rsquare]);
            AvgParam(i).t0 = mean([single_beat(i).beat.t0]);
            AvgParam(i).tend = mean([single_beat(i).beat.tend]);

        end  

        writetable([table((1:length(AvgParam))'),struct2table( AvgParam)], fullfile(path,'Calcium_measurements_AvgParam.xlsx'),'Sheet',1 );
        writetable(table({single_beat.name}.'), fullfile(path,'Calcium_measurements_AvgParam.xlsx'),'Sheet',1 );

        s2 = array2table([N_beat' std_beat_param],'VariableNames',[{'N_beat'};fieldnames(AvgParam)]);
        writetable(s2, fullfile(path,'Calcium_measurements_AvgParam.xlsx'),'Sheet',2 );

        warning off
        
    end

    if exist('non_regular_beat')==1
        nonRegBeat_all = [];
        for K = 1:length(non_regular_beat)
            cell_traces_len(K) = length(non_regular_beat(K).cell_trace);
        end
        for K = 1:length(non_regular_beat)
            nonRegBeat = [];
            v = zeros(2,max(cell_traces_len));%length(non_regular_beat(K).cell_trace));
            v(1,1:length(non_regular_beat(K).cell_trace)) = non_regular_beat(K).cell_trace;
            v(2,1:length(non_regular_beat(K).beat_trace)) = non_regular_beat(K).beat_trace;
            nonRegBeat = table([[{'cell ',num2str(non_regular_beat(K).cell_n)};...
                {'beat ',num2str(non_regular_beat(K).beat_n)}],[{'EADs ',num2str(non_regular_beat(K).EAD)};...
                {'DADs ',num2str(non_regular_beat(K).DAD)}]],v);
            nonRegBeat_all = [nonRegBeat_all; nonRegBeat];
        end

        writetable(nonRegBeat_all, fullfile(path,'Calcium_Traces_NonRegularBeats.xlsx'));
    end
    
end

if exist('non_regular_beat')==1
    irregular_cell_number = length(unique([non_regular_beat.cell_n]));
    regular_cell_number = tot_cell_number - irregular_cell_number;
else
    regular_cell_number = tot_cell_number;
end
    

%%

if exist('non_regular_beat')
    c = [non_regular_beat.cell_n];
    b = [non_regular_beat.beat_n];
else 
    c = [ ];
    b = [ ];   
end

for i = 1:length(single_traces)
    bj = c == i;
    bbj = b(bj);
   
    lengths = [];
    for j = 1:length(single_traces(i).data)
        lengths(j,1) = length(single_traces(i).data(j).beats);
    end
    min_length(i) = min(lengths);
    
    matrix = [];
    for j = 1:length(single_traces(i).data)
        
        if ismember(j,bbj)
            
        else
            single_traces(i).data(j).beats_cut  = single_traces(i).data(j).beats(1:min_length(i));
            matrix = [matrix; single_traces(i).data(j).beats_cut']; 
        end
        
    end
    single_traces(i).mean(1,:) = mean(matrix,1);
    
end


figure('Name','Segmented Beats'), clf
for i=1:length(single_traces)
    subplot(5,ceil(length(single_traces)/5),i)
    for j=1:length(single_traces(i).data)
        if exist('single_traces(i).data(j).beats_cut')
        plot(single_traces(i).data(j).beats_cut), hold on
        xticklabels('')
        yticklabels('')
        end
    end
    plot(single_traces(i).mean,'k'), hold on

    title(strcat('cell',num2str(i)))
end

remove =questdlg('Based on the figure, would you like to remove any cell?','Filter undesired averaged traces',...
'Yes','No ', 'Yes');
switch remove
    case 'Yes'
        prompt = {'\bf \fontsize{12} Please enter the cell number as shown in the sub-plot title (eg [1 2 3]):',...
        };
        dlgtitle = 'Remove cells';
        dims = [1 88];
        definput = {'','', '', '','','','',''};
        opts.Interpreter = 'tex';
        RC = [];
        RC = inputdlg(prompt,dlgtitle,dims,definput, opts);
        rem_cell = str2num(RC{1});

        single_traces(rem_cell)=[];
end
names = {single_traces.name}';

%%
for i = 1:length(single_traces)
    if isempty(single_traces(i).mean)
        
    else
    [P,X] = findpeaks(single_traces(i).mean,'SortStr','descend');
    half_v = single_traces(i).mean(X(1):end);
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
    [fitresult, gof, xData, yData] = fit_tau(x, half_v); %(xx, half_vv)
    single_traces(i).tau = ((1/(-fitresult.b)).*(1/AcquisitionFrequency))*1000; %ms
    single_traces(i).a = fitresult.a;
    single_traces(i).c = fitresult.c;
    single_traces(i).rsquare = gof.rsquare;
    end
    
end

%%
b = 1;
if isempty(single_traces) == 0

    for file = 1:length(single_traces)
        atLeast1Beat = 0;
        if isempty(single_traces(file).mean) == 0
           atLeast1Beat = 1; 
            if AcquisitionFrequency>100
                YY = movmean(single_traces(file).mean,4);
            else
                YY = single_traces(file).mean;
            end

        [CalciumPeak, CTD, CTD90, CTD50, CTD10, Time_C_peak, T_C_relax, AvgCalciumTrace, DoverD0,...
            yAVGCalciumTrace_baseline, Time_to_90a, Time_to_50a, Time_to_10a, Time_to_10Relax,...
            Time_to_50Relax, Time_to_90Relax, CalciumMagnitude,CalciumTDPre,CalciumTDPost] = parameters_opposite(YY,AcquisitionFrequency);

        try
        CD(file,1) = CTD;
        CD90(file,1) = CTD90;
        CD50(file,1) = CTD50;
        CD10(file,1) = CTD10;
        Ton (file,1) =  Time_C_peak;
        Toff (file,1) =  T_C_relax;   
        AVGCalcTrace(file).Avg = single_traces(file).mean';   

        Fmax_F0(file,1) = DoverD0;
        Baseline(file,1) = yAVGCalciumTrace_baseline;

        T10on(file,1) = Time_to_90a;
        T50on(file,1) = Time_to_50a;
        T90on(file,1) =Time_to_10a;
        T10off(file,1) = Time_to_10Relax;
        T50off(file,1) = Time_to_50Relax;
        T90off(file,1) = Time_to_90Relax; 

        CaAmplitude(file,1) = CalciumMagnitude;

        tau(file,1) = single_traces(file).tau;
        fit_a(file,1) = single_traces(file).a;
        fit_c(file,1) = single_traces(file).c;
        fit_rsquare(file,1) = single_traces(file).rsquare;

        t0(file,1) = CalciumTDPre;
        tend(file,1) = CalciumTDPost;

        BR(file,1) = single_traces(file).BR;
        BBdistance(file,1) = single_traces(file).BBdistance* 1000;
        NPeriod(file,1) = single_traces(file).Nperiod;

        catch
            disp('>>>>>> error occurred in computing the parameters. Trace saved in "error_parameters"')
            error_parameters(b).data = single_traces(file).mean';
            error_parameters(b).name = single_traces(file).name;
            error_parameters(b).id = file;
            b = b+1;

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
            NPeriod(file,1) = NaN;
            SNR(file,1) = NaN;

        end

        end

    end
else
    atLeast1Beat = 0;
end
%%
if_error_param = exist('error_parameters');

BBdistance_eb = zeros(length(single_traces),1);
NPeriod_eb = zeros(length(single_traces),1);

if atLeast1Beat == 1
    
Cell_n = 1:length(single_traces); 
Cell_n = Cell_n';

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

measurements_table = table( Cell_n,...
  Baseline, Fmax_F0, CD, CD90, CD50, ...
  CD10, Ton, Toff, ...
  T10on, T50on, T90on, T10off, ...
  T50off, T90off, tau, fit_a, fit_c,fit_rsquare,BR, t0, tend, BBdistance, NPeriod, BBdistance_eb,NPeriod_eb,CaAmplitude);


[SNR] = CalculateSC_SNR(single_traces, AcquisitionFrequency, cutoff, PacingFrequency); 
SNRN = cell2mat({SNR.SignalToNoiseRatio})';
Selection = measurements_table.Cell_n';
Signal_to_Noise_Ratio = table(SNRN);
Signal_to_Noise_Ratio.Properties.VariableNames = {'Signal_to_Noise_Ratio'};
measurements_table = horzcat(measurements_table , Signal_to_Noise_Ratio);

warning('off');

if exist('error_parameters')
    measurements_table([error_parameters.id],:) = [];
    AVGCalcTrace([error_parameters.id]) = [];
    names ([error_parameters.id]) = [];
    single_traces_with_errorparam = single_traces;
    single_traces([error_parameters.id])=[];
end


if (exist('extra_beat')) == 1
    extra_beat_marix = zeros(length(extra_beat),size(measurements_table,2));
    extra_beat_marix(:,end-2) = [extra_beat.BBdistance]';
    extra_beat_marix(:,end-1) = [extra_beat.Nperiod]';
    extra_beats_tab = array2table(extra_beat_marix);
    extra_beats_tab.Properties.VariableNames = measurements_table.Properties.VariableNames;

    writetable(vertcat(measurements_table,extra_beats_tab), fullfile(path,'Calcium_measurements.xlsx') );
    FileNames_table = table(vertcat(names,extra_beat.name));
    writetable(FileNames_table, fullfile(path,'Calcium_measurements.xlsx') );

else
    writetable(measurements_table, fullfile(path,'Calcium_measurements.xlsx') );
    FileNames_table = table(names);
    writetable(FileNames_table, fullfile(path,'Calcium_measurements.xlsx') );
end


T = [];
lengths = [];
for i = 1:length(AVGCalcTrace)
    lengths(i) = length(AVGCalcTrace(i).Avg);
end
L = min(lengths);
LL = max(lengths);
for i = 1:length(AVGCalcTrace)
    AVGCalcTrace_cut(i).Avg = AVGCalcTrace(i).Avg(1:L);
    AVGCalcTrace_ext(i).Avg = [AVGCalcTrace(i).Avg;AVGCalcTrace(i).Avg(end)*ones(LL-length(AVGCalcTrace(i).Avg),1)];

    T = [T, AVGCalcTrace_ext(i).Avg];

end

TT = num2cell (T);
Tab = vertcat(names',TT);
writetable(cell2table(Tab),fullfile(path,'Calcium_Traces.xlsx'))

end

if if_extra_beat == 1
    for i=1:length(extra_beat)
        Tab2 = vertcat(extra_beat(i).name,num2cell(extra_beat(i).data));
        writetable(cell2table(Tab2),fullfile(path,'Calcium_Traces_extra_beat.xlsx'),'Sheet',i)
    end
end

if if_skipped == 1
    for i=1:length(skipped)
        Tab3 = vertcat(skipped(i).name,num2cell(skipped(i).data));
        writetable(cell2table(Tab3),fullfile(path,'Calcium_Traces_noisy.xlsx'),'Sheet',i)
    end
end


if if_error == 1
    for i=1:length(errors)
        Tab5 = vertcat(errors(i).name,num2cell(errors(i).data));
        writetable(cell2table(Tab5),fullfile(path,'Calcium_Traces_errors.xlsx'),'Sheet',i)
    end
end

if if_error_param == 1
    regular_cell_number = regular_cell_number - length(error_parameters);
    for i=1:length(error_parameters)
        Tab6 = vertcat(error_parameters(i).name,num2cell(error_parameters(i).data));
        writetable(cell2table(Tab6),fullfile(path,'Calcium_Traces_error_param.xlsx'),'Sheet',i)
    end
end
warning('off');


%% Kymographs

if AcquisitionFrequency>100
    
else
    mkdir (path,'Kymographs')
    kpath = [path,'/Kymographs/'];

    for  l = 1:length(single_traces)

        tr = {single_traces(l).data.beats}.';
        FT = cat(1,tr{:});
        FullTrace = (FT-min(FT))/(max((FT-min(FT))));
        sizey = round((size(FullTrace,1))/7.5);
        Kymo = ones(sizey,(size(FullTrace,1))) .* FullTrace';

        grayKymo = uint8(floor(Kymo* 255));
        IndexedKymo = ind2rgb(grayKymo, autumn);
        FileName = [kpath,single_traces(l).name, '.jpeg'];

        imwrite(IndexedKymo, FileName, 'jpeg');

    end
    
end

fprintf (1, '>>>   Calcium tiff stack used to extract parameters.  \n')
fprintf (1, '>>>   Times are given in msec.  \n')
fprintf (1, '>>>   Parameters saved in excel spreadsheet "Calcium_measurements" in the folder "calcium_analysis".  \n')
fprintf (1, '>>>   Calcium Traces saved in excel spreadsheet "Calcium_Traces" in the folder "calcium_analysis".  \n')
fprintf (1, '>>>   Calcium Traces with extra beat(s) saved in excel spreadsheet "Calcium_Traces_extra_beat".  \n')
fprintf (1, '>>>   Calcium Traces with irregular peaks saved in excel spreadsheet "Calcium measurements NonRegularBeats".  \n')
fprintf (1, '>>>   Calcium Traces not analysed because noisy saved in excel spreadsheet "Calcium_Traces_noisy".  \n')
fprintf (1, '>>>   Calcium Traces not analysed because of an error occurred saved in excel spreadsheet "Calcium_Traces_errors".  \n')
fprintf (1, '>>>   Calcium Traces not analysed because of an error occurred when computing parameters saved in excel spreadsheet "Calcium_Traces_error_param".  \n')
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

