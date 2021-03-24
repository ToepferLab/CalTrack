
clearvars -except PacingFrequency AcquisitionFrequency main_folder cutoff scripts_folder atLeast1Beat non_regular_beat regular_cell_number irrregular_cell_number
res_folder = cd;

disp('Non Regular Beats Analysis:')

if exist('Calcium_Traces.xlsx')
    [ca_traces] = xlsread('Calcium_Traces.xlsx');
%     N_cell = size(ca_traces,2);
else
%     N_cell = 0;
end

Number_of_cell_extra_beat = 0;
Number_of_cell_non_regular_decay = 0;


if_extra_beat = exist ('Calcium_Traces_extra_beat.xlsx', 'file');
if_non_regular_decay = exist ('Calcium_Traces_NonRegularBeats.xlsx', 'file');

%%
if if_extra_beat == 2
    figure('Name','Calcium Traces - extra beat(s)')
    s_names_eb = sheetnames('Calcium_Traces_extra_beat.xlsx');
    
    for i=1:length(s_names_eb)
        
        subplot(4,ceil(length(s_names_eb)/4),i)
        [ca_traces_extra_beat,ca_traces_txt_extra_beat,ca_traces_raw_extra_beat] = xlsread('Calcium_Traces_extra_beat.xlsx',i);
        ca_eb(i).data = ca_traces_extra_beat;
        x = [1:size(ca_traces_extra_beat,1)]*(1/AcquisitionFrequency)*1000;
        plot(x,ca_traces_extra_beat)
        xlabel('Time (ms)')
        grid on
        title(strcat('Process',ca_traces_txt_extra_beat{2, 1}(9:end)))

    end

    Number_of_cell_extra_beat = length(s_names_eb)

    cd(scripts_folder)
    
    
    z = 1;
    BBdistance = [];

    CL = 1/PacingFrequency; %s
    offset_time = 0.1*CL;

    for i = 1:length(ca_eb)

        if AcquisitionFrequency>100
            y_original = ca_eb(i).data;
            DT = round((1/100)*AcquisitionFrequency);
            y = movmean(y_original(1:DT:end),8);
            AF = 1/((1/AcquisitionFrequency)*DT);
        else
            y = ca_eb(i).data;
            y_original = y;
            AF = AcquisitionFrequency;
            DT = 1;
        end
        offset = ceil(offset_time*AF);


        [~,~,~,pr] = findpeaks(diff(y)); 
        [peaks,location,~,~] = findpeaks(diff(y),'MinPeakProminence',0.5*max(pr)); 

%         y = y(location(1)-offset:end); 
%         IDD = location(1)-offset;
%         y_original = y_original(IDD+(IDD-1)*(DT-1):end);

        [p2,l2,w2,pr2] = findpeaks(y);
        [peaks2,location2,width2,prominence2] = findpeaks(y,'MinPeakProminence',0.5*max(pr2));

        peaks_time = location2/AcquisitionFrequency;

        BBdistance(i).dt = (diff(peaks_time));
        BBdistance(i).extra = find(BBdistance(i).dt<cutoff);  

%         peaks = [];
%         location = [];
%         [p,l,w,pr] = findpeaks(diff(y));
%         [peaks,location,width,prominence] = findpeaks(diff(y),'MinPeakProminence',0.5*max(pr)); 


        for j = 1:length(peaks)-1

            if (location(j)-offset)>=1
                
                ID1 = location(j)-offset;
                ID2 = location(j+1)-offset; 
                id1 = ID1+(ID1-1)*(DT-1);
                id2 = ID2+(ID2-1)*(DT-1); 
                temp = y_original ( id1 : id2 );

            else
                ID1 = location(j)-(location(j)-1); 
                ID2 = location(j+1)-(location(j)-1); 
                id1 = ID1+(ID1-1)*(DT-1); 
                id2 = ID2+(ID2-1)*(DT-1);
                temp = y_original( id1 : id2 );
            end

            regular_beat(z).data(j).beats = temp;

        end

        j = j+1;
        ID1 = location(j)-offset ; 
        id1 = ID1+(ID1-1)*(DT-1);
        temp2 = y_original(id1 : end);
        temp22 = y(ID1 : end);

        temp33 = [temp22; temp22(1)]; 
        temp3 = [temp2; temp2(1)]; 

        [p_temp,l_temp,w_temp,pr_temp] = findpeaks(temp33); 
        [pp_temp,ll_temp,ww_temp,prpr_temp] = findpeaks(temp33,'MinPeakProminence',0.2*max(pr_temp)); 

        if length(pp_temp)>1
            ID = (ll_temp(2)-offset);
            id = ID+(ID-1)*(DT-1);
            single_traces(z).data(j).beats = temp3(1:id);
        elseif (temp2(end)<temp2(1))
            regular_beat(z).data(j).beats = temp2;
        end

        regular_beat(z).original = y_original;
        z=z+1;

    end


    for i = 1:length(regular_beat)
        v = [];
        vv = [];
        v = [BBdistance(i).extra];
        vv = unique([v;v+1]);
        k = 1;

        for j = 1:length(regular_beat(i).data)
            if ismember(j,vv) == 0
                regular_beat(i).keep(k).beats = regular_beat(i).data(j).beats;
                k = k+1;
            end
        end

    end
    
    

    m = 1;
    for i = 1:length(regular_beat)
        if exist('regular_beat(i).keep') == 0
            ii(m)=i;
            m = m+1;
        end
    end

    if exist('ii')
        regular_beat(ii) = [];
    end


    k = 1;
    
    if length(regular_beat)>0
        
        for i = 1:length(regular_beat)
            l = 1;
            for jj = 1:length(regular_beat(i).keep)
                yy = regular_beat(i).keep(jj).beats;

                % fit
                [PP,XX] = findpeaks(yy,'SortStr','descend');
                half_vv = yy(XX(1):end);
                xx = 1:length(half_vv);
                CL = 1/PacingFrequency; %s
                n_time = 0.2*CL; % 20% of CL
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

                [fitresult2, gof2, xData2, yData2] = fit_tau(xx2, half_vv2');
                Sbeat(i).beat(jj).tau = ((1/(-fitresult2.b)).*(1/AcquisitionFrequency))*1000; %ms
                Sbeat(i).beat(jj).a = fitresult2.a;
                Sbeat(i).beat(jj).c = fitresult2.c;
                Sbeat(i).beat(jj).rsquare = gof2.rsquare; 

                pks = [];
                if length(half_vv3) > 3
                    [pks,locs,w,p] = findpeaks(half_vv3,'MinPeakHeight',bsl+0.05*bsl);                   
                end


                if (Sbeat(i).beat(jj).rsquare < 0.8) || (isempty(pks) == 0) % EADs or DADs
                    non_regular_peak(k).cell_n = i;
                    non_regular_peak(k).beat_n = jj;
                    non_regular_peak(k).beat_trace = yy;
                    irrPeaks (i,jj) = jj;
                    k = k+1;

                else
                    irrPeaks(i,jj) = 0;
                    regular_beat(i).keep2(l).beats = yy;
                    l = l+1;
                end

            end

        end


        for i = 1:length(regular_beat)  
            lengths = [];
            for j = 1:length(regular_beat(i).keep2)
                lengths(j,1) = length(regular_beat(i).keep2(j).beats);
            end
            min_length(i) = min(lengths);

            matrix = [];
            for j = 1:length(regular_beat(i).keep2)
                regular_beat(i).keep2(j).beats_cut  = regular_beat(i).keep2(j).beats(1:min_length(i));
                matrix = [matrix; regular_beat(i).keep2(j).beats_cut']; 
            end

            if size(matrix,1)>1
                regular_beat(i).mean(1,:) = mean(matrix);
            else
                regular_beat(i).mean(1,:) = matrix';
            end

        end


        for i = 1:length(regular_beat)
            [P,X] = findpeaks(regular_beat(i).mean,'SortStr','descend');
            half_v = regular_beat(i).mean(X(1):end);
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
            [fitresult, gof, xData, yData] = fit_tau(xx, half_vv);
            regular_beat(i).tau = ((1/(-fitresult.b)).*(1/AcquisitionFrequency))*1000; %ms
            regular_beat(i).a = fitresult.a;
            regular_beat(i).c = fitresult.c;
            regular_beat(i).rsquare = gof.rsquare;

        end

        b = 1;
        for file = 1:length(regular_beat)

            [CalciumPeak, CTD, CTD90, CTD50, CTD10, Time_C_peak, T_C_relax, AvgCalciumTrace, DoverD0,...
                yAVGCalciumTrace_baseline, Time_to_90a, Time_to_50a, Time_to_10a, Time_to_10Relax,...
                Time_to_50Relax, Time_to_90Relax, CalciumMagnitude,CalciumTDPre,CalciumTDPost] = parameters(regular_beat(file).mean,AcquisitionFrequency);

            try
            CD(file,1) = CTD;
            CD90(file,1) = CTD90;
            CD50(file,1) = CTD50;
            CD10(file,1) = CTD10;
            Ton (file,1) =  Time_C_peak;
            Toff (file,1) =  T_C_relax;   
            AVGCalcTrace(file).Avg = AvgCalciumTrace;   
            Fmax_F0(file,1) = DoverD0;
            Baseline(file,1) = yAVGCalciumTrace_baseline;
            T10on(file,1) = Time_to_90a;
            T50on(file,1) = Time_to_50a;
            T90on(file,1) =Time_to_10a;
            T10off(file,1) = Time_to_10Relax;
            T50off(file,1) = Time_to_50Relax;
            T90off(file,1) = Time_to_90Relax; 
            CaAmplitude(file,1) = CalciumMagnitude;
            tau(file,1) = regular_beat(file).tau;
            fit_a(file,1) = regular_beat(file).a;
            fit_c(file,1) = regular_beat(file).c;
            fit_rsquare(file,1) = regular_beat(file).rsquare;
            t0(file,1) = CalciumTDPre;
            tend(file,1) = CalciumTDPost;

            catch
                disp('>>>>>> error occurred in computing the parameters. Trace saved in "error_parameters"')
                error_parameters(b).data = regular_beat(file).mean';
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
                BBDistance(file,1) = NaN;
                NPeriod(file,1) = NaN;

            end

        end

        if_error_param = exist('error_parameters');

        K = 1:length(regular_beat); 
        K = K';
        for i = 1:length(regular_beat)
            removed_extra_beats(i,1) = length(regular_beat(i).data)-length(regular_beat(i).keep2);
            original_n_beats(i,1)= length(regular_beat(i).data);
        end

        [a,b]=find(irrPeaks>0);
        removed_IrrPeaks = zeros(length(regular_beat),1);
        for i = 1:length(a)
            removed_IrrPeaks(a(i),:) = length(find(irrPeaks(a(i),:)>0));
        end

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


        measurements_table = table( K,original_n_beats,removed_extra_beats,removed_IrrPeaks,...
          Baseline, Fmax_F0, CD, CD90, CD50, ...
          CD10, Ton, Toff, ...
          T10on, T50on, T90on, T10off, ...
          T50off, T90off, tau, fit_a, fit_c,fit_rsquare, t0, tend);


        if exist('error_parameters')
            measurements_table([error_parameters.id],:) = [];
        end

        cd(res_folder)
        writetable(measurements_table, fullfile('Calcium_measurements_extraBeats.xlsx') );
        if exist('ii')
            trace_entirely_removed = length(ii);
        end

        event_number = round((Number_of_cell_extra_beat)/(Number_of_cell_extra_beat+regular_cell_number) * 100);
        if exist('trace_entirely_removed')
            writetable(table(trace_entirely_removed,event_number,'VariableNames',{'traces_entirely_removed','event_number %'}), fullfile('Calcium_measurements_extraBeats.xlsx'),'Sheet',2 );
        end
        
    else
        disp('no regular beats in the extra beat trace(s)')
    end
    
    
else
    
    disp('no traces with extra beats found')
    
end

%%

cd(res_folder)
clearvars -except if_non_regular_decay PacingFrequency AcquisitionFrequency main_folder cutoff scripts_folder res_folder N_cell Number_of_cell_extra_beat Number_of_cell_non_regular_decay atLeast1Beat non_regular_beat regular_cell_number irrregular_cell_number

% if_non_regular_decay = exist ('Calcium_Traces_NonRegularBeats.xlsx', 'file');

if if_non_regular_decay == 2
    
    [ca_non_regular_decay,ca_non_regular_decay_txt,ca_non_regular_decay_raw] = xlsread('Calcium_Traces_NonRegularBeats.xlsx');
    
    j=1;
    for i = 1:2:size(ca_non_regular_decay,1)
        cell_n(j,1) = str2num(cell2mat(ca_non_regular_decay_txt(i+1,2)));
        cell_n(j,2) = str2num(cell2mat(ca_non_regular_decay_txt(i+2,2)));
        
        figure
        subplot(2,1,1)
        plot(ca_non_regular_decay(i,:)), grid on
        title(['cell ',num2str(cell_n(j,1))])
        subplot(2,1,2)
        plot(ca_non_regular_decay(i+1,find(ca_non_regular_decay(i+1,:)>0))), grid on
        title(['beat ',num2str(cell_n(j,2))])
        
        j=j+1;
    end
    cell_n(:,3) = [non_regular_beat.Nbeats];

    cell_nonReg = unique(cell_n(:,1));
    Number_of_cell_non_regular_decay = length(cell_nonReg)


    ca_multiBeat = ca_non_regular_decay(1:2:end,:);
    [un,un_id,c] = unique (cell_n(:,1));
    ca_multiBeat_un = ca_multiBeat(un_id,:);

    for i =1: length(un)
        idx = find(cell_n(:,1)==un(i));
        beats(i).n = cell_n(idx,2);
        beats(i).num = length(beats(i).n);
        cell_name = cell_n(un_id,1);
    end

    if atLeast1Beat ~= 0
    z = 1;
    BBdistance = [];

    CL = 1/PacingFrequency; %s
    offset_time = 0.1*CL;

    for i = 1:size(ca_multiBeat_un,1)
        regular_beat(z).name = ['cell',num2str(cell_name(i))];

        if AcquisitionFrequency>100
            y_original = ca_multiBeat_un(i,:);
            DT = round((1/100)*AcquisitionFrequency);
            y = movmean(y_original(1:DT:end),8);
            AF = 1/((1/AcquisitionFrequency)*DT);
        else
            y = ca_multiBeat_un(i,:);
            y_original = y;
            AF = AcquisitionFrequency;
            DT = 1;
        end
        offset = ceil(offset_time*AF);

        [p,l,w,pr] = findpeaks(diff(y));
        [peaks,location,width,prominence] = findpeaks(diff(y),'MinPeakProminence',0.5*max(pr)); %0.7
       

        [p2,l2,w2,pr2] = findpeaks(y);
        [peaks2,location2,width2,prominence2] = findpeaks(y,'MinPeakProminence',0.5*max(pr2)); %0.7


        for j = 1:length(peaks)-1

            if (location(j)-offset)>=1
                ID1 = location(j)-offset;
                ID2 = location(j+1)-offset;
                id1 = ID1+(ID1-1)*(DT-1);
                id2 = ID2+(ID2-1)*(DT-1);
                temp = y_original ( id1 : id2 );

            else
                ID1 = location(j)-(location(j)-1);
                ID2 = location(j+1)-(location(j)-1);
                id1 = ID1+(ID1-1)*(DT-1);
                id2 = ID2+(ID2-1)*(DT-1); 
                temp = y_original( id1 : id2 );
            end

            regular_beat(z).data(j).beats = temp;

        end

        j = j+1;
        ID1 = location(j)-offset ;
        id1 = ID1+(ID1-1)*(DT-1);
        temp2 = y_original(id1 : end);
        temp22 = y(ID1 : end);

        temp33 = [temp22'; temp22(1)];
        temp3 = [temp2'; temp2(1)];

        [p_temp,l_temp,w_temp,pr_temp] = findpeaks(temp33); 
        [pp_temp,ll_temp,ww_temp,prpr_temp] = findpeaks(temp33,'MinPeakProminence',0.2*max(pr_temp)); 

        regular_beat(z).original = y_original;
        z=z+1;
    
    end

    for i = 1:length(regular_beat)
        regular_beat(i).Noriginal = length(regular_beat(i).data);
        regular_beat(i).data([beats(i).n]) = [];
    end


    for i = 1:length(regular_beat)  
        lengths = [];
        
        if length(regular_beat(i).data) == 0
        
        else
            for j = 1:length(regular_beat(i).data)
                lengths(j,1) = length(regular_beat(i).data(j).beats);
            end
            min_length(i) = min(lengths);

            matrix = [];
            for j = 1:length(regular_beat(i).data)
                regular_beat(i).data(j).beats_cut  = regular_beat(i).data(j).beats(1:min_length(i));
                matrix = [matrix; regular_beat(i).data(j).beats_cut]; 
            end

            if size(matrix,1)>1
                regular_beat(i).mean(1,:) = mean(matrix);
            else
                regular_beat(i).mean(1,:) = matrix';
            end
            
        end
        
    end


    cd(main_folder)
    cd(scripts_folder)


    for i = 1:length(regular_beat)
        
        if exist('regular_beat(i).mean') %isempty(regular_beat(i).mean)
            regular_beat(i).a = 0;
            regular_beat(i).c = 0;
            regular_beat(i).rsquare = 0;
            
        else           
            try   
            [P,X] = findpeaks(regular_beat(i).mean,'SortStr','descend');
            half_v = regular_beat(i).mean(X(1):end);
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
            [fitresult, gof, xData, yData] = fit_tau(xx, half_vv);
            regular_beat(i).tau = ((1/(-fitresult.b)).*(1/AcquisitionFrequency))*1000; %ms
            regular_beat(i).a = fitresult.a;
            regular_beat(i).c = fitresult.c;
            regular_beat(i).rsquare = gof.rsquare;
            
            catch
                regular_beat(i).mean = [];
            end
        
        end
        
    end

    b = 1;
    for file = 1:length(regular_beat)
        
        if isempty(regular_beat(file).mean)
            
            CD(file,1) = 0;
            CD90(file,1) = 0;
            CD50(file,1) = 0;
            CD10(file,1) = 0;
            Ton (file,1) =  0;
            Toff (file,1) =  0;   
            AVGCalcTrace(file).Avg = 0;   
            Fmax_F0(file,1) = 0;
            Baseline(file,1) = 0;
            T10on(file,1) = 0;
            T50on(file,1) = 0;
            T90on(file,1) =0;
            T10off(file,1) = 0;
            T50off(file,1) = 0;
            T90off(file,1) = 0; 
            CaAmplitude(file,1) = 0;
            tau(file,1) = 0;
            fit_a(file,1) = 0;
            fit_c(file,1) = 0;
            fit_rsquare(file,1) = 0;
            t0(file,1) = 0;
            tend(file,1) = 0;
            BR(file,1) = 0;
            BBDistance(file,1) = 0;
            NPeriod(file,1) = 0;
            
        else

            [CalciumPeak, CTD, CTD90, CTD50, CTD10, Time_C_peak, T_C_relax, AvgCalciumTrace, DoverD0,...
                yAVGCalciumTrace_baseline, Time_to_90a, Time_to_50a, Time_to_10a, Time_to_10Relax,...
                Time_to_50Relax, Time_to_90Relax, CalciumMagnitude,CalciumTDPre,CalciumTDPost] = parameters(regular_beat(file).mean,AcquisitionFrequency);

            try
            CD(file,1) = CTD;
            CD90(file,1) = CTD90;
            CD50(file,1) = CTD50;
            CD10(file,1) = CTD10;
            Ton (file,1) =  Time_C_peak;
            Toff (file,1) =  T_C_relax;   
            AVGCalcTrace(file).Avg = AvgCalciumTrace;   

            Fmax_F0(file,1) = DoverD0;
            Baseline(file,1) = yAVGCalciumTrace_baseline;

            T10on(file,1) = Time_to_90a;
            T50on(file,1) = Time_to_50a;
            T90on(file,1) =Time_to_10a;
            T10off(file,1) = Time_to_10Relax;
            T50off(file,1) = Time_to_50Relax;
            T90off(file,1) = Time_to_90Relax; 

            CaAmplitude(file,1) = CalciumMagnitude;

            tau(file,1) = regular_beat(file).tau;
            fit_a(file,1) = regular_beat(file).a;
            fit_c(file,1) = regular_beat(file).c;
            fit_rsquare(file,1) = regular_beat(file).rsquare;

            t0(file,1) = CalciumTDPre;
            tend(file,1) = CalciumTDPost;

            catch
                disp('>>>>>> error occurred in computing the parameters."')
                error_parameters(b).data = regular_beat(file).mean';
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
                Time_to_90_Relax(file,1) = NaN; 
                T90off(file,1) = NaN;
                tau(file,1) = NaN;
                fit_a(file,1) = NaN;
                fit_c(file,1) = NaN;
                fit_rsquare(file,1) = NaN;
                t0(file,1) = NaN;
                tend(file,1) = NaN;
                BR(file,1) = NaN;
                BBDistance(file,1) = NaN;
                NPeriod(file,1) = NaN;

            end
        end

    end

    if_error_param = exist('error_parameters');

    K = 1:length(regular_beat); 
    K = K';

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

    removed_beats = [beats.num]';

    original_n_beats = [regular_beat.Noriginal]';
    measurements_table = table( cell_name,original_n_beats,removed_beats,...
      Baseline, Fmax_F0, CD, CD90, CD50, ...
      CD10, Ton, Toff, ...
      T10on, T50on, T90on, T10off, ...
      T50off, T90off, tau, fit_a, fit_c,fit_rsquare, t0, tend);


    if exist('error_parameters')
        measurements_table([error_parameters.id],:) = [];
    end

    cd(res_folder)
    writetable(measurements_table, fullfile('Calcium_measurements_IrregularBeats.xlsx') );
    event_number = round((Number_of_cell_non_regular_decay)/(Number_of_cell_non_regular_decay+regular_cell_number) * 100);
    writetable(table(event_number), fullfile('Calcium_measurements_IrregularBeats.xlsx'),'Sheet', 2 );

    end

else
    disp('no traces with irregular peaks found')


end


% if N_cell ~= 0
RegularCellNumber = regular_cell_number

EventNumber = round((Number_of_cell_extra_beat + Number_of_cell_non_regular_decay)/(Number_of_cell_extra_beat+Number_of_cell_non_regular_decay+regular_cell_number) * 100)

% end
disp('EventNumber is NonRegularCellNumber/(NonRegularCellNumber+RegularCellNumber). RegularCellNumber means all beats in the trace are regular.')
disp('RegularCellNumber may be different from the number of cells reported in Calcium_measurements.xlsx and Number_of_automatically_analysed_cell because Calcium_measurements.xlsx may contain measurements taken after excluding irregular beats.')
%% function to fit the calcium traces from peak to end using a bi-exp a*exp(b*x)+c*exp(d*x) with d=0;

function [fitresult, gof,xData, yData] = fit_tau(x, half_v)
             
[xData, yData] = prepareCurveData( x, half_v );

ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf 0];
opts.Upper = [Inf Inf Inf 0];

[fitresult, gof] = fit( xData, yData, ft, opts );

end