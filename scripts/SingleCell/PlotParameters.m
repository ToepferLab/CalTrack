
AF = AcquisitionFrequency;

[ca_measurements,ca_measurements_txt,ca_measurements_raw] = xlsread('Calcium_measurements.xlsx');
[ca_traces,ca_traces_txt,ca_traces_raw] = xlsread('Calcium_Traces.xlsx');

BBdistance_extra_beat = ca_measurements(:,end-3);
Nperiod_extra_beat = ca_measurements(:,end-2);
SNR = ca_measurements(:,end);

BBdistance = ca_measurements(:,end-5);
Nperiod = ca_measurements(:,end-4);

ca_measurements = ca_measurements(:,1:end-6);


if_extra_beat = exist ('Calcium_Traces_extra_beat.xlsx', 'file');

if if_extra_beat == 2
    
    figure('Name','Calcium Traces - extra beat(s)')
    s_names_eb = sheetnames('Calcium_Traces_extra_beat.xlsx');
    for i=1:length(s_names_eb)
        
        subplot(4,ceil(length(s_names_eb)/4),i)
        [ca_traces_extra_beat,ca_traces_txt_extra_beat,ca_traces_raw_extra_beat] = xlsread('Calcium_Traces_extra_beat.xlsx',i);
        x = [1:size(ca_traces_extra_beat,1)]*(1/AF)*1000;
        plot(x,ca_traces_extra_beat)
        xlabel('Time (ms)')
        grid on
        title(strcat('Process',ca_traces_txt_extra_beat{2, 1}(9:end)))

    end
    
    ca_measurements = ca_measurements(1:end-length(s_names_eb),:);
    Number_of_cell_extra_beat = length(s_names_eb)

end

Number_of_automatically_analysed_cell = size(ca_measurements,1)


if_noisy = exist ('Calcium_Traces_noisy.xlsx', 'file');

if if_noisy == 2
    
    figure('Name','Calcium Traces - noisy')
    s_names_n = sheetnames('Calcium_Traces_noisy.xlsx');
    for i=1:length(s_names_n)
        
        subplot(4,ceil(length(s_names_n)/4),i)
        [ca_traces_noisy,ca_traces_txt_noisy,ca_traces_raw_noisy] = xlsread('Calcium_Traces_noisy.xlsx',i);
        x = [1:size(ca_traces_noisy,1)]*(1/AF)*1000;
        plot(x,ca_traces_noisy)
        xlabel('Time (ms)')
        grid on
        title(strcat('Process',ca_traces_txt_noisy{2, 1}(9:end)))

    end
    Number_of_cell_noisy = length(s_names_n)

end


if_error_param = exist ('Calcium_Traces_error_param.xlsx', 'file');

if if_error_param == 2
    
    figure('Name','Calcium Traces - error param')
    s_names_erb = sheetnames('Calcium_Traces_error_param.xlsx');
    
    for i=1:length(s_names_erb)
        [ca_traces_error_param,ca_traces_txt_error_param,ca_traces_raw_error_param] = xlsread('Calcium_Traces_error_param.xlsx',i);
        subplot(4,ceil(length(s_names_erb)/4),i)
        x = [1:size(ca_traces_error_param,1)]*(1/AF)*1000;
        plot(x,ca_traces_error_param)
        xlabel('Time (ms)')
        grid on
        title(strcat('Process',ca_traces_txt_error_param{2, 1}(9:end)))
    end
    Number_of_cell_error_param = length(s_names_erb)
end


if_error = exist ('Calcium_Traces_errors.xlsx', 'file');

if if_error == 2
    
    figure('Name','Calcium Traces - non-specific error')
    s_names_er = sheetnames('Calcium_Traces_errors.xlsx');
    
    for i=1:length(s_names_er)
        [ca_traces_error,ca_traces_txt_error,ca_traces_raw_error] = xlsread('Calcium_Traces_errors.xlsx',i);
        subplot(4,ceil(length(s_names_er)/4),i)
        x = [1:size(ca_traces_error,1)]*(1/AF)*1000;
        plot(x,ca_traces_error)
        xlabel('Time (ms)')
        grid on
        title(strcat('Process',ca_traces_txt_error{2, 1}(9:end)))
    end
    Number_of_cell_error = length(s_names_er)
end


%%
n = size(ca_measurements,2); 
param_names = {'Baseline','DoverD0','Mean Calcium CD','Mean Calcium CD90',...
    'Mean Calcium CD50','Mean Calcium CD10Max','Time CalciumPeak','Time CalciumRelax',...
    'Time to 10Contr','Time to 50Contr','Time to 90Contr','Time to 10 Relax',...
    'Time to 50 Relax','Time to 90 Relax','tau','fit a','fit c','fit rsquare','BR','t0','tend','BBtime_mean','Nperiod'};


figure

subplot(1,2,1)
x = [1:size(ca_traces,1)]*(1/AF)*1000;
plot(x,ca_traces)
xlabel('Time (ms)')
ylabel('Calcium Traces')
grid on
title({'each trace is';'a different cell';'(mean over';'multiple beats)'})

subplot(1,2,2)
for i = 1:size(ca_traces,2)
    plot(x,(ca_traces(:,i)-min(ca_traces(:,i)))/(max(ca_traces(:,i))-min(ca_traces(:,i)))), hold on
end
xlabel('Time (ms)')
ylabel('Normalised Calcium Traces')
grid on


for i=1:size(ca_traces,2)
    y = ca_traces(:,i);
    figure('Name',(strcat('Process',ca_traces_txt{2, i}(9:end))))
    subplot(1,3,1)
    plot(x,y), hold on, grid on

    % fit
    [M,id] = max(y);
    xx = x(id:end);
    a = ca_measurements(:,16);
    c = ca_measurements(:,17);
    b = 1./(ca_measurements(:,15)*(-1));

    half_v = y(id:end);
    Xx = 1:length(half_v);
    CL = 1/PacingFrequency; %s
    n_time = 0.2*CL; % 20% of CL
    n = ceil(n_time*AcquisitionFrequency);            
    bsl = mean(half_v(end-n:end));
    bsln = ((ones (length(half_v),1) .*bsl))';
    [xint, yint] = intersections(Xx, half_v, Xx, bsln, ' robust');
    maxx = floor(min(xint))+ 3;
    if maxx > length (half_v)
        maxx = length (half_v);
    end
    half_vv = half_v(1:maxx);
    xx = 1:length(half_vv);
    
    yy = a(i)*exp(b(i)*x(1:end))+c(i);
    hold on
    plot(x(1:end)+x(id),yy)

    v = a(i)*(1/exp(1))+c(i);
    plot(x,ones(1,length(x))*v)
    plot([ca_measurements(i,7)+ca_measurements(i,end-1)+ca_measurements(i,15) ca_measurements(i,7)+ca_measurements(i,end-1)+ca_measurements(i,15)],[min(y) max(y)])
%     legend('calcium trace','fit (y=a*exp(b*x)+c)','a/e+c','\tau')
    xlabel('Time (ms)')

    subplot(1,3,2)
    plot(x,y), hold on, grid on
    plot(x,ones(1,length(x))*ca_measurements(i,1))
    plot([ca_measurements(i,end-1) ca_measurements(i,end-1)],[min(y) max(y)])
    plot([ca_measurements(i,end) ca_measurements(i,end)],[min(y) max(y)])
    plot([ca_measurements(i,7)+ca_measurements(i,end-1) ca_measurements(i,7)+ca_measurements(i,end-1)],[min(y) max(y)])
%     legend('calcium trace','baseline','t_0','t_e_n_d','peak (from t_0)')
    xlabel('Time (ms)')

    subplot(1,3,3)
    plot(x,y), hold on, grid on
    plot([ca_measurements(i,9)+ca_measurements(i,end-1) ca_measurements(i,9)+ca_measurements(i,end-1)],[min(y) max(y)])
    plot([ca_measurements(i,10)+ca_measurements(i,end-1) ca_measurements(i,10)+ca_measurements(i,end-1)],[min(y) max(y)])
    plot([ca_measurements(i,11)+ca_measurements(i,end-1) ca_measurements(i,11)+ca_measurements(i,end-1)],[min(y) max(y)])

    plot([ca_measurements(i,12)+ca_measurements(i,7)+ca_measurements(i,end-1) ca_measurements(i,12)+ca_measurements(i,7)+ca_measurements(i,end-1)],[min(y) max(y)])
    plot([ca_measurements(i,13)+ca_measurements(i,7)+ca_measurements(i,end-1) ca_measurements(i,13)+ca_measurements(i,7)+ca_measurements(i,end-1)],[min(y) max(y)])
    plot([ca_measurements(i,14)+ca_measurements(i,7)+ca_measurements(i,end-1) ca_measurements(i,14)+ca_measurements(i,7)+ca_measurements(i,end-1)],[min(y) max(y)])
%     legend('calcium trace','t10 contr (from t_0)','t50 contr (from t_0)','t90 contr (from t_0)','t10 relax (from peak)','t50 relax (from peak)','t90 relax (from peak)')
    xlabel('Time (ms)')

end
