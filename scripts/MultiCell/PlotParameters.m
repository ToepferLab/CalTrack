
AF = AcquisitionFrequency;

s_names = sheetnames('Calcium_Traces.xlsx');
for i = 1:length(s_names)
    ca_measurements(i).data = readtable('Calcium_measurements_forPlot.xlsx','Sheet',i);
    ca_traces(i).data = readtable('Calcium_Traces.xlsx','Sheet',s_names(i));
end


figure(100),clf
for i = 1:length(ca_traces)

    subplot(3,ceil(length(ca_traces)/3),i)
    video = table2array(ca_traces(i).data);
    for j = 1:size(video,2)
        cell = video(:,j);
        plot((cell-min(cell))/(max(cell)-min(cell))), hold on, grid on
        title(strcat('video',num2str(i)))
    end
end



for i = 1:length(ca_traces)
    figure('Name',strcat('video',num2str(i),' (fit)')), clf
    video = table2array(ca_traces(i).data);
    cell = [];
    for j = 1:size(video,2)
        subplot(3,ceil(size(video,2)/3),j)
        cell = video(:,j);
        x = [1:length(cell)]*(1/AF)*1000;
        plot(x,cell), hold on, grid on
        title(strcat('cell',num2str(j)))
        
        param = ca_measurements(i).data{j,3:end-2};
        

        [M,id] = max(cell);
        xx = x(id:end);
        a = param(:,16);
        c = param(:,17);
        b = 1./(param(:,15)*(-1));
        yy = a*exp(b*x)+c;
        hold on
        plot(x+x(id),yy)
        
        v = a*(1/exp(1))+c; 
        plot(x,ones(1,length(x))*v)
        plot([param(7)+param(end-1)+param(15) param(7)+param(end-1)+param(15)],[min(cell) max(cell)])
        xlabel('Time (ms)')        
    end
%     legend('calcium trace','fit (y=a*exp(b*x)+c)','y0/e+c','\tau')

end


for i = 1:length(ca_traces)
    figure('Name',strcat('video',num2str(i),' (baseline, ton, toff, peak)')), clf
    video = table2array(ca_traces(i).data);
    cell = [];
    for j = 1:size(video,2)
        subplot(3,ceil(size(video,2)/3),j)
        cell = video(:,j);
        x = [1:length(cell)]*(1/AF)*1000;
        plot(x,cell), hold on, grid on
        title(strcat('cell',num2str(j)))
        xlabel('Time (ms)') 
        param = ca_measurements(i).data{j,3:end-2};
        
        plot(x,ones(1,length(x))*param(1))
        plot([param(end-1) param(end-1)],[min(cell) max(cell)])
        plot([param(end) param(end)],[min(cell) max(cell)])
        plot([param(7)+param(end-1) param(7)+param(end-1)],[min(cell) max(cell)])
               
    end
%         legend('calcium trace','baseline','t_o_n','t_o_f_f','peak (from t_o_n)')

end


for i = 1:length(ca_traces)
    figure('Name',strcat('video',num2str(i),' (contr and relax times)')), clf
    video = table2array(ca_traces(i).data);
    cell = [];
    for j = 1:size(video,2)
        subplot(3,ceil(size(video,2)/3),j)
        cell = video(:,j);
        x = [1:length(cell)]*(1/AF)*1000;
        plot(x,cell), hold on, grid on
        title(strcat('cell',num2str(j)))
        xlabel('Time (ms)') 
        param = ca_measurements(i).data{j,3:end-2};
        
        plot([param(9)+param(end-1) param(9)+param(end-1)],[min(cell) max(cell)])
        plot([param(10)+param(end-1) param(10)+param(end-1)],[min(cell) max(cell)])
        plot([param(11)+param(end-1) param(11)+param(end-1)],[min(cell) max(cell)])

        plot([param(12)+param(7)+param(end-1) param(12)+param(7)+param(end-1)],[min(cell) max(cell)])
        plot([param(13)+param(7)+param(end-1) param(13)+param(7)+param(end-1)],[min(cell) max(cell)])
        plot([param(14)+param(7)+param(end-1) param(14)+param(7)+param(end-1)],[min(cell) max(cell)])
                 
    end
%         legend('calcium trace','t10 contr (from t_o_n)','t50 contr (from t_o_n)','t90 contr (from t_o_n)','t10 relax (from peak)','t50 relax (from peak)','t90 relax (from peak)')

end

