function [video, Er, Sk, extra_beat] = BeatSegmentation(ca_matrix)

    global AcquisitionFrequency
    global cutoff
    global PacingFrequency

    video = struct.empty;
    Sk = struct.empty;
    extra_beat = struct.empty;
    Er = struct.empty;

    k = 1;
    e = 1;
    m = 1;
    for i = 1:length(ca_matrix)  

        z = 1;

        for c = 1:size(ca_matrix(i).data,2)  


            if AcquisitionFrequency>100
                y = ca_matrix(i).data(:,c);
                DT = round((1/100)*AcquisitionFrequency);  
                y_movmean = movmean(y(1:DT:end),5);
                AF = 1/((1/AcquisitionFrequency)*DT);  
            else
                y = ca_matrix(i).data(:,c);
                y_movmean = movmean(y,5);
                AF = AcquisitionFrequency;
                DT = 1;
            end        


            X = (length(y_movmean))/AF;
            Nbeats_expected = round(X*PacingFrequency);


            y_movmean_diff = diff(y_movmean);
            yy = movmean(y_movmean_diff,10);

            skip = 0;


            [~,~,~,pr] = findpeaks(yy); 
            [peaks,location,~,~] = findpeaks(yy,'MinPeakProminence',0.5*max(pr)); 


            [~,~,~,pr2] = findpeaks(y_movmean);
            [peaks2,location2,~,~] = findpeaks(y_movmean,'MinPeakProminence',0.5*max(pr2)); 

            peaks_time = location2/AF;
            peaks_period = (diff(peaks_time));
            is_extra_beat = find(peaks_period<cutoff);


            if (length(peaks)-length(peaks2))>1 % 

                disp('>>>>>> noisy signal detected and excluded. Saved in "skipped"')

                skip = 1;
                Sk(k).data = y;
                Sk(k).Video_id = i;
                Sk(k).Video_name = ca_matrix(i).name;
                Sk(k).Cell_id = c;
                k=k+1;

            elseif isempty(is_extra_beat) == 0

                disp('>>>>>> extra beat(s) detected and excluded. Saved in "extra_beat"')
                extra_beat(m).data = y;
                extra_beat(m).Video_id = i;
                extra_beat(m).Video_name = ca_matrix(i).name;
                extra_beat(m).Cell_id = c;
                extra_beat(m).BBdistance = mean(peaks_period);
                extra_beat(m).Nperiod = length(peaks_period); 
                m=m+1;
                skip = 1;

            end

            CL = 1/PacingFrequency; 
            offset_time = 0.1*CL; 

            offset = ceil(offset_time*AcquisitionFrequency);

            if skip == 0
                try

                for j = 1:length(peaks)-1 

                    if (location(j)-offset)>=1
                        ID1 = location(j)-offset; 
                        ID2 = location(j+1)-offset; 
                        id1 = ID1+(ID1-1)*(DT-1); 
                        id2 = ID2+(ID2-1)*(DT-1);
                        temp = y ( id1 : id2 );              

                    else
                        ID1 = location(j)-(location(j)-1); 
                        ID2 = location(j+1)-(location(j)-1);
                        id1 = ID1+(ID1-1)*(DT-1); 
                        id2 = ID2+(ID2-1)*(DT-1);
                        temp = y( id1 : id2 ); 
                    end

                    video(i).single_traces(z).data(j).beats = temp;

                end

                video(i).single_traces(z).Nbeats = length(video(i).single_traces(z).data);

                len=0;
                for j=1:length(video(i).single_traces(z).data)
                    len = len +length(video(i).single_traces(z).data(j).beats);
                end
                video(i).single_traces(z).time = len/AcquisitionFrequency; %time in sec of the whole trace
                video(i).single_traces(z).BR = video(i).single_traces(z).Nbeats/(video(i).single_traces(z).time); % beats per sec
                video(i).single_traces(z).BBdistance = mean(peaks_period);
                video(i).single_traces(z).Nperiod = length(peaks_period); 
                video(i).single_traces(z).original = y; 
                video(i).name = ca_matrix(i).name;

                z=z+1;

                catch
                    disp('>>>>>> error occurred. Traced saved in "errors"')
                    Er(e).data = y;
                    Er(e).Video_id = i;
                    Er(e).Video_name = ca_matrix(i).name;
                    Er(e).Cell_id = c;
                    e = e+1;
                end
            end
        end
    end
end

