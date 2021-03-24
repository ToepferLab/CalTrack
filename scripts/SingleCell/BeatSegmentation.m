function [single_traces, skipped, extra_beat, errors] = BeatSegmentation(NAME, Y, param_single_beat)

    % BeatSegmentation_SingleCell.m

    global AcquisitionFrequency
    global cutoff
    global PacingFrequency


    Names = NAME;


    single_traces = struct.empty;
    skipped = struct.empty;
    extra_beat = struct.empty;
    errors = struct.empty;

    k = 1;
    m = 1;
    z = 1;
    e = 1;

    CL = 1/PacingFrequency;
    offset_time = 0.1*CL;

    for i = 1:length(Y)

        if AcquisitionFrequency>100
            y_original = Y{i};
            DT = round((1/100)*AcquisitionFrequency);
            y1 = movmean(Y{i},24);
            y = y1(1:DT:end);
            AF = 1/((1/AcquisitionFrequency)*DT);
        else
            y = Y{i};
            y_original = y;
            AF = AcquisitionFrequency;
            DT = 1;
        end
        
        offset = ceil(offset_time*AF);

        skip = 0;

        y_diff = movmean(diff(y),8);

        [~,~,~,pr] = findpeaks(y_diff);
        [peaks,location,~,~] = findpeaks(y_diff,'MinPeakProminence',0.5*max(pr));

        [~,~,~,pr2] = findpeaks(y);
        [peaks2,location2,~,~] = findpeaks(y,'MinPeakProminence',0.5*max(pr2));

        peaks_time = location2/AF;
        peaks_period = (diff(peaks_time));
        is_extra_beat = find(peaks_period<cutoff);  

        if (length(peaks)-length(peaks2))>1
            
            disp('>>>>>> noisy signal detected and excluded.')        
            skip = 1;
            skipped(k).data = y_original;
            skipped(k).name = Names{i,1};
            skipped(k).id = i;
            k=k+1;

        elseif isempty(is_extra_beat) == 0

            disp('>>>>>> extra beat(s) detected and excluded.')
            extra_beat(m).data = y_original;
            extra_beat(m).name = Names{i,1};
            extra_beat(m).id = i;
            extra_beat(m).BBdistance = mean(peaks_period);
            extra_beat(m).Nperiod = length(peaks_period); 
            m=m+1;
            skip = 1;
            
        end

        if skip == 0
            
            try

            single_traces(z).name = Names{i,1};  

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

                single_traces(z).data(j).beats = temp;

            end

            j = j+1;
            ID1 = location(j)-offset ; 
            id1 = ID1+(ID1-1)*(DT-1); 
            temp2 = y_original(id1 : end);
            temp22 = y(ID1 : end);

            temp33 = [temp22; temp22(1)]; 
            temp3 = [temp2; temp2(1)]; 

            [~,~,~,pr_temp] = findpeaks(temp33);
            [pp_temp,ll_temp,~,~] = findpeaks(temp33,'MinPeakProminence',0.2*max(pr_temp));

            if (AcquisitionFrequency>100 || param_single_beat == 'y')
                % 
            else

                if length(pp_temp)>1 
                    ID = (ll_temp(2)-offset);
                    id = ID+(ID-1)*(DT-1);
                    single_traces(z).data(j).beats = temp3(1:id);
                elseif (temp2(end)<temp2(1))
                    single_traces(z).data(j).beats = temp2;
                end
                
            end

            single_traces(z).Nbeats = length(single_traces(z).data);

            len=0;
            for j=1:length(single_traces(z).data)
                len = len +length(single_traces(z).data(j).beats);
            end
            single_traces(z).time = len/AcquisitionFrequency; 
            single_traces(z).BR = single_traces(z).Nbeats/single_traces(z).time; 
            single_traces(z).BBdistance = mean(peaks_period);
            single_traces(z).Nperiod = length(peaks_period);
            single_traces(z).original = y_original;
            z=z+1;

            catch
                disp('>>>>>> error occurred in beat segmentation.')
                errors(e).data = y_original;
                errors(e).name = Names{i,1}; 
                e = e+1;
            end
            
        end

    end

end
    



