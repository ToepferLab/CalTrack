function [video, Er, Sk, extra_beat] = BeatSegmentation2(ca_matrix,T0)

    %MultiCell

    global AcquisitionFrequency
    global cutoff
    global PacingFrequency

    video = struct.empty;
    Sk = struct.empty;
    extra_beat = struct.empty;
    Er = struct.empty;

    % k = 1;
    % e = 1;
    % m = 1;

    CL = 1/PacingFrequency; %s
    offset_time = 0.1*CL; 
    CL_frames = CL * AcquisitionFrequency;

    for i = 1:length(ca_matrix)

        t0 = T0(i,2);
        video(i).name = ca_matrix(i).name;

        for c = 1:size(ca_matrix(i).data,2)

            y = ca_matrix(i).data(:,c);


        [~,~,~,pr] = findpeaks(movmean(y,5));
        [peaks,location,~,~] = findpeaks(movmean(y,5),'MinPeakProminence',0.5*max(pr)); 

        for j = 0:length(peaks)-2
            video(i).single_traces(c).data(j+1).beats = y(t0+j*CL_frames : t0+(j+1)*CL_frames);
        end



        peaks_time = location/AcquisitionFrequency;
        peaks_period = (diff(peaks_time));
        video(i).single_traces(c).Nbeats = length(peaks);
        video(i).single_traces(c).time = length(y)/AcquisitionFrequency; %time in sec of the whole trace
        video(i).single_traces(c).BR = video(i).single_traces(c).Nbeats/video(i).single_traces(c).time; % beats per sec
        video(i).single_traces(c).BBdistance = mean(peaks_period);
        video(i).single_traces(c).Nperiod = length(peaks_period);
        video(i).single_traces(c).original = y;


        end
    end
end

