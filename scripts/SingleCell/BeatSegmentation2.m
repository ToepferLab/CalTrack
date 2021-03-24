function [single_traces, skipped, extra_beat, errors] = BeatSegmentation2(NAME, Y, T0)

    % BeatSegmentation_SingleCell.m

    global AcquisitionFrequency
%     global cutoff
    global PacingFrequency


    Names = NAME;

    single_traces = struct.empty;
    skipped = struct.empty;
    extra_beat = struct.empty;
    errors = struct.empty;


    CL = 1/PacingFrequency; %s
%     offset_time = 0.1*CL;
    CL_frames = CL * AcquisitionFrequency;


    for i = 1:length(Y)

        y = Y{i};
        t0 = T0(i,2);

        [~,~,~,pr] = findpeaks(y);
        [peaks,location,~,~] = findpeaks(y,'MinPeakProminence',0.5*max(pr));


        for j = 0:length(peaks)-2
            single_traces(i).data(j+1).beats = y(t0+j*CL_frames : t0+(j+1)*CL_frames);
        end
        j = j+1;
        single_traces(i).data(j+1).beats = y(t0+j*CL_frames : end);

        single_traces(i).name = Names{i,1}; 
        peaks_time = location/AcquisitionFrequency;
        peaks_period = (diff(peaks_time));
        single_traces(i).Nbeats = length(peaks);
        single_traces(i).time = length(y)/AcquisitionFrequency;
        single_traces(i).BR = single_traces(i).Nbeats/single_traces(i).time;
        single_traces(i).BBdistance = mean(peaks_period);
        single_traces(i).Nperiod = length(peaks_period);
        single_traces(i).original = y;

    end


end

