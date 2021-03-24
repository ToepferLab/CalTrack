function [single_traces] = PhotoBleachCorrection(single_traces)

    global PacingFrequency
    global AcquisitionFrequency

    CL = 1/PacingFrequency; %s
    n_time = 0.2*CL;

    n = ceil(n_time*AcquisitionFrequency); 
    for i = 1:length(single_traces)

        cll = struct2cell(single_traces(i).data);
        sqcll = squeeze(cll);
        selected_trace = cat(1,sqcll{:});

        for K = 1:size(cll,3)
            single_trace = cll{:,:,K};
            peak = max(single_trace);
            Baseline = mean(single_trace(end-n:end));
            Magn = peak - Baseline;
            val95 = Baseline + (0.05 * Magn); 
            BaselineValuePts = (single_trace<val95);
            BaselineVals = (single_trace .* BaselineValuePts)';
            BaselineValues(K).data = BaselineVals;
        end

        baseline_trace = horzcat(BaselineValues.data);
        vs = baseline_trace>0;
        y = 1:size (baseline_trace,2);
        ycords  = (vs .*y)';
        baseline_trace = baseline_trace';
        ys = baseline_trace(baseline_trace(:,1) > .0,:);
        xs = ycords(baseline_trace(:,1) > .0,:);

        t = 1:length(selected_trace);
        [pt,st,mut] = polyfit(xs,ys,2);
        trend = polyval(pt,t,[],mut);

        adj = selected_trace - trend';  
        Corrected_Trace =  adj + (selected_trace(1)- adj(1));  

        single_traces(i).BaselineCorrectedTraces = Corrected_Trace;
        
    end
    
end

