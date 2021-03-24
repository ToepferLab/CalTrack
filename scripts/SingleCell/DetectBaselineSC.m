function [baseline] = DetectBaselineSC(single_traces,AcquisitionFrequency,PacingFrequency)

    CL = 1/PacingFrequency; %s
    n_time = 0.2*CL;

    n = ceil(n_time*AcquisitionFrequency); 

    for i = 1:length(single_traces)

        cll = {single_traces(i).data.beats}.';
%         selected_trace = vertcat(1,cll{:});

        for K = 1:size(cll,1)
            single_trace = cll{K,:};
            if size (single_trace,1)> n
                peak = max(single_trace);
                Baseline = nanmean(single_trace(end-n:end));
                Magn = peak - Baseline;
                val95 = Baseline + (0.03 * Magn); 
                BaselineValuePts = (single_trace<val95);
                BaselineVals = (single_trace .* BaselineValuePts)';
                BaselineValues(K).data = BaselineVals;
                TransientValuePts = (single_trace>val95);
                TransientVals = (single_trace .* TransientValuePts);
                TransientValues(K).data = TransientVals;
            end
        end

        baseline_trace = horzcat(BaselineValues.data);
%         vs = baseline_trace>0;
%         y = 1:size (baseline_trace,2);
%         ycords  = (vs .*y)';
        baseline_trace = baseline_trace';
        ys = baseline_trace(baseline_trace(:,1) > .0,:);

        Transient_trace = vertcat(TransientValues.data);
%         tvs = Transient_trace>0;
%         ty = 1:size (Transient_trace,2);
%         tycords  = (tvs .*ty)';
%         tys = Transient_trace(Transient_trace(:,1) > .0,:);
        ts = Transient_trace(Transient_trace(:,1) > .0,:);

        baseline(i).BaselineTraces = ys ;
        baseline(i).TransientTraces = ts ;
        baseline(i).name = single_traces(i).name;
        clearvars  BaselineValues
        
    end

end
    


