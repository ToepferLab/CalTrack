function [SNR] = CalculateSNR(single_tracesN, AcquisitionFrequency, ~, PacingFrequency)

    [baselineN] = DetectBaseline(single_tracesN,AcquisitionFrequency,PacingFrequency);

    for h = 1:length (baselineN)
        MeanBaselineValueN = mean(nonzeros(baselineN(h).BaselineTraces));
        VarBaselineN = std(nonzeros(baselineN(h).BaselineTraces));
        MeanTransientValueN = mean(nonzeros(baselineN(h).TransientTraces));
        SNRN = MeanTransientValueN/VarBaselineN;

        SNR(h).AvgNoise = VarBaselineN;
        SNR(h).AvgSignal = MeanTransientValueN;
        SNR(h).SignalToNoiseRatio = SNRN;
    end
end
