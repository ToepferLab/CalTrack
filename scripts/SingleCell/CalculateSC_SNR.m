function [SNR] = CalculateSC_SNR(ca_matrix, AcquisitionFrequency, ~, PacingFrequency)

%     CaTraces = ca_matrix.data;
%     NAMEs = (char(ca_matrix.name));

    [baselineN] = DetectBaselineSC(ca_matrix,AcquisitionFrequency,PacingFrequency);

    for h = 1:length (baselineN)
    %     MeanBaselineValueN = mean(nonzeros(baselineN(h).BaselineTraces));
        VarBaselineN = std(nonzeros(baselineN(h).BaselineTraces));
        MeanTransientValueN = mean(nonzeros(baselineN(h).TransientTraces));
        SNRN = MeanTransientValueN/VarBaselineN;

        SNR(h).name = baselineN.name;
        SNR(h).AvgNoise = VarBaselineN;
        SNR(h).AvgSignal = MeanTransientValueN;
        SNR(h).SignalToNoiseRatio = SNRN;
    end
end
