function [CalciumPeak, CTD, CTD90, CTD50, CTD10, Time_C_peak, T_C_relax, AvgCalciumTrace, DoverD0, yAVGCalciumTrace_baseline, ...
    Time_to_90a, Time_to_50a, Time_to_10a, Time_to_10Relax, Time_to_50Relax, Time_to_90Relax, CalciumMagnitude,CalciumTDPre,CalciumTDPost] = parameters(AverageTrace,AcquisitionFrequency)

    global PacingFrequency

    frametime = 1/AcquisitionFrequency; 


    apdLevel10 = 0.10; 
    apdLevel50 = 0.50; 
    apdLevel90 = 0.90;
    apdLevel100 = 1.00; 

    CL = 1/PacingFrequency; %s
    n_time = 0.2*CL;

    n = ceil(n_time*AcquisitionFrequency); 

    yAVGCalciumTrace_baseline = mean(AverageTrace(end-n:end)); 

    CalciumPeak =  max(AverageTrace);
    CalciumMagnitude = CalciumPeak - yAVGCalciumTrace_baseline;
    DoverD0 = CalciumPeak/yAVGCalciumTrace_baseline;

    CalciumTraceXval1 = (1:size(AverageTrace,2));  
    CalciumTraceXval2 = CalciumTraceXval1 .* frametime;
    CalciumTraceXvals = CalciumTraceXval2' ; 
    AvgCalciumTrace = AverageTrace';

    CalciumPeakLocation = find (AvgCalciumTrace == max(AvgCalciumTrace));
    CalPeak = CalciumTraceXvals(CalciumPeakLocation);
    Calcurvexpre = CalciumTraceXvals(1:CalciumPeakLocation);
    Calcurvexpost = CalciumTraceXvals(CalciumPeakLocation:end);
    Calcurveypre = (AvgCalciumTrace(1:CalciumPeakLocation))';
    Calcurveypost = (AvgCalciumTrace(CalciumPeakLocation:end))';

    CTDninetyfivevalue = yAVGCalciumTrace_baseline + (0.03 * (CalciumMagnitude));
    CTDvalue = CTDninetyfivevalue;
    CalciumMagnitude2 = CalciumPeak - CTDvalue;

    CTDninetyvalue = CTDvalue + (apdLevel10 * (CalciumMagnitude2));
    CTDfiftyvalue = CTDvalue + (apdLevel50 * (CalciumMagnitude2));
    CTDtenvalue = CTDvalue + (apdLevel90 * (CalciumMagnitude2));


    CTDlinePre = (ones (size(Calcurvexpre,1),1) .*CTDvalue);
    [CalciumTDPreX, CalciumTDPreY] = intersections(Calcurvexpre, Calcurveypre, Calcurvexpre, CTDlinePre, ' robust');
    CalciumTDPre = max(CalciumTDPreX);
    CTDlinePost = (ones (size(Calcurvexpost,1),1) .*CTDvalue);
    [CalciumTDPostX, CalciumTDPostY] = intersections(Calcurvexpost, Calcurveypost, Calcurvexpost, CTDlinePost, ' robust');
    CalciumTDPost = min (CalciumTDPostX);
    clear CTDvalue
    CTD = CalciumTDPost - CalciumTDPre;

    CTDNinetylinePre = (ones (size(Calcurvexpre,1),1) .*CTDninetyvalue);
    [CalciumTDPreX90, CalciumTDPreY90] = intersections(Calcurvexpre, Calcurveypre, Calcurvexpre, CTDNinetylinePre, ' robust');
    Calcium90TDPre = max(CalciumTDPreX90);
    CTDNinetylinePost = (ones (size(Calcurvexpost,1),1) .*CTDninetyvalue);
    [CalciumTDPostX90, CalciumTDPostY] = intersections(Calcurvexpost, Calcurveypost, Calcurvexpost, CTDNinetylinePost, ' robust');
    Calcium90TDPost = min (CalciumTDPostX90);
    Time_to_90a = Calcium90TDPre - CalciumTDPre;
    Time_to_90Relax = Calcium90TDPost - CalPeak;
    clear CTDninetyvalue
    CTD90 = Calcium90TDPost - Calcium90TDPre;

    CTDFiftylinePre = (ones (size(Calcurvexpre,1),1) .*CTDfiftyvalue);
    [CalciumTDPreX50, CalciumTDPreY] = intersections(Calcurvexpre, Calcurveypre, Calcurvexpre, CTDFiftylinePre, ' robust');
    Calcium50TDPre = max(CalciumTDPreX50);
    CTDFiftylinePost = (ones (size(Calcurvexpost,1),1) .*CTDfiftyvalue);
    [CalciumTDPostX50, CalciumTDPostY] = intersections(Calcurvexpost, Calcurveypost, Calcurvexpost, CTDFiftylinePost, ' robust');
    Calcium50TDPost = min (CalciumTDPostX50);
    Time_to_50a = Calcium50TDPre - CalciumTDPre; 
    Time_to_50Relax = Calcium50TDPost - CalPeak;
    clear CTDfiftyvalue
    CTD50 = Calcium50TDPost - Calcium50TDPre;

    CTDTenlinePre = (ones (size(Calcurvexpre,1),1) .*CTDtenvalue);
    [CalciumTDPreX10, CalciumTDPreY] = intersections(Calcurvexpre, Calcurveypre, Calcurvexpre, CTDTenlinePre, ' robust');
    Calcium10TDPre = max(CalciumTDPreX10);
    CTDTenlinePost = (ones (size(Calcurvexpost,1),1) .*CTDtenvalue);
    [CalciumTDPostX10, CalciumTDPostY] = intersections(Calcurvexpost, Calcurveypost, Calcurvexpost, CTDTenlinePost, ' robust');
    Calcium10TDPost = min (CalciumTDPostX10);
    Time_to_10a = Calcium10TDPre - CalciumTDPre ;
    Time_to_10Relax = Calcium10TDPost - CalPeak;
    clear CTDtenvalue 
    CTD10 = Calcium10TDPost - Calcium10TDPre;


    T_C_peak = (CalPeak - CalciumTDPre);

    Time_C_peak = T_C_peak;
    T_C_relax = ( CalciumTDPost - CalPeak ) ;

end

