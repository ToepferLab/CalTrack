function [CalciumPeak, CTD, CTD90, CTD50, CTD10, Time_C_peak, T_C_relax, AvgCalciumTrace, DoverD0, yAVGCalciumTrace_baseline, ...
    Time_to_90a, Time_to_50a, Time_to_10a, Time_to_10Relax, Time_to_50Relax, Time_to_90Relax, CalciumMagnitude,CalciumTDPre,CalciumTDPost] = parameters(AverageTrace,AcquisitionFrequency)

global PacingFrequency


 frametime = 1/AcquisitionFrequency; 
 
% if AcquisitionFrequency>100
%  AverageTrace = movmean(AverageTrace,20);
% end

 apdLevel10 = 0.10; 
 apdLevel50 = 0.50; 
 apdLevel90 = 0.90;
 apdLevel100 = 1.00; 

 CL = 1/PacingFrequency; %s
 n_time = 0.2*CL; % 20% of CL

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
 CalPeak = CalciumTraceXvals(CalciumPeakLocation); % time to peak (sec)
 Calcurvexpre = CalciumTraceXvals(1:CalciumPeakLocation);
 Calcurvexpost = CalciumTraceXvals(CalciumPeakLocation:end);
 Calcurveypre = (AvgCalciumTrace(1:CalciumPeakLocation))';
 Calcurveypost = (AvgCalciumTrace(CalciumPeakLocation:end))';

%  CTDvalue =  yAVGCalciumTrace_baseline ;
 
% new magnitude based of the 95% of the signal
 CTDninetyfivevalue = yAVGCalciumTrace_baseline + (0.03 * (CalciumMagnitude)); % 0.03
 CTDvalue = CTDninetyfivevalue;
 CalciumMagnitude2 = CalciumPeak - CTDvalue;
 
 CTDninetyvalue = CTDvalue + (apdLevel10 * (CalciumMagnitude2));
 CTDfiftyvalue = CTDvalue + (apdLevel50 * (CalciumMagnitude2));
 CTDtenvalue = CTDvalue + (apdLevel90 * (CalciumMagnitude2));


 
%%%%% CTD
 CTDlinePre = (ones (size(Calcurvexpre,1),1) .*CTDvalue);
 [CalciumTDPreX, CalciumTDPreY] = intersections(Calcurvexpre, Calcurveypre, Calcurvexpre, CTDlinePre, ' robust');
 CalciumTDPre = max(CalciumTDPreX);
 CTDlinePost = (ones (size(Calcurvexpost,1),1) .*CTDvalue);
 [CalciumTDPostX, CalciumTDPostY] = intersections(Calcurvexpost, Calcurveypost, Calcurvexpost, CTDlinePost, ' robust');
%  if length(CalciumTDPostX)>3
%     CalciumTDPost = mean (CalciumTDPostX(1:3)); %min
%  else
     CalciumTDPost = min (CalciumTDPostX);
%  end
 
 clear CTDvalue
 CTD = CalciumTDPost - CalciumTDPre;
 

%%%%% CTD 90
 CTDNinetylinePre = (ones (size(Calcurvexpre,1),1) .*CTDninetyvalue);
 [CalciumTDPreX90, CalciumTDPreY90] = intersections(Calcurvexpre, Calcurveypre, Calcurvexpre, CTDNinetylinePre, ' robust');
 Calcium90TDPre = min(CalciumTDPreX90); %min
 if Calcium90TDPre < CalciumTDPre
     Calcium90TDPre = CalciumTDPre;
 end
 CTDNinetylinePost = (ones (size(Calcurvexpost,1),1) .*CTDninetyvalue);
 [CalciumTDPostX90, CalciumTDPostY90] = intersections(Calcurvexpost, Calcurveypost, Calcurvexpost, CTDNinetylinePost, ' robust');
 Calcium90TDPost = max (CalciumTDPostX90);
 if Calcium90TDPost>CalciumTDPost
     Calcium90TDPost = CalciumTDPost;
 end
 Time_to_90a = Calcium90TDPre - CalciumTDPre;
 Time_to_90Relax = Calcium90TDPost - CalPeak;
 clear CTDninetyvalue
 CTD90 = Calcium90TDPost - Calcium90TDPre;

 
%%%%% CTD 50
 CTDFiftylinePre = (ones (size(Calcurvexpre,1),1) .*CTDfiftyvalue);
 [CalciumTDPreX50, CalciumTDPreY50] = intersections(Calcurvexpre, Calcurveypre, Calcurvexpre, CTDFiftylinePre, ' robust');
 Calcium50TDPre = min(CalciumTDPreX50);
 CTDFiftylinePost = (ones (size(Calcurvexpost,1),1) .*CTDfiftyvalue);
 [CalciumTDPostX50, CalciumTDPostY50] = intersections(Calcurvexpost, Calcurveypost, Calcurvexpost, CTDFiftylinePost, ' robust');
 Calcium50TDPost = max (CalciumTDPostX50);
 Time_to_50a = Calcium50TDPre - CalciumTDPre; 
 Time_to_50Relax = Calcium50TDPost - CalPeak;
 clear CTDfiftyvalue
 CTD50 = Calcium50TDPost - Calcium50TDPre;
 
 
 %%%%% CTD 10
 CTDTenlinePre = (ones (size(Calcurvexpre,1),1) .*CTDtenvalue);
 [CalciumTDPreX10, CalciumTDPreY10] = intersections(Calcurvexpre, Calcurveypre, Calcurvexpre, CTDTenlinePre, ' robust');
 Calcium10TDPre = min(CalciumTDPreX10);
 CTDTenlinePost = (ones (size(Calcurvexpost,1),1) .*CTDtenvalue);
 [CalciumTDPostX10, CalciumTDPostY10] = intersections(Calcurvexpost, Calcurveypost, Calcurvexpost, CTDTenlinePost, ' robust');
 Calcium10TDPost = max (CalciumTDPostX10);
 Time_to_10a = Calcium10TDPre - CalciumTDPre ;
 Time_to_10Relax = Calcium10TDPost - CalPeak;
 clear CTDtenvalue 
 CTD10 = Calcium10TDPost - Calcium10TDPre;
 
 
%%%%% Peak & Relaxation times
%  T_C_peak = ((CalciumPeakLocation .* frametime) - CalciumTDPre);
 T_C_peak = (CalPeak - CalciumTDPre);

 Time_C_peak = T_C_peak;
%  T_C_relax = ( CalciumTDPost -(CalciumPeakLocation .* frametime) ) ;
 T_C_relax = ( CalciumTDPost - CalPeak ) ;
 

% figure, plot(CalciumTraceXvals,AvgCalciumTrace), hold on
%  plot( Calcurvexpre, CTDlinePre,Calcurvexpost, CTDlinePost)
%  plot(CalciumTDPreX, CalciumTDPreY,'*',CalciumTDPostX, CalciumTDPostY,'*')
%  
%  plot(Calcurvexpre, CTDNinetylinePre,Calcurvexpost, CTDNinetylinePost)
%  plot(CalciumTDPreX90, CalciumTDPreY90,'*',CalciumTDPostX90, CalciumTDPostY90,'*')
%  
%  plot(Calcurvexpre, CTDFiftylinePre,Calcurvexpost, CTDFiftylinePost)
%  plot(CalciumTDPreX50, CalciumTDPreY50,'*',CalciumTDPostX50, CalciumTDPostY50,'*')
%  
%  plot(Calcurvexpre, CTDTenlinePre,Calcurvexpost, CTDTenlinePost)
%  plot(CalciumTDPreX10, CalciumTDPreY10,'*',CalciumTDPostX10, CalciumTDPostY10,'*')
 
%  ContractionRate (file,:)= ContrRate;
%  TheoreticalContactionrRate (file,:)= TheoreticalContrRate;
end

