
function [] = FR_to_calcium(results_folder)

    global scripts_folder

    load('Calcium_Traces.mat');

    files = dir('*.xlsx');
    if ismember('CalibrationCurve.xlsx',{files.name})

        [param_calib_curve]= (xlsread('CalibrationCurve.xlsx'));

        f = @(fr) (-param_calib_curve(1) + param_calib_curve(2) *fr.^param_calib_curve(3));

    elseif ismember('CalibrationPoints.xlsx',{files.name})

        [points_calib_curve]= xlsread('CalibrationPoints.xlsx');
        ca = points_calib_curve(:,1);
        fr = points_calib_curve(:,2);  
        [fitresult, ~] = createFit(fr, ca);
        f = @(fr) (fitresult.a + fitresult.b *fr.^fitresult.c);

    end

    cd(scripts_folder)

    for i = 1:length(Calcium_Traces)
        FR = [];
        FR = Calcium_Traces(i).data;
        Calcium_Traces_quantitative(i).data = f(FR);
        Calcium_Traces_quantitative(i).name = Calcium_Traces(i).name;
    end

    cd(results_folder)
    save('Calcium_Traces_quantitative.mat','Calcium_Traces_quantitative')

end

function [fitresult, gof] = createFit(fr, ca)

    [xData, yData] = prepareCurveData( fr, ca );

    ft = fittype( 'a+b*x^c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [0.0682503973370361 0.99380184306144 0.605633843187886];

    [fitresult, gof] = fit( xData, yData, ft, opts );

end
