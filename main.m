clear, clc, close all 

set(0,'defaultlinelinewidth',3);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultTextFontSize',18);

% tic

% -------------------- MAIN SCRIPT FOR CALTRACK -------------------- %

global AcquisitionFrequency
global cutoff
global PacingFrequency
global scripts_folder

main_folder  = cd;

MSG = questdlg('Which type of data are you running?','Please select function',...
    'One cell per video','Multiple cells per video', 'One cell per video');

% user input


if isfile('savedAnalysisParameters.mat')== 0
    
    prompt = {'\bf \fontsize{12} Please enter acquisition frequency (fps):',...
        '\bf \fontsize{12} Please enter pacing frequency (Hz):'...
        '\bf \fontsize{12} Please enter the number of frames to be removed at the beginning of the movie:'...
        '\bf \fontsize{12} Please enter the number of frames to be analysed (if all, enter 0):'...
        '\bf \fontsize{12} Please enter the cutoff for extra beat detection (as a % of the pacing rate):'...
        '\bf \fontsize{12} Is the signal to noise ratio good enough to allow single beat computation? (y/n):'...  
        '\bf \fontsize{12} Would you like to extract traces from videos? (y/n):'...
        '\bf \fontsize{12} Would you like to provide the stimulus frames to segment transients? (y/n):'...
        '\bf \fontsize{12} Would you like to convert fluorescence in [Ca]? (y/n):'...
    };

    dlgtitle = 'Analysis parameters';
    dims = [1 88];
    definput = {'','', '', '','','','','',''};
    opts.Interpreter = 'tex';
    answers = inputdlg(prompt,dlgtitle,dims,definput, opts);


    AcquisitionFrequency = str2num(answers{1,1});
    PacingFrequency = str2num(answers{2,1});
    N = str2num(answers{3,1});
    M = str2num(answers{4,1});
    cutoff = str2num(answers{5,1});
    cutoff = (1/PacingFrequency)-(cutoff/100)*(1/PacingFrequency);
    param_single_beat =answers{6,1};
    conversion = answers{7,1};
    segmentation = answers{8,1};
    quantitative_data = answers{9,1};


    savedAnalysisParameters.AcquisitionFrequency = AcquisitionFrequency;
    savedAnalysisParameters.PacingFrequency = PacingFrequency;
    savedAnalysisParameters.N = N;
    savedAnalysisParameters.M = M;
    savedAnalysisParameters.cutoff = cutoff;
    savedAnalysisParameters.param_single_beat = param_single_beat;
    savedAnalysisParameters.conversion = conversion;
    savedAnalysisParameters.segmentation = segmentation;
    savedAnalysisParameters.quantitative_data = quantitative_data;

    prompt = {'\bf \fontsize{12} Would you like to save the analysis parameters for later (y/n)?',...
    };

    dlgtitle = 'Save analysis parameters';
    dims = [1 88];
    definput = {''};
    opts.Interpreter = 'tex';
    answers = inputdlg(prompt,dlgtitle,dims,definput, opts);
    switch answers{1}
        case 'y'
            save('savedAnalysisParameters','savedAnalysisParameters')
    end
    
else
    
    prompt = {'\bf \fontsize{12} Would you like to load previoulsy saved analysis parameters (y/n)?',...
    };

    dlgtitle = 'Load analysis parameters';
    dims = [1 88];
    definput = {''};
    opts.Interpreter = 'tex';
    answers = inputdlg(prompt,dlgtitle,dims,definput, opts);
    
    switch answers{1}
        case 'y'
            load('savedAnalysisParameters.mat')
            AcquisitionFrequency = savedAnalysisParameters.AcquisitionFrequency;
            PacingFrequency = savedAnalysisParameters.PacingFrequency;
            N = savedAnalysisParameters.N;
            M = savedAnalysisParameters.M;
            cutoff = savedAnalysisParameters.cutoff;
            param_single_beat = savedAnalysisParameters.param_single_beat;
            conversion = savedAnalysisParameters.conversion;
            segmentation = savedAnalysisParameters.segmentation;
            quantitative_data = savedAnalysisParameters.quantitative_data;
        case 'n'
            prompt = {'\bf \fontsize{12} Please enter acquisition frequency (fps):',...
            '\bf \fontsize{12} Please enter pacing frequency (Hz):'...
            '\bf \fontsize{12} Please enter the number of frames to be removed at the beginning of the movie:'...
            '\bf \fontsize{12} Please enter the number of frames to be analysed (if all, enter 0):'...
            '\bf \fontsize{12} Please enter the cutoff for extra beat detection (as a % of the pacing rate):'...
            '\bf \fontsize{12} Is the signal to noise ratio good enough to allow single beat computation? (y/n):'...  
            '\bf \fontsize{12} Would you like to extract traces from videos? (y/n):'...
            '\bf \fontsize{12} Would you like to provide the stimulus frames to segment transients? (y/n):'...
            '\bf \fontsize{12} Would you like to convert fluorescence in [Ca]? (y/n):'...
            };

            dlgtitle = 'Analysis parameters';
            dims = [1 88];
            definput = {'','', '', '','','','','',''};
            opts.Interpreter = 'tex';
            answers = inputdlg(prompt,dlgtitle,dims,definput, opts);


            AcquisitionFrequency = str2num(answers{1,1});
            PacingFrequency = str2num(answers{2,1});
            N = str2num(answers{3,1});
            M = str2num(answers{4,1});
            cutoff = str2num(answers{5,1});
            cutoff = (1/PacingFrequency)-(cutoff/100)*(1/PacingFrequency);
            param_single_beat =answers{6,1};
            conversion = answers{7,1};
            segmentation = answers{8,1};
            quantitative_data = answers{9,1};


            savedAnalysisParameters.AcquisitionFrequency = AcquisitionFrequency;
            savedAnalysisParameters.PacingFrequency = PacingFrequency;
            savedAnalysisParameters.N = N;
            savedAnalysisParameters.M = M;
            savedAnalysisParameters.cutoff = cutoff;
            savedAnalysisParameters.param_single_beat = param_single_beat;
            savedAnalysisParameters.conversion = conversion;
            savedAnalysisParameters.segmentation = segmentation;
            savedAnalysisParameters.quantitative_data = quantitative_data;

            prompt = {'\bf \fontsize{12} Would you like to save the analysis parameters for later (y/n)?',...
            };

            dlgtitle = 'Save analysis parameters';
            dims = [1 88];
            definput = {''};
            opts.Interpreter = 'tex';
            answers = inputdlg(prompt,dlgtitle,dims,definput, opts);
            switch answers{1}
                case 'y'
                    save('savedAnalysisParameters','savedAnalysisParameters')
            end
    end
end


switch MSG
    case 'One cell per video'
        scripts_folder = ([main_folder,'/scripts/SingleCell']);
        addpath(scripts_folder)
        cd(scripts_folder)
    case 'Multiple cells per video'
        scripts_folder = ([main_folder,'/scripts/MultiCell']);
        addpath(scripts_folder)
        cd(scripts_folder)
end    


%% scripts

atLeast1Beat = 1;

% CalTrack

if conversion == 'y'
    CalTrack
end


% DataAnalysis

switch MSG
    case 'One cell per video'
    
    if exist('results_folder')~=0

        if quantitative_data == 'y'
            copyfile *.xlsx results_folder
            cd(results_folder)
            FR_to_calcium(results_folder)
        else
            cd(results_folder)

        end

    else

        fprintf(1,'\n')
        disp('Please select the folder with the calcium traces')
        results_folder = uigetdir('folder with calcium data?');
        cd(results_folder)

        if quantitative_data == 'y'
            FR_to_calcium(results_folder)
        end

    end
    
    case 'Multiple cells per video'
    
        if exist('results_folder')~=0

        cd(results_folder)

        else

            fprintf(1,'\n')
            disp('Please select the folder with the calcium traces')
            results_folder = uigetdir('folder with calcium data?');
            cd(results_folder)

        end
    
end

DataAnalysis

% PlotParameters

if atLeast1Beat~=0
    show_plots =questdlg('Would you like to plot the results?','Result plotting',...
        'Yes','No ', 'Yes');
    if show_plots == 'Yes'  
        PlotParameters
    end
end

% Analysis of non-regular traces (non-adherence to pacing or non-regular
% single beat decay)
switch MSG
    case 'One cell per video'
    nonRegularAnalysis =questdlg('Would you like to run a separate analysis for non-regular traces?','analysis',...
        'Yes','No ', 'Yes');

    if nonRegularAnalysis == 'Yes'  
        NonRegularBeatsAnalysis
    end
end

cd(main_folder)
rmpath(scripts_folder)

%toc