function [ ] = convert_to_tif()

    Maindatapath = uigetdir('D:/','Select DATA directory for files');
    Maindatapath = ([Maindatapath '/']);
    cd (Maindatapath);

    prompt = {'\bf \fontsize{12} Please enter the format of your videos (eg .avi):',...
        };
    dlgtitle = 'Video format';
    dims = [1 88];
    definput = {''};
    opts.Interpreter = 'tex';
    prt = inputdlg(prompt,dlgtitle,dims,definput, opts);

    frmt = char (prt{1,1} );
    frmt = frmt(1:end);
    fprintf(1,'\n')
    disp (Maindatapath)

    FilePattern = fullfile(Maindatapath,['*',frmt]);
    fileList = dir(FilePattern);
    for item = 1:length(fileList)
       sample_name{item} = fileList(item).name(1:end-4);
       fileList(item).sample_name = fileList(item).name(1:end-4);
    end
    fprintf (1, '>>> Files found %s\n')
    fprintf(1,'\n') 

    unique_file_names = sample_name;
    Nums = numel (unique_file_names);
    fprintf (1, '>>> Number of files to convert:%s\n')
    disp (Nums)
    fprintf(1,'\n') 
    for file=1:length(unique_file_names')
        FileName1 = [unique_file_names{file} frmt];
        FileName = FileName1(1:end);
        fopen (FileName1);
        fprintf (1, '>>> Now reading :        %s\n', FileName)
        [PFinalImage_c1] = single_stack_loader([Maindatapath, unique_file_names{file} frmt]);

    Write_filename = unique_file_names{file};
    WFN = ([Write_filename '.tif']);

    for K=1:length(PFinalImage_c1(1, 1, :))
       imwrite(PFinalImage_c1(:, :, K), WFN, 'WriteMode', 'append');
    end


    end

    fprintf(1,'>>> Videos stored as .tif files in\n' )
    disp (Maindatapath)
    fprintf(1,'\n')

end

