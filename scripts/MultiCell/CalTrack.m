
% CalTrack_MultiCell

script_folder = cd;
addpath(cd);

disp('Please select the movies directory');

rootPath = uigetdir(path, 'Please select the movies directory');
cd(rootPath)

mkdir calcium_analysis


minArea = 45;

minElongationParameter = 3;

mov = questdlg('Do you have .m4v files?','Please select video file format',...
    'Yes','No ', 'Yes');

switch mov
    
    case 'Yes'
    videoPaths = dir('*.m4v');

    for video = 1: length(videoPaths)

        fileName = videoPaths(video).name;
        fprintf(' >>> analyzing video:  %s \n', fileName)
        cap = [];
        cap = VideoReader(fileName); 

        h = cap.Height;
        w = cap.Width;
        n = cap.NumFrames;

        V = (zeros(h,w,n));
        vidFrame_g = [];

        for i = 1:n
            vidFrame = read(cap,i);
            vidFrame_g(:,:) = rgb2gray(vidFrame);
            tif_file = fullfile(rootPath,strcat(fileName,'_TimeLapse.tif'));
            V(:,:,i) = double(vidFrame_g)/255;
        end

        videos(video).data = V;
        
    end

    case 'No ' 

    tif = questdlg('Do you have .tif files?','Please select video file format',...
    'Yes','No ', 'Yes');

    switch tif
        case'Yes'
        disp('Please select the bio-format directory');
        bf_path=uigetdir('D:/','Select bf directory');
        addpath(bf_path);
        videoPaths = dir('*.tif*');

        for item = 1:length(videoPaths)
           sample_name{item} = videoPaths(item).name(1:end -4);
           videoPaths(item).sample_name = videoPaths(item).name(1:end-4);
        end
        
        unique_file_names = (unique(sample_name));

        for file=1:length(unique_file_names)
            FileName = [unique_file_names{file}, '.tif'];
            fprintf (1, '>>> Now reading :        %s\n', FileName)
            name=strcat(unique_file_names{file},'.tif');
            videos(file).data = single_stack_loader(name);
            videos(file).name = name;

        end
    
        case 'No '
        disp('Please select the video folder to convert to tif');
        convert_to_tif()
        
    end
    
end


for video = 1: length(videoPaths)

    V = videos(video).data;
    Calcium(video).name = videos(video).name;

    I(:,:) = rescale(mean(V,3));
    
    eqI = adapthisteq(I);
    gaus2 = imgaussfilt(eqI, 2.5);
    gaus1 = imgaussfilt(eqI, 0.2);
    imgf = gaus2 - gaus1;
    neg = adapthisteq(adapthisteq(imgf));
    negim = imbinarize(neg);
    complim = imcomplement(negim);    
    imag = eqI .* complim; 
 
    eq_eqI = adapthisteq (imag);
    
    eqI_b = imbinarize(eq_eqI);
    
    SE = strel('disk',1);
    ProjFgMask = imerode(eqI_b,SE);
    sz = 8*(size(ProjFgMask));
    im = imresize(ProjFgMask,sz);
    SE = strel ('disk',5);
    IM = imdilate(im,SE);
    rs_im = imresize(IM, size(ProjFgMask));
    
    fgmc = imregionalmax(ProjFgMask);
    fgm = imclearborder(fgmc);

    I2 = labeloverlay(I,fgm);
%     I3 = bwlabel(fgm);

    labeled = labelmatrix(bwconncomp(fgm));

    stats = regionprops(labeled,'all');

    X = zeros(size(labeled));
    Y = zeros(size(labeled));
    count = 0;
%     j=1;
    for i=1:length(stats)
            mj = [stats(i).MajorAxisLength];
            mn = [stats(i).MinorAxisLength];
            ar = [stats(i).Area];
            if (ar > minArea) && (mj/mn > minElongationParameter)

                count = count+1;
                for l=1:size(stats(i).PixelList,1)
                    X(stats(i).PixelList(l,2),stats(i).PixelList(l,1)) = count;
                    Y(stats(i).PixelList(l,2),stats(i).PixelList(l,1)) = 0.5;
                end
                
            else

                for l=1:size(stats(i).PixelList,1)
                    Y(stats(i).PixelList(l,1),stats(i).PixelList(l,2)) = 0.25;
                end

            end
    end

    figure('Name',strcat('video',num2str(video))), clf
    map = get(gca,'ColorOrder');
    map = [map; map; map; map];
    subplot(1,2,1)
    I3 = labeloverlay(I,X,'Colormap',map);

    imshow(I3)
    title({'Selected Cells';'Superimposed on Original Image'})   
    hold on
    s = regionprops(X,'Centroid');
    for k = 1:numel(s)
        c = s(k).Centroid;
        text(c(1), c(2), sprintf('%d', k), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle','Color','w','LineWidth',12);
    end
    hold off

    for i=1:size(V,3)
        v(:,:) = V(:,:,i);
        Calcium(video).mean_intensities(:,i) = struct2array(regionprops(X,v,'MeanIntensity'))';
    end
    
    figure(video)
    subplot(1,2,2)
    for i=1:size(Calcium(video).mean_intensities,1)
            plot(Calcium(video).mean_intensities(i,:)), hold on
            
            text(200+10*i, Calcium(video).mean_intensities(i,1), sprintf('%d', i), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle','Color','k','LineWidth',12);

    end
    title('Mean Intensity for Selected Cells')

end

remove =questdlg('Would you like to remove any cell?','Filter undesired cells',...
    'Yes','No ', 'Yes');

switch remove
    case 'Yes'
    prompt = {'\bf \fontsize{12} Please enter the video number for those videos that contain the cells to be removed (from figures title): e.g. [1 2 3] (if more than one cell from the same video, than repeat the video number)',...
    '\bf \fontsize{12} Please enter the cell number for the cells to be removed (subplots title): e.g. 1 2 3 (their order needs to match the videos order and when you want to remove more than 1 cell from the same video, report the cell number from the largest to the smallest'};
    dlgtitle = 'Remove cells';
    dims = [1 120];
    definput = {'','', '', '','','','',''};
    opts.Interpreter = 'tex';
    RC = inputdlg(prompt,dlgtitle,dims,definput, opts);
    rem_vid = str2num(RC{1});
    rem_cell = str2num(RC{2});

    r = [rem_vid' rem_cell'];
      
    for i=1:length(rem_vid)
        Calcium(rem_vid(i)).mean_intensities(rem_cell(i),:)=[];
    end
end

j=1;
for i=1:length(Calcium)
    if (isempty(Calcium(i).mean_intensities)==1)
        id(j)=i;
        j=j+1;
    end
end

if exist('id') == 1
    Calcium(id) = [];
end


cd('calcium_analysis')
switch mov
    case 'Yes'
        for i=1:size(videoPaths,1)
            Cell = Calcium(i).mean_intensities';
            name = {videoPaths.name}';
            writematrix(Cell,fullfile('CalciumData_MultiCell.xlsx'),'Sheet',Calcium(i).name)
        end
end

switch tif
    case 'Yes'
    for i=1:size(videoPaths,1)
        Cell = Calcium(i).mean_intensities';
        name = {videoPaths.name}';
        writematrix(Cell,fullfile('CalciumData_MultiCell.xlsx'),'Sheet',Calcium(i).name)
    end
end

warning('off')

results_folder = pwd;
disp(strcat('Done! Results saved in: ',pwd))
