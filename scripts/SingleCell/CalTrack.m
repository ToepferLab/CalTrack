
% CalTrack_SingleCell

script_folder = cd;
addpath(script_folder);

disp('Please select the bioformats directory');
bf_path=uigetdir('Select bioformats directory'); %  bio-format directory
addpath(bf_path);

frametime = 1/AcquisitionFrequency;

% convert vsi to tif

convert = questdlg('Do you need to convert videos into .tif?','Convert to tif',...
    'Yes','No ', 'Yes');

switch convert
    
    case 'Yes'
    disp('Please select the video data folder');
    convert_to_tif()
    Maindatapath = cd;
    
    case 'No '
    disp('Please select the .tif data directory');
    Maindatapath=uigetdir('Select DATA (.tif) directory');
    cd (Maindatapath);
    addpath(Maindatapath);

end

mkdir calcium_analysis

FilePattern = fullfile(Maindatapath,'*.tif');
fileList = dir(FilePattern);
 
for item = 1:length(fileList)
   sample_name{item} = fileList(item).name(1:end -4);
   fileList(item).sample_name = fileList(item).name(1:end-4);
end
unique_file_names = (unique(sample_name));

for file=1:length(unique_file_names)
    FileName = [unique_file_names{file}, '.tif'];
    fprintf (1, '>>> Now reading :        %s\n', FileName)
    name=strcat(unique_file_names{file},'.tif');
    [IMM(file).PFinalImage_c1] = single_stack_loader(name);
end

for file=1:length(unique_file_names)

    PFinalImage_c1 = IMM(file).PFinalImage_c1;
    img1 = graythresh(PFinalImage_c1);
    IM1 = imbinarize (PFinalImage_c1,img1);
    IM = mean (IM1,3);

    mask1 = imbinarize(IM);
    mask2 = bwareaopen(mask1, 500);
    mask3 = bwmorph(mask2, 'bridge');
    mask4 = imfill(mask3,'holes');
    mask_kernel = strel('sphere', 10);
    FinalMask = imdilate(mask4, mask_kernel);
    
    PFI =  double(PFinalImage_c1);
    MaskedCalcium =   FinalMask .* PFI;
    Calcium = squeeze(sum(sum(MaskedCalcium)));
   
    Calcium_Traces(file).data = Calcium;
    Calcium_Traces(file).name = unique_file_names{file};
    
end 

cd('calcium_analysis')
save('Calcium_Traces.mat','Calcium_Traces')
disp(strcat('Done! Results saved in folder "calcium_analysis" in : ',pwd))
results_folder = pwd;
cd(script_folder);