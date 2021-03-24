function xlsx_to_mat()

[file,path] = uigetfile('*.xlsx');

names = sheetnames(fullfile(path,file));
for i = 1:length(names)
    data = xlsread(fullfile(path,file),names{i});
    data_matrix(i).time = data(:,1);
    data_matrix(i).calcium = data(:,2:4:end);
end

% Calcium_Traces.mat 

cd(path)
for i = 1:length(data_matrix)
    Calcium_Traces = [];
    folder_mat = ['Calcium_Traces_',names{i}];
    mkdir(folder_mat)
    for j = 1:size(data_matrix(i).calcium,2)
        Calcium_Traces(j).name = ['cell',num2str(j)];
        Calcium_Traces(j).data = data_matrix(i).calcium(:,j);
        save(fullfile(folder_mat,'Calcium_Traces.mat'),'Calcium_Traces')
    end
end
