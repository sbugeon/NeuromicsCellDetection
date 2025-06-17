% import cell detection for brainbow dataset for neuroReg
ColorLab = {'R','G','B','RG','RB','GB'}; % colors contained in the 6 labels from ex vivo stack
ScaleF = 3.53699395067473;

Animal_ID = 'brainbowZ79';
MainPath = ['D:\invivoReg\',Animal_ID]; % path where the registration data will be saved
SliceF_bb = fullfile(MainPath,'alignement in vivo ex vivo\invivo');
StackF_bb = fullfile(MainPath,'alignement in vivo ex vivo\exvivo');

StackCC = readtable(fullfile(StackF_bb,'coor.csv'));
StackCC.X = StackCC.X / ScaleF;
StackCC.Y = StackCC.Y / ScaleF;

load run_info.mat
ZStackNames = dir(filepathZ); % names of z-stack files from different session
ZStackNames = ZStackNames(3:end);

for i = 1:length(ZStackNames)
    pt_list_vol0 = [];
    ZStackFileName = ZStackNames(i).name;
    filename_output = [ZStackFileName,'.mat'];
    filename_output_cur = [ZStackFileName,'_curated.mat'];
    F = fullfile(filepathZ_process,filename_output);
    s = strfind(ZStackFileName,'.tif');
    C = ZStackFileName(s-1);
    ThisL = find(contains(ColorLab,C));
%     ThisL = find(strcmp(ColorLab,C));
    ThisC = StackCC(ismember(StackCC.label,ThisL),1:3);
    pt_list_vol0(1,:) = ThisC.Y';
    pt_list_vol0(2,:) = ThisC.Z';
    pt_list_vol0(3,:) = ThisC.X';
    
    load(F);
    
    % Visualize
    pt_area = ones(1,size(pt_list_vol0,2));
    assignin('base','data_temp_749',dataZ_mid);
    neuroReg.PlotSlices('data_temp_749','y',pt_list_vol0,pt_area,5,[]);
    pt_list_vol = pt_list_vol0;
    save(fullfile(filepathZ_process,filename_output_cur),'pt_list_vol')
    
end
Channel = Channel_reg;

Scale = 470/512;
for i=1:size(slice_files,1)
    cf = (slice_files.Channel == Channel);
    slice_files_selected = slice_files(cf==1,:);
    filepath = fullfile(slice_files_selected.RunPath{i,1});
    SliceDataName = slice_files_selected.DataName{i,1}; % without .tif
    SliceN = strrep(slice_files.FileName{i},'.tif','');
    MatPath = fullfile(filepath,[SliceDataName,'.mat']);
    load(MatPath)
    T = readtable(fullfile(SliceF_bb,['centroids_',SliceN,'.csv']));
    
    figure
    neuroReg.plotData2(data_slice);
    hold on
    %      scatter(pt_list_slice(1,:),pt_list_slice(2,:),'or')
    pt_list_slice0=[];
    pt_list_slice0(1,:) = T.centroid_x'*Scale;
    pt_list_slice0(2,:) = size(data_slice.value,1) - T.centroid_y'*Scale;
    scatter(pt_list_slice0(1,:),pt_list_slice0(2,:),'ob')
    pt_list_slice = pt_list_slice0;
    pt_area_slice = ones(size(pt_list_slice,2),1);
    save(MatPath,'SliceInfo','data_slice',...
        'pt_list_slice','pt_area_slice','Option','Option_detect2','ROI_limX','ROI_limY');
end


% RGB=[];
% for i=1:6
%     gg = StackCC.coor.label == i;
%     RGB(i,:) = mean([StackCC.coor.R(gg) StackCC.coor.G(gg) StackCC.coor.B(gg)]);
% end
% RGB = RGB./max(RGB);
% figure
% bar(RGB)

%% for JC
% - dans coor.mat, X et Y ont la meme valeur, du coup j'ai retrouve le
% scaling factor a partir de X, et je l'ai applique a Y en reprenant le
% xlsx d'origin

% - j'ai essaye de retrouver les couleurs correspondant aux labels:
% {'R','G','B','RG','RB','GB'}, est-ce correct?

