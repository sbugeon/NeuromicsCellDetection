addpath(genpath('C:\Users\bugeon\Documents\GitHub\NeuromicsCellDetection\'))
addpath('C:\Users\bugeon\Documents\MATLAB\GUI Layout Toolbox\layout') % https://fr.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
%% create folders for this animal
clear
% run every time you had new sections
% ================================================================
Animal_ID = 'jm033';
MainPath = ['D:\invivoReg\',Animal_ID]; % path where the registration data will be saved
% ================================================================

mkdir(fullfile(MainPath,'2DSlice')) % folder where to put the slice images
mkdir(fullfile(MainPath,'ZStack','Raw_tiff')) % folder where to put the raw zstacks images
mkdir(fullfile(MainPath,'ZStack','processed'))
cd(MainPath)
% Configure some parameters for slice files
% --------- Set the data file folder -------------
XStep = 1; % Specify X resolution in microns
ZStep = 1; % Specify Z resolution in microns
slice_file_source = fullfile(MainPath,'2DSlice');
SliceNames = dir(slice_file_source);
SliceNames = {SliceNames.name};
% --------- Create the information table ------------
for dd = 3:length(SliceNames)
    i=dd-2;
    FileName{i,1} = SliceNames{dd};
    Channel{i,1} = 'Red';
    [~,name_temp,~] = fileparts(FileName{i,1});
    DataName{i,1} = ['Slice_',name_temp];
    RunPath{i,1} = fullfile(pwd,DataName{i,1});
end
Channel = categorical(Channel);
Channel_reg = 'Red';
slice_files = table(DataName,Channel,FileName,RunPath);
filepathZ = fullfile(MainPath,'ZStack','Raw_tiff');
filepathZ_process = fullfile(MainPath,'ZStack','processed');

CellPoseFolder = 'D:\Data-Analysis\jm33\Slices';

save('run_info.mat','slice_files','slice_file_source','XStep','ZStep','Channel_reg',...
    'filepathZ','filepathZ_process','CellPoseFolder');
fprintf('[run_config] Environment configuration done at %s\n',...
    datestr(now))
fprintf('[run_config] run_info.mat saved at %s\n',pwd)
%% Preprocess z-stack (identifies cell positions in 3D and returns the point cloud coordinates)
clear
load run_info.mat
% ================================================================
XStepZstack = 1; % Specify X resolution in microns (X is medial to lateral for right hemisphere)
YStepZstack = 1; % Specify Y resolution in microns (Y is posterior to anterior)
ZStepZstack = 1; % Specify Z resolution in microns (Z is deep to superficial)
ZStackNames = dir(filepathZ); % names of z-stack files from different session
ZStackNames = ZStackNames(3:end);
% Set parameters for 3D cell detection
Option_detect3.Sigma = [1 1 1]*4; % size of gaussian filter (microns)
Option_detect3.Res0 = 0.04; % sensitivity of cell detection
Option_detect3.SizeLimit = [100,100000]; % size interval of cells
Option_detect3.MedianFilterSize = [1,1,1]*2; % size of median filter
% ================================================================

for i = 1
    %---------  Load Z-stack tif file data and convert into data.mat
    ZStackFileName = ZStackNames(i).name;
    % folder where data files will be saved
    Full_name = fullfile(filepathZ,ZStackFileName);
    filename_output = [ZStackFileName,'.mat'];
    
    data = neuroReg.loadTiff(XStepZstack, ZStepZstack, Full_name,YStepZstack); %original data
    % --------- Detect cells 3D
    % Apply median filter
    dataZ_mid = neuroReg.medfilt3(data,Option_detect3);
    % Detect cells 3D
    [pt_list_vol, pt_area_vol] = neuroReg.detectCells3(dataZ_mid,Option_detect3);
    % Save data
    ZStackInfo.FileName = ZStackFileName;
    ZStackInfo.RunDateTime = datestr(now);
    save(fullfile(filepathZ_process,filename_output),'-v7.3','ZStackInfo','dataZ_mid','Option_detect3','pt_list_vol','pt_area_vol');
    disp(datestr(now));
    fprintf([ZStackFileName,' saved at the folder ']);
    disp(filepathZ_process);
end
%% manual curation of zstack
clear
load run_info.mat

ZStackNames = ZStackNames(3:end);
for i = 1
    %---------  Load Z-stack tif file data and convert into data.mat
    ZStackFileName = ZStackNames(i).name;
    filename_output = [ZStackFileName,'.mat'];
    GUI_StackCuration(fullfile(filepathZ_process,filename_output))
    
%     GUI_StackCuration_light(fullfile(filepathZ_process,filename_output))
end
%%       Preprocess the slice data.
%       - Load data
%       - Detect cells
%       - Manual curation
% preprocess_slice: Detect cells
% Steps:
% - circle an area where the match should roughly be, then press enter
% - Click the picture to delete unwanted points, then press enter
% - Click on points to add and press enter again

% ================================================================
% ------------ Load environment information -------------
clear;
close all
load run_info.mat
% ------------ Settings, filepath and options -----------
Option_detect2.Res0 = 0.06; % threshold for the cell detection
Slice2Run = 1;
% ================================================================

Channel = Channel_reg;
% options for 2D cell detection
Option_detect2.Sigma = [1 1]*2;
Option_detect2.Threshold = 0;
Option_detect2.SizeLimit = [100,10000];
Option_detect2.MedianFilterSize = [2 ,2];
% Set Options for BW cell detection
Option = neuroReg.setOption(Option_detect2);
Option.StepX = 10;
Option.StepD = 8;
Option.Integ = 15;
Option.MagicNumber = 2;
Option.CellRadius = 7;
Option.TransTol = 50;
Option.AngleTol = 4;
Option.MaxPeakNum = 200;
Option.SigmaRender = 15;
for i = Slice2Run % slice number
    % ------------ Set current slice --------------------
    cf = (slice_files.Channel == Channel);
    slice_files_selected = slice_files(cf==1,:);
    SliceFileName = slice_files_selected.FileName{i,1};
    SliceDataName = slice_files_selected.DataName{i,1}; % without .tif
    filepath = fullfile(slice_files_selected.RunPath{i,1});
    source_filePathName = fullfile(slice_file_source,SliceFileName);
    disp([SliceFileName, ' selected.']);
    % -------------- Load slice ----------------------
    dataS = neuroReg.loadTiff(XStep,ZStep,source_filePathName);
    data_slice.value = squeeze(dataS.value(:,1,:));
    data_slice.x = dataS.x;
    data_slice.y = dataS.z;
    disp([SliceFileName, ' Loaded.']);
    
    % first outline the surface of the slice
    neuroReg.plotData2(data_slice);
    title('Circle the potential area for match center...')
    d = drawpolygon();
    ROI_limX = d.Position(:,1);
    ROI_limY = d.Position(:,2);
    
    [pt_list_slice, pt_area_slice] = neuroReg.detectCells2(data_slice,Option_detect2);
    h = gca;
    % then manual curation, adding and deleting cells
    [pt_list_slice,pt_area_slice] = neuroReg.addDelCells2(h,pt_list_slice,pt_area_slice);
    
    SliceInfo.FileName = SliceDataName;
    SliceInfo.RunDateTime = datestr(now);
    disp(Option);
    MatFilePath = fullfile(filepath,[SliceDataName,'.mat']);
    mkdir(filepath)
    save(MatFilePath,'SliceInfo','data_slice',...
        'pt_list_slice','pt_area_slice','Option','Option_detect2','ROI_limX','ROI_limY');
    fprintf('||||||||||||||||||||||||||||||||||||||||||||\n');
    fprintf(['Slice_',SliceDataName,' saved at the folder\n']);
    fprintf('||||||||||||||||||||||||||||||||||||||||||||\n');
    disp(filepath);
    disp(datestr(now));
    close all;
end

%% manual slice curation
clear;
close all
load run_info.mat

GUI_CurateSliceNeuroReg(slice_files.DataName,strrep(slice_file_source,'2DSlice',''))
%% Prepare for correlation calculation
clear;
load run_info.mat

Slice2Run = 1;
Stack2Run = 1;

% ------------- Adjust parameters --------------
% ================================================================
AngleRange = [-20 -13 8;-6 0 7;-14 -10 5];% AngleRange = [Alpha_start Alpha_end Alpha_points; Beta_start...; Gamma_...]
% AngleRange = [-20 -13 8;-2 -2 1;-13 -13 1];
% set a range for the slice position
Option.DepthRange = [0 Inf]; % if there is any assumption on which depth this section is
% or give it for each slice
SlicePosPath = []; % if none, set to []
Range = 50; % tolerance range around the assumed slice position

Option.StepX = 10; % smaller value will give more accurate matches, but are slower
Option.StepD = 20; % smaller value will give more accurate matches, but are slower
Option.Integ =  25; % how much to integrate pixels around the plane for the stack = slice thickness

ScaleF_Y = 1; % Set to 1 usually!!!!!!!!
ScaleF_X = 1; % Set to 1 usually!!!!!!!!

% ================================================================
Option.CellRadius = 10;
Option.MagicNumber = 1.5;
Option.MaxPeakNum = 1000;
Option.TransTol = 100;%50;
Option.AngleTol = 4;%2;
Option.MismatchPen = 1; %1 for neg values on stack img, 0 for neg values on slice img
%
Option.Sigma = [1 1]*2;
Option.Res0 = 0.03;
Option.Threshold = 0;
Option.SizeLimit = [10,Inf];
Option.MedianFilterSize = [2 ,2];
%-------------  Load data -------------
Channel = Channel_reg; % 'Red' for red channel, 'Green' for green channel
% Specify the storage path for preprocessed ZStack files
cf = (slice_files.Channel == Channel);
slice_files_selected = slice_files(cf,:);
SlicesName = slice_files_selected.DataName;
ZStackNames = dir(filepathZ); %names of z-stack files from different session
ZStackNames = ZStackNames(3:end);
fprintf('Preparation ready.\n')
% Run the FFT correlation calculation, find best matches
for j = Stack2Run % loop through stacks
    for i = Slice2Run %loop through sections
        % ----- Set Path -------
        this_slice = SlicesName{i};
        this_ZStack = ZStackNames(j).name;
        this_slice_path = fullfile(pwd,this_slice);
        this_ZStack_file = fullfile(filepathZ_process,[this_ZStack,'.mat']);
        this_ZStack_Proc_file = fullfile(filepathZ_process,[this_ZStack,'_curated.mat']);
        ZStack_file = load(this_ZStack_file);
        Option_detect3 = ZStack_file.Option_detect3;
        dataZ_mid = ZStack_file.dataZ_mid;
        pt_area_vol = ZStack_file.pt_area_vol;
        pt_list_vol = ZStack_file.pt_list_vol;
        if isfile(this_ZStack_Proc_file)
            clear pt_list_vol
            load(this_ZStack_Proc_file)
        end
        
        Slice_file = load(fullfile(this_slice_path,[this_slice,'.mat']));
        % set options
        pt_area_slice = Slice_file.pt_area_slice;
        pt_list_slice = Slice_file.pt_list_slice;
        try
            load(fullfile(this_slice_path,[this_slice,'_curated.mat']));
        catch
            fprintf('\n no curation')
        end
        data_slice = Slice_file.data_slice;
        % find slice approximate position if given
        if ~isempty(SlicePosPath)
            SlicePos = loadSlicePos(SlicePosPath,SlicesName{i}(7:end));
            Option.DepthRange = -[SlicePos+Range SlicePos-Range];
        end
        
        Option.Hist3Flag = 0;
        this_result_path = fullfile(this_slice_path,[this_ZStack,'_xcc']);
        [~] = mkdir(this_result_path);
        RunInfo.SliceName = this_slice;
        RunInfo.ZStackName = this_ZStack;
        RunInfo.AngleRange = AngleRange;
        
        % ----- Calculation ----------
        fprintf('*******************\n');
        disp(this_result_path);
        RunInfo.StartTime = datestr(now);
        disp(datestr(now));
        fprintf('Calculation started...\n')
        data_slice_bw = neuroReg.bwCell2(data_slice,Option,[]);
        h = figure;
        subplot(3,1,1);
        neuroReg.plotData2(data_slice);
        subplot(3,1,2);
        neuroReg.plotData2(data_slice);
        hold on;
        scatter(pt_list_slice(1,:),pt_list_slice(2,:),'r');
        hold off;
        subplot(3,1,3);
        neuroReg.plotData2(data_slice_bw);
        saveas(h, fullfile(this_result_path,'Slice.tif'))
        pt_area_slice = ones(length(pt_list_slice),1)*mean(pt_area_slice);
        % rotate z-stack and find best peaks from FFT
        
        % rescale...
        pt_list_vol(2,:) = pt_list_vol(2,:)*ScaleF_X;
        pt_list_vol(1,:) = pt_list_vol(1,:)*ScaleF_Y;
        
        dataZ_mid.value = imresize3(dataZ_mid.value,[size(dataZ_mid.value,1)*ScaleF_Y  ...
            size(dataZ_mid.value,2)*ScaleF_X size(dataZ_mid.value,3)]);
        dataZ_mid.x = 1:size(dataZ_mid.value,1);
        dataZ_mid.y = 1:size(dataZ_mid.value,2);
        
        % cross-correlation
        tic
        TransTable = neuroReg.rotationCorr3(pt_list_slice,pt_area_slice,...
            pt_list_vol,data_slice,AngleRange,Option,Slice_file.ROI_limX,Slice_file.ROI_limY);
        toc
        % find score for pt cloud registration for each peak
        [ScoreStack , ScoreSlice] = neuroReg.make_icp(TransTable,pt_list_slice, ...
            pt_list_vol,data_slice,dataZ_mid,15,Option,Slice_file.ROI_limX,Slice_file.ROI_limY);
        if ~isempty(TransTable.Intensity)
        TransTable.IntensityXCorr = mat2gray(TransTable.Intensity)+0.1;
        TransTable.ScoreStack = mat2gray(ScoreStack)+0.1;
        TransTable.ScoreSlice = mat2gray(ScoreSlice)+0.1;
        TransTable.Intensity = sqrt(2*TransTable.IntensityXCorr.^2 + TransTable.ScoreStack.^2 + TransTable.ScoreSlice.^2);
        TransTable.Intensity = TransTable.IntensityXCorr + TransTable.ScoreStack + TransTable.ScoreSlice;
        TransTable = sortrows(TransTable,'Intensity','descend');
        end
        RunInfo.EndTime = datestr(now);
        disp(RunInfo.EndTime);
        fprintf('Calculation completed.\n')
        filepathname = fullfile(this_result_path,'result_norm.mat');
        save(filepathname,'RunInfo','TransTable',...
            'Option','Option_detect3',...
            'pt_list_vol','pt_area_vol','pt_list_slice','pt_area_slice');
        disp(RunInfo);
        t1 = datetime(RunInfo.StartTime);
        t2 = datetime(RunInfo.EndTime) ;
        t = t2-t1;
        fprintf('Duration ');
        disp(t)
        fprintf('Option\n ');
        disp(Option);
        fprintf('Result saved!\n---------------\n');
        disp(datestr(now));
        fprintf('Mission Completed!\n*******************\n');
        close all;
    end
end

%% Visualisation of results, adjustement of matches and save figures
clear
load run_info.mat
%
Slice2Run = 1;%1:size(slice_files,1);
Stack2Run = 1;
%
Channel = Channel_reg; % 'Red' for red channel, 'Green' for green channel
ZStackNames = dir(filepathZ); %names of z-stack files from different session
ZStackNames = ZStackNames(3:end);
cf = (slice_files.Channel == Channel);
slice_files_selected = slice_files(cf,:);
SlicesName = slice_files_selected.DataName; % i

ScaleF_Y = 1; % Set to 1 usually!!!!!!!!
ScaleF_X = 1; % Set to 1 usually!!!!!!!!

Visualization = 1; % 0 for point cloud, 1 for image overlay(slower)

for j = Stack2Run% loop through stacks
    for i = Slice2Run% loop through sections
        this_slice = SlicesName{i};
        this_ZStack = ZStackNames(j).name;
        this_slice_path = fullfile(pwd,this_slice);
        this_ZStack_file = fullfile(filepathZ_process,[this_ZStack,'.mat']);
        this_ZStack_Proc_file = fullfile(filepathZ_process,[this_ZStack,'_curated.mat']);
        this_result_path = fullfile(this_slice_path,[this_ZStack,'_xcc']);
        fprintf('Result in %s',this_result_path);
        filepathname = fullfile(this_result_path,'result_norm.mat');
        load(this_ZStack_file); % load zstack
        if isfile(this_ZStack_Proc_file)
            load(this_ZStack_Proc_file)
        end
        
        fprintf('%s Loaded.\n',this_ZStack);
        load(fullfile(this_slice_path,[this_slice,'.mat'])); % load slice
        fprintf('%s Loaded.\n',this_slice);
        load(filepathname); % load result
        fprintf('Result Loaded.\n');
        
        Slice_file = load(fullfile(this_slice_path,[this_slice,'.mat']));
         
        Option.Visualization = Visualization;
        
        % Set DataSets structure to record all the raw data (and binarized slice)
        DataSets.dataZ = dataZ_mid;
        DataSets.data_slice = data_slice;
        pt_area_slice_norm = ones(length(pt_area_slice),1)*mean(pt_area_slice);
        x_lim_slice = [data_slice.x(1),data_slice.x(end)];
        y_lim_slice = [data_slice.y(1),data_slice.y(end)];
        data_slice_bw_low = neuroReg.renderCell2(x_lim_slice,y_lim_slice,Option.StepX,pt_list_slice,pt_area_slice_norm,Option);
        nf = isnan(data_slice_bw_low.value);
        data_slice_bw_low.value = data_slice_bw_low.value - mean(data_slice_bw_low.value(~nf))*Option.MagicNumber;
        data_slice_bw_low.value(nf) = 0;
        DataSets.data_slice_bw_low = data_slice_bw_low;
        
        dataZ_mid.value = imresize3(dataZ_mid.value,[size(dataZ_mid.value,1)*ScaleF_Y  ...
            size(dataZ_mid.value,2)*ScaleF_X size(dataZ_mid.value,3)]);
        dataZ_mid.x = 1:size(dataZ_mid.value,1);
        dataZ_mid.y = 1:size(dataZ_mid.value,2);
        DataSets.dataZ = dataZ_mid;
        
%         TransTable = keepGoodMatches(TransTable,Slice_file.ROI_limX,Slice_file.ROI_limY);
        
        neuroReg.VisTransform3(TransTable,DataSets,pt_list_vol,pt_list_slice,[],Option,this_result_path,'Match_found.mat');
        fprintf('VisTransform\n');
        close all hidden;
        
        % ----- Saves the match found -----
        new_filepathname = fullfile(this_result_path,'Match_found.mat');
        load(new_filepathname);
        h = figure(10086);
        if isempty(TransTable)==0
            neuroReg.plotTransform(h,TransTable(1,:),DataSets,pt_list_vol,pt_list_slice,Option,this_result_path,[],[],1);
        end
        fprintf('*******************\n');
        
        disp(datestr(now));
        fprintf('Mission Completed!\n*******************\n');
    end
end
close all
%% recover all angles and translation for good matches
clear
load run_info.mat
Channel = Channel_reg; % 'Red' for red channel, 'Green' for green channel
ZStackNames = dir(filepathZ); %names of z-stack files from different session
ZStackNames = ZStackNames(3:end);
cf = (slice_files.Channel == Channel);
slice_files_selected = slice_files(cf,:);
SlicesName = slice_files_selected.DataName; % i
filepathZ_proc = fullfile(filepathZ_process,'Proc');% folder where to find raw tiff zstacks
AllTables = table();
for j =1 %  loop through stacks
    for i = 1:length(SlicesName)% loop through sections
        
        this_slice = SlicesName{i};
        this_ZStack = ZStackNames(j).name;
        this_slice_path = fullfile(pwd,this_slice);
        this_ZStack_file = fullfile(filepathZ_process,[this_ZStack,'.mat']);
        this_ZStack_Proc_file = fullfile(filepathZ_process,[this_ZStack,'.mat_curated.mat']);
        this_result_path = fullfile(this_slice_path,[this_ZStack,'_xcc']);
        
        try
            new_filepathname = fullfile(this_result_path,'Match_found.mat');
            load(new_filepathname);
        catch
            TransTable = table();
        end
        if ~isempty(TransTable)
            TransTable.slice =   {this_slice};
            AllTables = [AllTables;TransTable];
        end
        
    end
end
