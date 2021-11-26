addpath('C:\Users\Stephane\Downloads\NeuroReg_pipeline_2\NeuroReg_pipeline\NeuroReg2')
addpath('C:\Users\Stephane\Downloads\NeuroReg_pipeline_2\NeuroReg_pipeline\NeuroReg2\polyhedraCut')
addpath('C:\Users\Stephane\OneDrive - University College London\Documents\MATLAB\toolbox\layout')
%% Configure some parameters for slice files
clear;
close all;
% --------- Set the data file folder -------------
Animal_ID = 'SB026';
XStep = 1; % Specify X resolution in microns (X is medial to lateral for rigth hemisphere)
ZStep = 1; % Specify Z resolution in microns (Z is deep to superficial)
slice_file_source = fullfile(pwd,'2DSlice');% Specify the path for slice files storage
slice_positions = 1:10; % Specify order of slices
fmt = ['C%d-',Animal_ID,'-%02d.tif']; % The format of the tif files.
% --------- Create the information table ------------
for i = 1:length(slice_positions)
    FileName{i,1} = sprintf(fmt,1,i);
    Channel{i,1} = 'Red';
    [~,name_temp,~] = fileparts(FileName{i,1});
    DataName{i,1} = ['Slice_',name_temp];
    RunPath{i,1} = fullfile(pwd,DataName{i,1});
end
Channel = categorical(Channel);
Channel_reg = 'Red';
slice_files = table(DataName,Channel,FileName,RunPath);
filepathZ = fullfile(pwd,'\ZStack\Raw_tiff\');% folder where to find raw tiff zstacks
filepathZ_process = fullfile(pwd,'\ZStack\processed\');

save('run_info.mat','slice_files','slice_file_source','XStep','ZStep','Channel_reg',...
    'filepathZ','filepathZ_process');
fprintf('[run_config] Environment configuration done at %s\n',...
    datestr(now))
fprintf('[run_config] run_info.mat saved at %s\n',pwd)

%% Preprocess z-stack (identifies cell positions in 3D and returns the point cloud coordinates)
clear;
load run_info.mat
XStepZstack = 1; % Specify X resolution in microns (X is medial to lateral for right hemisphere)
YStepZstack = 1; % Specify Y resolution in microns (Y is posterior to anterior)
ZStepZstack = 1; % Specify Z resolution in microns (Z is deep to superficial)
ZStackNames = dir(filepathZ); % names of z-stack files from different session
ZStackNames = ZStackNames(3:end);
% Set parameters for 3D cell detection
Option_detect3.Sigma = [1 1 1]*4; % size of gaussian filter (microns)
Option_detect3.Res0 = 0.01; % sensitivity of cell detection
Option_detect3.SizeLimit = [700,1500]; % size interval of cells
Option_detect3.MedianFilterSize = [1,1,1]*2; % size of median filter

for i = 1:length(ZStackNames)
    %---------  Load Z-stack tif file data and convert into data.mat
    ZStackFileName = ZStackNames(i).name;
    % folder where data files will be saved
    Full_name = fullfile(filepathZ,ZStackFileName);
    data = neuroReg.loadTiff(XStepZstack, ZStepZstack, Full_name,YStepZstack); %original data
    filename_output = [ZStackFileName,'.mat'];
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

%%       Preprocess the slice data.
%       - Load data
%       - Detect cells
%       - Manual curation
% ================================================================
% ------------ Load environment information -------------
clear;
close all
load run_info.mat
% ------------ Settings, filepath and options -----------
% THRESHOLD = 0.2; %%%%%%%%%%%%%%% set threshold for thresholding red image
Channel = Channel_reg;
% options for 2D cell detection
Option_detect2.Sigma = [1 1]*2;
Option_detect2.Res0 = 0.05;
Option_detect2.Threshold = 0;
Option_detect2.SizeLimit = [100,10000];
Option_detect2.MedianFilterSize = [2 ,2];
% Set Options for BW cell detection
Option = neuroReg.setOption(Option_detect2);
Option.StepX = 10;
Option.StepD = 8;
Option.Integ = 16;
Option.MagicNumber = 2;
Option.CellRadius = 7;
Option.TransTol = 50;
Option.AngleTol = 4;
Option.MaxPeakNum = 200;
Option.SigmaRender = 15;

for i = 1
    
    % ------------ Set current slice --------------------
    cf = (slice_files.Channel == Channel);
    slice_files_selected = slice_files(cf==1,:);
    SliceFileName = slice_files_selected.FileName{1,1};
    SliceDataName = slice_files_selected.DataName{1,1}; % without .tif
    filepath = fullfile(slice_files_selected.RunPath{1,1});
    source_filePathName = fullfile(slice_file_source,SliceFileName);
    disp([SliceFileName, ' selected.']);
    % -------------- Load slice ----------------------
    dataS = neuroReg.loadTiff(XStep,ZStep,source_filePathName);
    data_slice.value = squeeze(dataS.value(:,1,:));
    data_slice.x = dataS.x;
    data_slice.y = dataS.z;
    disp([SliceFileName, ' Loaded.']);
    % preprocess_slice: Detect cells
    % Steps:
    % - circle an area where the match should roughly be
    % - Click the picture to select points to be excluded.
    % - Press enter to apply.
    % - Click on points to add and press enter again until satisfied.
    % - Select nothing and press enter to finish.
    
    % first outline the area where the peak should be
    neuroReg.plotData2(data_slice);
    title('Circle the potential area for match center...')
    [ROI_limX , ROI_limY] = getpts();
    % then automatic cell detection
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

%% Prepare for correlation calculation
clear;
load run_info.mat
% ------------- Adjust parameters --------------
AngleRange = [-6 6 13;-5 0 6;0 7 8];% AngleRange = [Alpha_start Alpha_end Alpha_points; Beta_start...; Gamma_...]
Option.DepthRange = [50 250]; % adjust the depth at which this slice should be
Option.StepX = 9; % larger value will give more accurate matches, but are slower
Option.StepD = 10;
Option.Integ =  20; % how much to integrate pixels around the plane for the stack
Option.CellRadius = 9;
Option.MagicNumber = 1.5;
Option.MaxPeakNum = 1000;
Option.TransTol = 50;
Option.AngleTol = 2;
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
filepathZ_proc = fullfile(filepathZ_process,'Proc');% folder where to find raw tiff zstacks
% Run the FFT correlation calculation, find best matches
for j =1 % loop through stacks
    for i =1 % loop through sections
        % ----- Set Path -------
        this_slice = SlicesName{i};
        this_ZStack = ZStackNames(j).name;
        this_slice_path = fullfile(pwd,this_slice);
        this_ZStack_file = fullfile(filepathZ_process,[this_ZStack,'.mat']);
        this_ZStack_Proc_file = fullfile(filepathZ_proc,[this_ZStack,'.mat']);
        ZStack_file = load(this_ZStack_file);
        Option_detect3 = ZStack_file.Option_detect3;
        dataZ_mid = ZStack_file.dataZ_mid;
        pt_area_vol = ZStack_file.pt_area_vol;
        pt_list_vol = ZStack_file.pt_list_vol;
        if isfile(this_ZStack_Proc_file)
            load(this_ZStack_Proc_file)
            pt_list_vol = SliceROI;
        end
        
        Slice_file = load(fullfile(this_slice_path,[this_slice,'.mat']));
        % set options
        
        pt_area_slice = Slice_file.pt_area_slice;
        pt_list_slice = Slice_file.pt_list_slice;
        data_slice = Slice_file.data_slice;
        
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
        tic
        pt_area_slice = ones(length(pt_area_slice),1)*mean(pt_area_slice);
        % rotate z-stack and find best peaks from FFT
        % cross-correlation
        TransTable = neuroReg.rotationCorr3(pt_list_slice,pt_area_slice,...
            pt_list_vol,data_slice,AngleRange,Option,Slice_file.ROI_limX,Slice_file.ROI_limY);
        % find score for pt cloud registration for each peak
        [ScoreStack , ScoreSlice] = neuroReg.make_icp(TransTable,pt_list_slice, ...
            pt_list_vol,data_slice,dataZ_mid,15);
        TransTable.IntensityXCorr = mat2gray(TransTable.Intensity)+0.1;
        TransTable.ScoreStack= mat2gray(ScoreStack)+0.1;
        TransTable.ScoreSlice= mat2gray(ScoreSlice)+0.1;
        TransTable.Intensity= sqrt(2*TransTable.IntensityXCorr.^2 + TransTable.ScoreStack.^2 + TransTable.ScoreSlice.^2);
        TransTable.Intensity= TransTable.IntensityXCorr + TransTable.ScoreStack + TransTable.ScoreSlice;
        TransTable = sortrows(TransTable,'Intensity','descend');
        toc
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
Channel = Channel_reg; % 'Red' for red channel, 'Green' for green channel
ZStackNames = dir(filepathZ); %names of z-stack files from different session
ZStackNames = ZStackNames(3:end);
cf = (slice_files.Channel == Channel);
slice_files_selected = slice_files(cf,:);
SlicesName = slice_files_selected.DataName; % i
filepathZ_proc = fullfile(filepathZ_process,'Proc');% folder where to find raw tiff zstacks

for j =1
    
    for i =1
        
        this_slice = SlicesName{i};
        this_ZStack = ZStackNames(j).name;
        this_slice_path = fullfile(pwd,this_slice);
        this_ZStack_file = fullfile(filepathZ_process,[this_ZStack,'.mat']);
        this_ZStack_Proc_file = fullfile(filepathZ_proc,[this_ZStack,'.mat']);
        this_result_path = fullfile(this_slice_path,[this_ZStack,'_xcc']);
        fprintf('Result in %s',this_result_path);
        filepathname = fullfile(this_result_path,'result_norm.mat');
        load(this_ZStack_file); % load zstack
        if isfile(this_ZStack_Proc_file)
            load(this_ZStack_Proc_file)
            pt_list_vol = SliceROI;
        end
        
        fprintf('%s Loaded.\n',this_ZStack);
        load(fullfile(this_slice_path,[this_slice,'.mat'])); % load slice
        fprintf('%s Loaded.\n',this_slice);
        load(filepathname); % load result
        fprintf('Result Loaded.\n');
        
        Option.Visualization = 1; % 0 for point cloud, 1 for image overlay
        
        % neuroReg.VisTransform
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
        neuroReg.VisTransform3(TransTable,DataSets,pt_list_vol,pt_list_slice,[],Option,this_result_path);
        fprintf('VisTransform\n');
        close all hidden;
        
        % ----- Saves the match found --------
        new_filepathname = fullfile(this_result_path,'Match_found.mat');
        load(new_filepathname);
        h = figure(10086);
        if isempty(TransTable)==0
            neuroReg.plotTransform(h,TransTable(1,:),DataSets,pt_list_vol,pt_list_slice,Option,this_result_path);
        end
        fprintf('*******************\n');
        
        disp(datestr(now));
        fprintf('Mission Completed!\n*******************\n');
    end
end