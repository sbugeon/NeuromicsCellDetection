function TransTable = rotationCorr3_piece(pt_list_slice,pt_area_slice,pt_list_vol,data_slice,AngleRange,Option,X,Y)
% rotationCorr3 returns the transformation parameters and the corresponding
% correlation functions.
% TransTable = ...
% rotationCorr3(pt_list_slice,pt_list_vol,data_slice,AngleRange,Option)
% The transfermation matrix is from volume to slice.
%
% ---- OUTPUT ----
% TransTable record the possible transformation parameters and the
% coorelation intensity.
% n-by-(1+6) table with column names:
% Intensity, Angle, Translation
% [{intensity},{alpha,beta,gamma},{tx,ty,yz}]
% The transformation is from volume to cut.
% Corr record the corresponding correlation functions.
% 1-by-n
%
% --- INPUT ----
% pt_list is 3xN array. Coordinations of the cells in 3D volume
% data_slice is the *2D* slice (with dataS.x, dataS.y and dataS.value)
% AngleRange should be a 3x3 matrix with the form
% [Alpha_Start, Alpha_End, Number; Beta_Start, Beta_End, Beta_Number; Gamma...]
% alpha: rotation around y-axis
% beta: rotation around z-axis
% gamma: rotation around x-axis
% Convention:
% R = Rx*Rz*Ry
% Y = R*X
% **********
% * OPTION *
% **********
% Option is a structure specifying the options used in calculating
% correlation function.
% Fields not specified will be set as default.
% Below is the description of each field:
% (I write it with enormous patience because tomorrow morning I will forget
% them all)
% ---------
% For peak record
% 'TransTol'    >> Translation tolerance to divide two peaks in XCC
% 'AngleTol'    >> Angle tolerance to divide two peaks in XCC
% 'MaxPeakNum'  >> Number of peaks in the output table
% ---------
% For image reconstruction from cell list
% 'Integ' >> Integration range of the reconstructed plane. Default: 10 (um)
% 'StepD' >> Off-plane translational resolution. Default: 5 (um)
% 'StepX' >> In-plane translational resolution. Smaller resolution leads to
% better result but also with more computation time. Default: 7 (um)
% 'CellRadius' >> The radius for each neuron, used in the render process.
% Default: 3 (um)
% ---------
% For cell detection in the slice data
% 'Sigma'   >> [sx,sy],Size of the gaussian filter. Default: [8,8] (um)
% 'Res0'    >> Intensity threshold for the cell. Change with the resolution
% of the input dataS. Trial with neuroReg.bwCell2 is recommended.
% Default: 0.03
% 'SizeLimit' >> Size limit of the cell. Default: [100,2000] (um2)
% ---------
% For correlation function calculation
% MagicNumber >> I = I - MageicNumber * Mean(I(:)) is calculated before
% going to correlation function Default:2
% PeakThreshold >> from 0 to 1. The peaks larger than
% PeakThreshold * PeakIntensity will be regarded as a possible match.
% Defualt: 0.8
% ---------
% By Han Peng, penghan1992@gmail.com
% 2017.07.18
%
%
% **********
% * TEST *
% **********
% pt_list = evalin('base','pt_list2');
% data = evalin('base','dataS_resam_gen');
% N = 21;                             % N: Number of angles
% alpha = linspace(-5,5,N);           % alpha: rotation around y-axis
% beta = linspace(-5,5,N);            % beta: rotation around z-axis
% gamma = linspace(-5,5,N);           % gamma: rotation around x-axis
% Integ = 10;                         % Integ: slice thickness
% stepD = 5;                          % stepD: slice step (recommend: Integ/2)
% stepX = 7;                          % stepX: resolution of the slice
% CellRadius = 3;
% MagicNumber = 2;                    % Intensity ratio for cell and background



%% Rotation Angles
fprintf('neuroReg.rotationCorr3 running...\n')
alpha1 = AngleRange(1,1);
alpha2 = AngleRange(1,2);
n_alpha = AngleRange(1,3);
beta1 = AngleRange(2,1);
beta2 = AngleRange(2,2);
n_beta = AngleRange(2,3);
gamma1 = AngleRange(3,1);
gamma2 = AngleRange(3,2);
n_gamma = AngleRange(3,3);
alpha = linspace(alpha1,alpha2,n_alpha);           % alpha: rotation around y-axis
beta = linspace(beta1,beta2,n_beta);            % beta: rotation around z-axis
gamma = linspace(gamma1,gamma2,n_gamma);           % gamma: rotation around x-axis

%% Option
% Peak Records
Option = neuroReg.setOption(Option); % Check input, set default value.
TRANSLATION_TOLERANCE2 = Option.TransTol^2;
ANGLE_TOLERANCE2 = Option.AngleTol^2;
MAX_PEAK_RECORD_NUM = Option.MaxPeakNum;
DepthRange = Option.DepthRange;
% For image reconstruction from cell list
stepD = Option.StepD; % stepD: slice step (recommend: Integ/2)
stepX = Option.StepX; % stepX: resolution of the slice
MagicNumber = Option.MagicNumber; % Intensity ratio for cell and background

%% reconstruct data_slice using point list only
x_lim_slice = [data_slice.x(1),data_slice.x(end)];
y_lim_slice = [data_slice.y(1),data_slice.y(end)];
data_slice_bw_low = neuroReg.renderCell2(x_lim_slice,y_lim_slice,stepX,pt_list_slice,pt_area_slice,Option);
nf = isnan(data_slice_bw_low.value);
if Option.MismatchPen == 0
    data_slice_bw_low.value = data_slice_bw_low.value - mean(data_slice_bw_low.value(~nf))*MagicNumber;
end
data_slice_bw_low.value(nf) = 0;
figure(2)
imshow(data_slice_bw_low.value,[])

figure(100860); cla;
subplot(2,1,1)
neuroReg.plotData2(data_slice);
subplot(2,1,2)
neuroReg.plotData2(data_slice_bw_low);
title('Render result from the input slice');

%% Rotation and Calculation of Correlation Function
h1 = waitbar(0,'Calculating correlation function... Have a coffee...',...
    'Name','Calculating XCC... ');
data_corr.x = data_slice_bw_low.x;
data_corr.y = data_slice_bw_low.y;
N = n_gamma * n_beta * n_alpha;
% Record the peak positions
TransTable = table;
TransTable.Intensity = zeros(0);
TransTable.Angles = zeros(0,3);
TransTable.Translation = zeros(0,3);
tic;
for i = 1:n_gamma
    for j = 1:n_beta
        for k = 1:n_alpha
            %% Rotate the cells from the volume to slice coordination
            [pt_list_rotated,~,x_lim,y_lim,d_lim] = neuroReg.rotateCells(pt_list_vol,alpha(k),beta(j),gamma(i));
            % x_lim: dim1, y_lim: dim3, d_lim: dim2
            pt_list_rotated_now = [pt_list_rotated(1,:);pt_list_rotated(3,:);pt_list_rotated(2,:)];
            data_now = neuroReg.renderCell3(x_lim,y_lim,d_lim,stepX,stepD,pt_list_rotated_now,Option);
            % Get the planes to calculate correlation function
            x_c(1,1) = mean(x_lim); % center position (x) of the 2d plane
            x_c(3,1) = mean(y_lim); % center position (z) of the 2d plane
            %% Calculate Cross correlation by FFT
            % x-z: in-plane directions. y: normal direction of the plane.
            Im1 = data_slice_bw_low.value;
            ImZ = data_now.value; % 3D
            value_temp = neuroReg.xcorrFFT(Im1,ImZ);
            %% After the distance loop, record the maxima position
            % y_c: coordination after transformation (Slice coordination)
            % x_c: coordiniation before transformation (Volume coordination)
            % M: the transform from Volume coordinate system to Slice
            % coordinate sytem.
            value_temp_bw = value_temp; 
            % ignore pixels outside of bounding box defined by X Y
            BB = ones(size(value_temp_bw));
            BB (round(X(1)/stepX):round(X(2)/stepX),round(Y(1)/stepX):round(Y(3)/stepX),:)=0;
            value_temp_bw(logical(BB)) = 0;
            value_temp_bw(value_temp<0.1*max(value_temp(:))) = 0;
            
            CC = bwconncomp(value_temp_bw);
            PeakNum = CC.NumObjects;
            Angles = [alpha(k),beta(j),gamma(i)];
            for i_pk = 1:PeakNum
                % Get the Maximum intensity from each object
                [Intensity,i_pos] = max(value_temp_bw(CC.PixelIdxList{i_pk}));
                % Get the Maximum position from each object
                [y_c_temp_x,y_c_temp_z,y_c_temp_d] = ind2sub(CC.ImageSize,CC.PixelIdxList{i_pk}(i_pos));
                x_c(2,1) = data_now.z(y_c_temp_d);
                value_temp = value_temp - min(min(value_temp));
                % keep peak if in the region of interest
                if inpolygon(y_c_temp_x * stepX, y_c_temp_z*stepX, X,Y)
                    y_c = [y_c_temp_x;0;y_c_temp_z];
                    y_c(1) = y_c(1)*stepX+data_slice_bw_low.x(1)-stepX;
                    y_c(3) = y_c(3)*stepX+data_slice_bw_low.y(1)-stepX;
                    % Get the translation vector (after rotation. from volune
                    % to slice)
                    Translation = (-x_c+y_c)';
                    Angles_diff = TransTable.Angles - Angles;
                    af = (sum(Angles_diff.^2,2)<ANGLE_TOLERANCE2);
                    Translation_diff = TransTable.Translation(af,:) - Translation;
                    tf = (sum(Translation_diff.^2,2) < TRANSLATION_TOLERANCE2);
                    TransTableLength = size(TransTable,1);
                    % Modify the TransTable
                    if sum(tf)==0 %&& x_c(2)<350 && x_c(2)>250
                        % If this peak does not exist, then add it to the table
                        TransTable(TransTableLength+1,:) = {Intensity,Angles,Translation};
                    else%if x_c(2)<350 && x_c(2)>250
                        % If this peak already exists, then compare the intensities
                        index_list = 1:TransTableLength;
                        index_tf = index_list(af);
                        index_tf = index_tf(tf);
                        Intensity_diff = Intensity-TransTable.Intensity(index_tf);
                        intensf = (Intensity_diff>0);
                        index_change = index_tf(intensf);
                        % If the intensity is larger, then change it
                        if ~isempty(index_change)
                            TransTable(index_change(1),:) = {Intensity,Angles,Translation};
                            if length(index_change)>1
                                % delete duplicated peaks
                                TransTable(index_change(2:end),:) = [];
                            end
                        end
                    end
                    % Peaks.Corr, Peaks.M
                else
                    continue
                end
            end
            % Keep peaks in the depth range
            InDepth = - TransTable.Translation(:,2)> DepthRange(1) & - TransTable.Translation(:,2)<DepthRange(2);
            TransTable(~InDepth,:)=[];
            % Rank the Peaks with intensity
            TransTable = sortrows(TransTable,'Intensity','descend');
            TransTableLength = size(TransTable,1);
            if TransTableLength>MAX_PEAK_RECORD_NUM
                % Only keep the peaks with large intensities.
                TransTable((MAX_PEAK_RECORD_NUM+1):end,:)=[];
            end
            %% Visualization after each cube-XCC
            [corr_now_max,ind_max] = max(value_temp(:));
            [~,~,y_c_temp_d] = ind2sub(size(value_temp),ind_max);
            x_c(2,1) = data_now.z(y_c_temp_d); % x_c: center position for the volume for XCC maxima
            % Visualization for test
            data_corr.value(:,:) = value_temp(:,:,y_c_temp_d);
            data_this_cut = neuroReg.renderCell3(x_lim,y_lim,[x_c(2,1),x_c(2,1)],stepX,stepD,pt_list_rotated_now,Option);
            % Visualize current result
            figure(100861);
            subplot(1,2,1);
            neuroReg.plotData2(data_this_cut);
            subplot(1,2,2);
            neuroReg.plotData2(data_corr);
            title(['d=',num2str(data_now.z(y_c_temp_d)),' Max\_XCC=',num2str(corr_now_max)]);
            %% Waitbar
            i_tot = sub2ind([n_alpha,n_beta,n_gamma],k,j,i);
            a = toc;
            waitbar(i_tot/N,h1,[...
                'Current angle: ',num2str([alpha(k),beta(j),gamma(i)]),...
                '  Time remaining: ',datestr(a/i_tot*(N-i_tot)/24/3600,'DD HH:MM:SS')]);
        end % end of alpha cycle
    end % end of beta cycle
end % end of gamma cycle
close(h1);


end