function [pt_list,pt_area] = detectCells3(data,Option)
% detectCells3 detects cell positions from a ZStack data.
% [pt_list,pt_area] = detectCells3(data,Option)
% ---------
% OUTPUT:
% pt_list: positions of the detected cells (3-by-N array)
% pt_area: volume of each detected cells (1-by-N array)
% ---------
% INPUT
% data: data.x, data.y, data.z, data.value (see doc neuroReg)
% Option:
%   Note: pass [] to the function to use default settings.
%   -----------------------------
%   Options and default settings:
%   Option.Detect3Mode = 'Red';
%       >> Specify the detection mode. Can be 'Red' or 'Green'
%   Option.Sigma = [1 1 1]*4;
%       >> Size for difference of Gaussian filter. Unit: um
%   Option.Res0 = 0.0035;
%       >> Threshold for difference of Gaussian filter
%   Option.SizeLimit = [100,4500];
%       >> Estimization of the cell volume. Unite: um^3
%   Option.Sensitivity = 1.3;
%       >> Used in Green Mode. The larger this number is, the less cells
%          will be detected.
%   Option.MedianFilterSize = [7,7,7];
%       >> Size for median filter
%   Option.NeighborSize = [30 30 50];
%       >> Neighbor size in adaptive thresholding. Green Mode only. Unit:
%       um.


%% Preset options
if ~isfield(Option,'Detect3Mode')
    % Specify the detection mode
    Option.Detect3Mode = 'Red';
end
if ~isfield(Option,'Sigma')
    % Size for difference of Gaussian filter. Unit: um
    Option.Sigma = [1 1 1]*4;
end
if ~isfield(Option,'Res0')
    % Threshold for difference of Gaussian filter
    Option.Res0 = 0.0035;
end
if ~isfield(Option,'SizeLimit')
    % Estimization of the cell volume. Unite: um^3
    Option.SizeLimit = [100,4500];
end
if ~isfield(Option,'Sensitivity')
    % Used in Green Mode. The larger this number is, the less cells will be
    % detected.
    Option.Sensitivity = 1.3;
end
if ~isfield(Option,'MedianFilterSize')
    % Size for median filter
    Option.MedianFilterSize = [7,7,7];
end
if ~isfield(Option,'NeighborSize')
    % Used in Green Mode. Neighbor size in adaptive thresholding.
    Option.NeighborSize = [50 50 50];
end
disp(Option);

%% Preprocessing
% Use difference of Gaussian detector
Res0 = Option.Res0;
SizeLimit = Option.SizeLimit;
Sigma = Option.Sigma;

stepX = mean(diff(data.x));
stepY = mean(diff(data.y));
stepZ = mean(diff(data.z));
SizeLimitPx = SizeLimit/stepX/stepY/stepZ;
v = data.value/ max(data.value(:));
NullIndx = v==0;
v(NullIndx) = mean(v(:))*ones(length(v(NullIndx)),1);

%% Gaussian
% um2vx
sigmaInd = [1 1 1];
sigmaInd(2) = Sigma(1)/stepX; % in imgaussfilt, sigmaInd(2) is for 1st dimension
sigmaInd(1) = Sigma(2)/stepY; % sigmaInd(1) is for 2nd dimension
sigmaInd(3) = Sigma(3)/stepZ;
tic
fprintf('Applying difference of gaussian detector...\n');
Ic = imgaussfilt3(v,sigmaInd);
Is = imgaussfilt3(v,sigmaInd*5);
res = Ic - Is;
fprintf('Difference of gaussian detector done...\n');
toc
% illumination correction
%  Res = adaptthresh3(res,Sensitivity);
%  im_bw = single((res-Res)>0);

%%
[pt_list, pt_area] = neuroReg.gaussfiltdetectCells3(res, Res0, data, SizeLimitPx, stepX, stepY, stepZ,Option);

end
