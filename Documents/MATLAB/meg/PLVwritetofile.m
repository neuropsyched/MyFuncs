% %% Write PLVscript Output to Excel Workbook
% % Separate into 5 sheets
% % delta: 1-4Hz (PLVr1to4)
% % theta: 4-8Hz (PLVr4to8)
% % alpha: 8-12Hz (PLVr8to12)
% % low-beta: 12-20Hz (PLVr12to16/PLVr16to20)
% % hi-beta:  20-28Hz
% clear, clc
% % Input Computer's data path
% dp = '/Users/markrichardson/Documents/Projects/MEG_Connectivity/';
% % ------------------------------------------------------
% % ------------------------------------------------------
% % load data from GUI
% % addpath(genpath(dp))
% [file,path,ext] = ea_uigetfile(dp,'Choose PLV File');
% load([char(path),filesep,[char(file),char(ext)]])

%% Organize PLV outputs
delta = PLVr1to4;
theta = PLVr4to8;
alpha = PLVr8to12;
lowBeta = zeros(n_sources,n_sources);
hiBeta = zeros(n_sources,n_sources);
for i = 1:n_sources
    for j = 1:n_sources
        lowBeta(i,j) = mean([PLVr12to16(i,j),PLVr16to20(i,j)]);
        hiBeta(i,j) = mean([PLVr20to24(i,j),PLVr24to28(i,j),PLVr28to32(i,j)]);
    end
end

% keep dp delta theta alpha lowBeta hiBeta PtId
%% Write Data to Xcel File (warning: output is .csv)
% Create Save path (sp)
% spath = path{1};
spath = fullfile(sp,PtId,'analysis');
if ~exist([spath '/xcel'],'dir')
    mkdir([spath '/xcel'])
end
% delta: 1-4Hz
filename = [spath '/xcel/' PtId '_delta.xlsx'];
csvwrite(filename,delta); %,'delta');
% theta: 4-8Hz
filename = [spath '/xcel/' PtId '_theta.xlsx'];
csvwrite(filename,theta); %,'theta');
% alpha: 8-12Hz
filename = [spath '/xcel/' PtId '_alpha.xlsx'];
csvwrite(filename,alpha); %,'alpha');
% low-beta: 12-20Hz
filename = [spath '/xcel/' PtId '_lowBeta.xlsx'];
csvwrite(filename,lowBeta); %,'lowBeta');
% hi-beta:  20-28Hz
filename = [spath '/xcel/' PtId '_hiBeta.xlsx'];
csvwrite(filename,hiBeta); %,'hiBeta');
