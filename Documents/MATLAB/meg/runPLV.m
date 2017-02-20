% function [PLV]=PLV(varargin)
%%
%% 1. import recordings
clear; close all; clc

% User inputs the matlab path (mp) for their computer:
mp = '/Users/markrichardson/Documents/MATLAB';
% User inputs the general path (gp) to Projects folder for their computer:
gp = '/Users/markrichardson/Documents/MATLAB/brainstorm_db/MEG_Connectivity';
% User inputs the save path (sp) for their computer:
sp = '/Users/markrichardson/Documents/Projects/MEG_Connectivity';
% ------------------------------------------------------
% ------------------------------------------------------
% % path should include brainstorm3 and brainstorm_db
% addpath(genpath([mp '/brainstorm3']))
% addpath(genpath([mp '/brainstorm_db']))
addpath([mp '/meg'])

% locate Subject folder
folder_name = uigetdir(gp);
PtId=folder_name(strfind(folder_name,'data/')+5:end);
folders=dir(folder_name);
folders=cat(2,{folders(:).name});
folders=folders(3:end);

% select PtId_Resting_raw
[s,~] = listdlg('PromptString','Select a recording:',...
    'SelectionMode','single',...
    'ListString',folders);
files=dir([folder_name '/' folders{s}]);
files=cat(2,{files(:).name});
% cd([folder_name '/' folders{s}])
try 
    if strfind(folders{s},'Resting'); recordtype = 'Resting'; 
    elseif strfind(folders{s},'Empty'); recordtype = 'EmptyRoom'; end
catch
    error; 
end
%% load Data from source model
idx = find(~cellfun(@isempty,strfind(files,'matrix')));
% if ~isempty
for i = 1:size(idx,2)
    load([folder_name '/' folders{s} '/' files{idx(1,i)}])
    if i == 1
        hippoValue = Value;
        hippoComment = Comment;
        hippoDescription = Description;
        clear Value Comment Description
    elseif i == 2
        cortexValue = Value;
        cortexComment = Comment;
        cortexDescription = Description;
        clear Value Comment Description
    end
end
% end
clear files folder_name folders i idx s


% vertically concatenate data
% rows 1:148 = cortex scouts 1:148 (see Description.cortex)
% rows 149:150 = hippocampus scouts 1 and 2
data = [cortexValue ; hippoValue];
% Commnet
Comment.cortex = cortexComment;
Comment.hippocampus = hippoComment;
% Description
Description = cortexDescription;
Description{end+1} = hippoDescription{1};
Description{end+1} = hippoDescription{2};
% remove Value and Time files from workspace to free up memory
clear hippoValue cortexValue hippoC* cortexC* hippoD* cortexD* Time

%% Find Phase Locking Values
% set wavelet parameters
srate = 200; %Hz
dsr = 4; % down sample by a factor of 4
bandfreq = [1 4;4 8;8 12;12 16;16 20;20 24;24 28; 28 32];
% bandfreq determines the bandpass filter upper and lower frequencies

% find other parameters
nrate = srate/dsr; % new sampling rate after downsampling by a factor of 4
n_sources=size(data,1);

% Outside loop to run through different bandpass frequencies
for i = 1:size(bandfreq,1)
    % bandpass filter and hilbert
    complex = hilbert(eegfilt(data,srate,bandfreq(i,1),bandfreq(i,2))');
    % divide complex by magnitude of complex and downsample by dsr
    complex = downsample(complex./abs(complex),dsr);
    % cutoff artifact at edges caused by bandpass filter
    complex = complex(10*nrate:size(complex,1)-190*nrate,:);
    % initiate PLV variable
    PLV = zeros(n_sources,n_sources);
    % Calculate PLV in parfor loop
    tic
    for j=1:n_sources
        PLV(j,:) = [zeros(1,j-1) abs(mean(bsxfun(@times,conj(complex(:,j)),complex(:,j:n_sources))))];
        if mod(j,500)==0
            disp([ 'channel ' num2str(j) '  out of ' num2str(n_sources)]);toc
        end
    end
    toc
    % gather PLV
    PLVr=PLV+PLV';
    eval(['PLV' num2str(bandfreq(i,1)) 'to' num2str(bandfreq(i,2)) ' = PLV;']);
    eval(['PLVr' num2str(bandfreq(i,1)) 'to' num2str(bandfreq(i,2)) ' = PLVr;']);
    clear PLV PLVr
end
clear i j

%% Plot Figures
% close all
% for i = 1:length(bandfreq)
%     figure; 
%     % pcolor(PLV);
%     eval(['pcolor(PLV' num2str(bandfreq(i,1)) 'to' num2str(bandfreq(i,2)) '); colorbar']);
%     eval(['title(''bandpass at ' num2str(bandfreq(i,1)) ' to ' num2str(bandfreq(i,2)) ' Hz'')']);
%     if i <5
%         set(gcf,'Position',[i*440-440 660 430 324])
%     elseif i >4
%         set(gcf,'Position',[(i-4)*440-440 250 430 324])
%     end
% %     figure;
% %     % PLVr = PLV+PLV';
% %     eval(['pcolor(PLVr' num2str(bandfreq(i,1)) 'to' num2str(bandfreq(i,2)) '); colorbar']);
% %     eval(['title(''bandpass at ' num2str(bandfreq(i,1)) ' to ' num2str(bandfreq(i,2)) ' Hz'')']);
% %     if i <5
% %         set(gcf,'Position',[i*460-460 640 430 324])
% %     elseif i >4
% %         set(gcf,'Position',[(i-4)*460-460 230 430 324])
% %     end
% end

%% Save Data and Write to File
% cd(fullfile(sp,PtId,'analysis'))
clear i j ans nAvg
save(fullfile(sp,PtId,'analysis',[PtId '_' recordtype '_05PLV_' date '.mat']))
PLVwritetofile
disp('DONE')