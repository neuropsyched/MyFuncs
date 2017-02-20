% function running_ICA()
% running_ICA
%% 1. import recordings
clear; close all; clc
% elapsed=0;
% User inputs the matlab path (mp) for their computer:
mp = '~/Documents/MATLAB';
% User inputs the general path (gp) for their computer:
% gp = '/Users/markrichardson/Documents/MATLAB/brainstorm_db/MEG_Connectivity';
gp = bst_get('BrainstormDbDir');
% User inputs the save path (sp) for their computer:
sp = '/Users/markrichardson/Documents/Projects/MEG_Connectivity';
% ------------------------------------------------------
% ------------------------------------------------------
% % path should include eeglab13_6_5b and eeglab13_6_5b/functions/octavefunc
% addpath(genpath([mp '/eeglab13_6_5b']))
% rmpath(genpath([mp '/eeglab13_6_5b/functions/octavefunc/']))
% addpath(genpath([mp '/brainstorm3']))
% addpath(genpath([mp '/brainstorm_db']))
% rmpath(genpath([mp '/brainstorm3/external/']))
% addpath(genpath([mp '/brainstorm3/external/mne']))
cd([mp '/meg'])
% locate Subject folder
% Example: KW081215
folder_name = uigetdir(gp);
PtId=folder_name(strfind(folder_name,'data/')+5:end);
folders=dir(folder_name);
folders=cat(2,{folders(:).name});
folders=folders(3:end);
% select PtId_Resting_raw, e.g. KW081215_Resting_raw
[s,~] = listdlg('PromptString','Select a recording:',...
    'SelectionMode','single',...
    'ListString',folders);
files=dir([folder_name '/' folders{s}]);
files=cat(2,{files(:).name});
% determine EmptyRoom v. Resting (recordtype)
% determine EmptyRoom v. Resting (recordtype)
try 
    if strfind(folders{s},'Resting'); recordtype = 'Resting'; 
    elseif strfind(folders{s},'Empty'); recordtype = 'EmptyRoom'; end
catch
    error; 
end
% Load and Save Original Data
% load '.../data/KW081215/KW081215/KW081215_Resting_raw/data_block001.mat'
load([ folder_name '/' folders{s} '/' files{find(~cellfun(@isempty,strfind(files,'data')))}],'F','Time','Events','ChannelFlag','History');
% load('.../data/KW081215/KW081215_Resting_raw/channel_vectorview306_acc1.mat','Channel')
load([ folder_name '/' folders{s} '/' files{find(~cellfun(@isempty,strfind(files,'channel')))}],'Channel');

% save original data
% locate original data 
tic
origstr=History{1,3};
origstr=origstr(strfind(origstr,': /')+2:strfind(origstr,'.fif')+3);
% save original data
orig = F;
vars = whos; % {vars(:).name}';
try
    save(fullfile(sp,PtId,'analysis',[PtId '_' recordtype '_01orig.mat']),vars(:).name);
catch
    mkdir(fullfile(sp,PtId,'analysis'))
    save(fullfile(sp,PtId,'analysis',[PtId '_' recordtype '_01orig.mat']),vars(:).name);
end
elapsed.save_orig = toc
% load /Users/markrichardson/Documents/Projects/MEG_Connectivity/CK053116/analysis/CK053116_Resting_orig.mat

%% 2. Remove Bad Channels &  High Pass Filter
% Identify bad channels from vector, ChannelFlag
% tic
badch=find(ChannelFlag~=1);
% Separate channels in to GRAD and MAG and REJECT bad channels
name=cat(1,{Channel(:).Name});
type=cat(1,{Channel(:).Type});
Mind=3:3:306; % index MAGnetometers
Gind=setdiff(1:306,Mind); % index GRADiometers

if ~isempty(badch) % if badch is not empty
    Gind(ismember(Gind,badch))=[];  % remove bad GRADiometers
    Mind(ismember(Mind,badch))=[];  % remove bad MAGnetometers
end
% Plot data
% eegplot(F,'srate',1000,'winlength',25)
% eegplot(F(Mind,:),'srate',1000,'winlength',25)
% eegplot(F(Gind,:),'srate',1000,'winlength',25)


% High Pass Filter
% Filter data
% A minimum-order highpass digital FIR filter with:
%       normalzied stopband frequency:  0.1pi rad/s
%       passband frequency:             1pi rad/s
%       passband ripple:                0.01 dB
%       stopband attenuation:           65 dB 
% Use the Parks-McClellan algorithm to design an equiripple FIR filter. 
% Equiripple filters have a frequency response that minimizes the maximum 
% ripple magnitude over all bands.
if ~exist('hpFilt','var')
    tic
    fprintf('Designing High Pass Filter...\n')
    hpFilt = designfilt('highpassfir',...   % Response type
    'PassbandFrequency',1, ...              % Frequency constraints
    'StopbandFrequency',0.1,...
    'StopbandAttenuation',65,...            % Magnitude constraints   
    'PassbandRipple',0.01, ...
    'DesignMethod','equiripple', ...        % Design method option
    'SampleRate',1000);                     % Sample rate
    %     elapsed = elapsed+toc;
    %     fprintf('Time Elapsed Designing High Pass Filter: %d\n',toc)
    %     fprintf('Total Time Elapsed: %d',elapsed)
    elapsed.hpFilt = toc
elseif exist('hpFilt','var') % if hpFilt was already loaded in
    %     fprintf('Designing High Pass Filter...')
    %     fprintf('Time Elapsed Designing High Pass Filter: %d\n',0)
    %     fprintf('Total Time Elapsed: %d\n',elapsed)
    elapsed.hpFilt = 0
end
% Filter the data forwards and backwards
tic
fprintf('Running High Pass Filter...\n')
F=filtfilt(hpFilt,F')'; %(setdiff(1:size(F,1),badch),:)')';
elapsed.filtfilt = toc
% elapsed = elapsed+toc;
% fprintf('Time Elapsed Running High Pass Filter: %d\n',0)
% fprintf('Total Time Elapsed: %d\n',elapsed)

% Remove Bad Segments
% inititate variables
if size(Events)>0
    %Time=round(Time,4,'decimal');
    %badseg=round(Events.times,4,'decimal');
    badseg=Events.times;
elseif isempty(Events)
    badseg = []
end

tic
fs=round(1/diff(Time(1:2))); %1000;
if fs~=1000
    error('check srate')
end

artifact=[];
for i = 1:size(badseg,2)
    artifact = [artifact,intersect(find(Time>badseg(1,i)==1),find(Time<badseg(2,i)==1))];    
end
F=F(:,setdiff(1:length(F),artifact));
Time=[0:1/fs:(length(F)-1)/fs]+Time(1);
elapsed.removebadsegs = toc

tic
clear ans orig i 
save(fullfile(sp,PtId,'analysis',[PtId '_' recordtype '_02removebadchfiltfilt.mat']))
elapsed.save_removebadchfiltfilt = toc
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% Downsample data
% keep sp PtId recordtype
% load(fullfile(sp,PtId,'analysis',[PtId '_' recordtype '_removebadchfiltfilt.mat']))

clear jump
jump(1,:) = find(diff(artifact)>1);
% jump(2,:) = [jump(2:end)-1,size(artifact,2)];
% jump = artifact(find(diff(artifact)>1));
if isempty(jump) && ~isempty(artifact)
    jump(1,1)  = 1;
%     jump(2,1)  = length(artifact);
end
if max(jump(1,:))>size(F,2)-30
    jump(1,end)=size(F,2)-30;
    sprintf('jump was adjusted because it was out of bounds of matrix F')
end
% Adjust jump
dAJump = artifact(jump)-sum(diff(jump)); %get index with artifact removed

tic
dF = zeros(size(F));
dF(Mind,:)=detrend(F(Mind,:)')';
dF(Gind,:)=detrend(F(Gind,:)')';
if ~isempty(jump)
    for i=1:size(dAJump,2)
    tmp=arrayfun(@(x) smooth(dF(x,dAJump(i)-30:dAJump(i)+30),...
        60,'rloess'),1:306,'UniformOutput',0);
    tmp=horzcat(tmp{:})';
    dF(1:306,dAJump(i)-30:dAJump(i)+30)=tmp;
    end
end
elapsed.smooth = toc
% create low pass FIR filter
% eegplot(dF(Mind,:),'srate',fs,'winlength',25,'plottitle','Magnetometers')
tic;
lpFilt = designfilt('lowpassfir','PassbandFrequency',45, ...
'StopbandFrequency',55,'PassbandRipple',0.01, ...
'StopbandAttenuation',65,'DesignMethod','equiripple','SampleRate',fs);
elapsed.lpfilt = toc

tic
d=5;
dfs=fs/d; %200;
dF=downsample(filtfilt(lpFilt,dF'),d)';
elapsed.Downsample = toc

%
% REJECT NEW 
% eegplot(dF(Mind,:),'srate',dfs,'winlength',25,'plottitle','Magnetometers','command',cmd)
% eegplot(dF(Gind,:),'srate',dfs,'winlength',25,'plottitle','Gradiometers')

% load('CK053116_Resting_03badsegs2.mat')
% save(fullfile(sp,PtId,'analysi s',[PtId '_' recordtype '_03badsegs2.mat']),'badsegs2')
% load(fullfile(sp,PtId,'analysis',[PtId '_' recordtype '_03badsegs2.mat']))
% dfs=200;
clear TMPREJ
dTime=downsample(Time,d);
cmd = ['[dFc, badsegs2] = eeg_rejbadsegs(TMPREJ,dF,dTime,dfs);'];
eegplot(dF(Mind,:),'srate',dfs,'winlength',25,'plottitle','Magnetometers','command',cmd)
hand = findobj('tag', 'EEGPLOT');
waitfor(hand)

while exist('TMPREJ','var')
    clear TMPREJ
    cmd = ['[dFc, badsegsplus] = eeg_rejbadsegs(TMPREJ,dFc,dTime,dfs);'];
    eegplot(dFc(Mind,:),'srate',dfs,'winlength',25,'plottitle','Magnetometers','command',cmd);
    hand = findobj('tag', 'EEGPLOT');
    waitfor(hand)
end
% badsegs2=sort(badsegs2,badsegsplus];

savetog = questdlg('savenow?');
if strcmp(savetog,'Yes')
    clear dF F
    save(fullfile(sp,PtId,'analysis',[PtId '_' recordtype '_03removebadchfiltfiltdownsampleclean.mat']))
    elapsed.save_removebadchfiltfiltdownsampleclean = toc
else
    return
end
% ,...
%         'PtId','dFc','Events','Mind','Gind','mp','sp','gp','Channel','Channelmat',...
%         'Time','History',...
%          'artifact','badch','badseg','badsegs2','TEMPREJ','dt','fs','dfs','origstr')
% clear dF F dt i jump artifact ans badseg badsegs2 TEMPREJ

%% Build the SSP

tic
% Use fieldtrip to read header
addpath(genpath([mp '/fieldtrip-20160627']))
    [header]=ft_read_header(origstr);
rmpath(genpath([mp '/fieldtrip-20160627']))

for i=1:length(header.orig.projs)
ChannelMat.Projector(i).Comment=header.orig.projs(i).desc;
ChannelMat.Projector(i).Components=cat(1,header.orig.projs(i).data.data',zeros(header.nChans-306,1));
ChannelMat.Projector(i).CompMask=1;
ChannelMat.Projector(i).Status=1;
ChannelMat.Projector(i).SingVal=[];
end
% addpath(genpath('/home/richardsonlab/Dropbox/Toolboxes/brainstorm3'))

Projector = process_ssp2('BuildProjector', ChannelMat.Projector, 1);
% rmpath(genpath('/home/richardsonlab/Dropbox/Toolboxes/brainstorm3'))

% Apply projector
% Remove bad channels from the projector 
if ~isempty(badch)
    Projector(badch,:) = 0;
    Projector(:,badch) = 0;
    Projector(badch,badch) = eye(length(badch));
end
% Apply projector

if exist('dFc','var')
    F = Projector * dFc;
    clear dFc
else
    F = Projector * dF;
end
elapsed.SSP = toc

% look at components
hand = figure
subplot(2,2,1:2)
[~, ~, ~, ~, EXPLAINEDM] = pca(F(Mind,:)');
plot(cumsum(EXPLAINEDM)); title('Magnetometers PCA')
hold on
line([0 length(EXPLAINEDM)],[95 95],'Color','r');
line([0 length(EXPLAINEDM)],[97.5 97.5],'Color','r','LineStyle','--');

subplot(2,2,3:4)
[~, ~, ~, ~, EXPLAINEDG] = pca(F(Gind,:)');
plot(cumsum(EXPLAINEDG));title('Gradiometers PCA')
line([0 length(EXPLAINEDG)],[95 95],'Color','r');
line([0 length(EXPLAINEDG)],[97.5 97.5],'Color','r','LineStyle','--');

% input number of ICA components 
waitfor(hand)
prompt = {'Enter number of ICAs for MAG','Enter number of ICAs for GRAD'};
dlg_title = 'ICA';
num_lines = 1;
def = {'25','100'}; % default iterations
explained = inputdlg(prompt,dlg_title,num_lines,def);
 
%% run ICA and save
% eegplot(F(Mind,:),'srate',dfs,'winlength',25,'plottitle','Magnetometers')
if strcmp(recordtype,'Resting')
    tic
    [weightsM,spheresM]=runica(F(Mind,:),'extended',1,'pca',str2num(explained{1}),'maxsteps',1e3);
    [weightsG,spheresG]=runica(F(Gind,:),'extended',1,'pca',str2num(explained{2}),'maxsteps',1e3);
    elapsed.ICA = toc
    %
    close all
    eegplot(weightsM*spheresM*F(Mind,:),'srate',dfs,'winlength',30,'dispchans',20,'plottitle','MagnetometersICA')
    hand = findobj('tag','EEGPLOT'); waitfor(hand)
    tmp = inputdlg('Enter space limited artifact channels');
    aChM = str2double(strsplit(tmp{1},' '));
    if ~isnan(aChM)
        [cleandataM]=ahmad_icaproj(weightsM,spheresM,F(Mind,:),aChM);
    end
    
    close all
    eegplot(weightsG*spheresG*F(Gind,:),'srate',dfs,'winlength',30,'dispchans',20,'plottitle','GradiometersICA')
    hand = findobj('tag','EEGPLOT'); waitfor(hand)
    tmp = inputdlg('Enter space limited artifact channels');
    aChG = str2double(strsplit(tmp{1},' '));
    if ~isnan(aChG)
        [cleandataG]=ahmad_icaproj(weightsG,spheresG,F(Gind,:),aChG);
    end
end
clear tmp t1 t2 i
save(fullfile(sp,PtId,'analysis',[PtId '_' recordtype '_04removebadchfiltfiltICAcleaned.mat']))

% Backproject
tic
if exist('cleandataM','var')
    F(Mind,:)=cleandataM;
end
if exist('cleandataG','var')
    F(Gind,:)=cleandataG;
end
Time = 0:1/dfs:(size(F,2)-1)/dfs;
Comment = sprintf('CleanICA (%0.2fs,%0.2fs)', Time(1),Time(end));
% History{3,3} = sprintf('    [%0.6f, %0.6f]', Time(1), Time(end));

% save([ folder_name '/' folders{s} '/' files{find(~cellfun(@isempty,strfind(files,'data')))}],'F', 'Time','ChannelFlag', 'Comment', 'DataType', 'Device', 'Events','History', 'nAvg')
save([ folder_name '/' folders{s} '/' files{find(~cellfun(@isempty,strfind(files,'data')))}],'F', 'Time','Comment','-append')
elapsed.Backproject = toc
% end

%%
% eeg_rejbadsegs command
% function [cleandata,badsegs] = eeg_rejbadsegs(TMPREJ,data,time,srate)
% 
% % SYNTAX:
% %   [dFc, badsegs2] = eeg_rejbadsegs(TMPREJ,F,dTime,dfs)
% 
% if isempty(TMPREJ(TMPREJ>1))
%     return
% else 
%     nsegs = size(TMPREJ,1);
% end
% 
% dt=[];
% cnt=1;
% for i = 1:nsegs
%     badsegs(1,cnt) = TMPREJ(cnt)/srate+time(1);
%     badsegs(2,cnt) = TMPREJ(cnt+1)/srate+time(1);
%     dt = [dt,intersect(find(time>badsegs(1,i)==1),find(time<badsegs(2,i)==1))];
%     cnt=cnt+1;
% end
% 
% cleandata=data(:,setdiff(1:length(data),dt));
% end
