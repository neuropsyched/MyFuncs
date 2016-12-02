function [ResponseBlockTime, ResponseTimes, ForceResponseBlock] = MotivationTaskResponseTrig(Force, BinSize, pre, post, sf)

MinISI = 2.5*sf;
threshold = 200;

%h = waitbar(0,'Calculating derivative...');

DelForce = zeros(size(Force));
nsteps = size(Force,1)-2*BinSize-1;
for i=(BinSize+1):(size(Force,1)-BinSize)
    DelForce(i,:) = (mean(Force(i:(i+BinSize)),1)-mean(Force((i-BinSize):i),1))/(BinSize/sf);
    %waitbar(i/nsteps,h);
end

%waitbar(0,h,'Imposing MinISI...');

trig = find(DelForce>threshold);
lasttrig = trig(1);
for i = 2:length(trig)
    if (trig(i)-lasttrig) < MinISI
        trig(i)=-1;
    else
        lasttrig = trig(i);
    end
    %waitbar(i/(length(trig)-1),h);
end
trig(trig==-1)=[];

%waitbar(0,h,'Generating Force matrix...');

for i = 1:length(trig)
    ForceResponseBlock(:,i) = Force((trig(i)-pre*sf):(trig(i)+post*sf));
    %waitbar(i/length(trig),h);
end

%close(h);

%figure; plot(Force); hold on; plot(DelForce, 'r');

ResponseBlockTime = -pre:(1/sf):post;
ResponseTimes = trig;