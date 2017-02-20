function [cleandata,badsegs] = eeg_rejbadsegs(TMPREJ,data,time,srate)

% SYNTAX:
%   [dFc, badsegs2] = eeg_rejbadsegs(TMPREJ,F,dTime,dfs)

if isempty(TMPREJ(TMPREJ>1))
    return
else 
    nsegs = size(TMPREJ,1);
end

dt=[];
for i = 1:nsegs
    badsegs(1,i) = TMPREJ(i,1)/srate+time(1);
    badsegs(2,i) = TMPREJ(i,2)/srate+time(1);
    dt = [dt,intersect(find(time>badsegs(1,i)==1),find(time<badsegs(2,i)==1))];
end

cleandata=data(:,setdiff(1:length(data),dt));