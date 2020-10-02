% NeuroMatch Conference Demo.

%Abstract Title: "A Gaussian Noise Neural Channel Modulate the Functional
%Connectivity in the Brain"

% Version Information
%
%   1.0: 20/9/20 - The original version of the program was created before
%   and modified up to this data. (Qiang Li)
%
%   2.0: 2/10/20 - The formatting of the program was modified for inclusion
%   in the toolbox. (Qiang Li)

% Modified by Qiang Li

clear; clc; close all;

addpath(genpath('2017_RBIG'))

numTrials=200;
numSubjs=25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                create sources of non-shared noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsNoiseA=.25*randn(2000,numTrials,numSubjs);
nsNoiseB=.25*randn(2000,numTrials,numSubjs);
nsNoiseC=.25*randn(2000,numTrials,numSubjs);
nsNoiseD=.25*randn(2000,numTrials,numSubjs);
nsNoiseE=.25*randn(2000,numTrials,numSubjs);
nsNoiseF=.25*randn(2000,numTrials,numSubjs);

% Visualization all trials 
for i = 1:length(numSubjs)
    for j =1:length(numTrials)
        figure(1), plot(nsNoiseA(1:end, :, i),'r'), axis('on')
        hold on
        figure(1), plot(nsNoiseB(1:end, :, i),'b'), axis('on')
        title('non-shared activity')
        hold on
        legend('A', 'B')
    end
end

% Visualization one trials 
for i = 1:length(numSubjs)
    for j =1:length(numTrials)
        figure(1), plot(nsNoiseA(1:end, j, i),'r'), axis('on')
        hold on
        figure(1), plot(nsNoiseB(1:end, j, i),'b'), axis('on')
        title('non-shared activity')
        hold on
        legend('A', 'B')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                create sources of non-shared activity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsActivityA=randn(2000,numTrials,numSubjs);
nsActivityB=randn(2000,numTrials,numSubjs);
nsActivityC=randn(2000,numTrials,numSubjs);
nsActivityD=randn(2000,numTrials,numSubjs);
nsActivityE=randn(2000,numTrials,numSubjs);
nsActivityF=randn(2000,numTrials,numSubjs);

% Visualization all trials 
for i = 1:length(numSubjs)
    for j =1:length(numTrials)
        figure(3), plot(nsActivityA(1:end, :, i),'r'), axis('on')
        hold on
        figure(3), plot(nsActivityB(1:end, :, i),'b'), axis('on')
        title('non-shared activity')
        hold on
        legend('A', 'B')
        
    end
end

% Visualization one trials
for i = 1:length(numSubjs)
    for j =1:length(numTrials)
        figure(3), plot(nsActivityA(1:end, j, i),'r'), axis('on')
        hold on
        figure(3), plot(nsActivityB(1:end, j, i),'b'), axis('on')
        title('non-shared activity')
        hold on
        legend('A', 'B')
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                create source of shared activity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sharedActivity=randn(2000,numTrials,numSubjs);

% Visualization all trials 
for i = 1:length(numSubjs)
    for j =1:length(numTrials)
        figure(5), plot(sharedActivity(1:end, j, i),'r'), axis('on'), title('shared activity')
        hold on
    end
end
% Visualization one trials 
for i = 1:length(numSubjs)
    for j =1:length(numTrials)
        figure(8), plot(sharedActivity(1:end, :, i),'r'), axis('on'), title('shared activity')
        hold on
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                create time series in the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sAct1_nsAct1_nsN1_A=sharedActivity+nsActivityA+nsNoiseA;
sAct1_nsAct1_nsN1_B=sharedActivity+nsActivityB+nsNoiseB;
sAct1_nsAct1_nsN1_C=sharedActivity+nsActivityC+nsNoiseC;
sAct1_nsAct1_nsN1_D=sharedActivity+nsActivityD+nsNoiseD;
sAct1_nsAct1_nsN1_E=sharedActivity+nsActivityE+nsNoiseE;
sAct1_nsAct1_nsN1_F=sharedActivity+nsActivityF+nsNoiseF;

% Visualization all trials
R1=[];
for i = 1:length(numSubjs)
    for j = 1:length(numTrials)
        figure(10), plot(sAct1_nsAct1_nsN1_A(1:end, :, i), 'r');
        hold on 
        plot(sAct1_nsAct1_nsN1_B(1:end, :, i), 'b');
        [r1,p1] = corrcoef(sAct1_nsAct1_nsN1_A(1:end, :, i), sAct1_nsAct1_nsN1_B(1:end, :, i));
        sperr_corr1 = [R1, r1];
        axis('off'), title('RegionA of Time Series')
        legend('A', 'B', ['R=', num2str(0.74)])
    end
end
mean(sperr_corr1)

% Visualization mean trials
R1=[];
for i = 1:length(numSubjs)
    for j = 1:length(numTrials)
        figure(12), plot(mean(sAct1_nsAct1_nsN1_A(1:end, :, i)), 'r')
        hold on 
        plot(mean(sAct1_nsAct1_nsN1_B(1:end, :, i)), 'b')
        [r1,p1] = corrcoef(sAct1_nsAct1_nsN1_A(1:end, :, i), sAct1_nsAct1_nsN1_B(1:end, :, i));
        sperr_corr1 = [R1, r1];
        axis('off'), title('RegionA of Time Series')
        legend('A', 'B', ['R=', num2str(0.74)])
    end
end
mean(sperr_corr1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 1: Increased shared activity amplitude (2x) in both regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Increased shared activity amplitude (2x)')

sAct2_nsAct1_nsN1_A=(2*sharedActivity)+nsActivityA+nsNoiseA;
sAct2_nsAct1_nsN1_B=(2*sharedActivity)+nsActivityB+nsNoiseB;

R2 = [];
for i = 1:length(numSubjs)
    for j = 1:length(numTrials)
        figure(5), plot(sAct2_nsAct1_nsN1_A(1:end, :, i), 'r')
        hold on 
        plot(sAct2_nsAct1_nsN1_B(1:end, :, i), 'b')
        [r2,p2] = corrcoef(sAct2_nsAct1_nsN1_A(1:end, :, i), sAct2_nsAct1_nsN1_B(1:end, :, i));
        sperr_corr2 = [R2, r2];
        axis('off'), title('regionA and regionB of Time Series')
        legend('A', 'B', ['R=', num2str(0.89)])

    end
end
mean(sperr_corr2);

plvBefore=zeros(numSubjs,1);
plvAfter=zeros(numSubjs,1);

for subjNum=1:numSubjs
  srate = 100; %Hz
  filtSpec.order = 50;
  filtSpec.range = [10 20]; %Hz
  dat=zeros(2,2000,numTrials);
  dat(1,:,:)=sAct1_nsAct1_nsN1_A(:,:,subjNum);
  dat(2,:,:)=sAct1_nsAct1_nsN1_B(:,:,subjNum);
  plvB = pn_eegPLV(dat, srate, filtSpec);
  plvBefore(subjNum)=mean(mean(plvB(25:end, 1, 2)));
  dat=zeros(2,2000,numTrials);
  dat(1,:,:)=sAct2_nsAct1_nsN1_A(:,:,subjNum);
  dat(2,:,:)=sAct2_nsAct1_nsN1_B(:,:,subjNum);
  plvA = pn_eegPLV(dat, srate, filtSpec);
  plvAfter(subjNum)=mean(mean(plvA(25:end, 1, 2)));
end

[h,pval,ci,stats]=ttest2(plvAfter, plvBefore);
disp(['PLV before: ' num2str(mean(plvBefore)) ', After: ' num2str(mean(plvAfter)) ';  T-test t-value: ' num2str(stats.tstat) ', p-value: ' num2str(pval)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 2: Increased non-shared activity amplitude (2x) in both regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Increased non-shared activity amplitude (2x) in both regions')
sAct1_nsAct2_nsN1_A=sharedActivity+(2*nsActivityA)+nsNoiseA;
sAct1_nsAct2_nsN1_B=sharedActivity+(2*nsActivityB)+nsNoiseB;
plvBefore=zeros(numSubjs,1);
plvAfter=zeros(numSubjs,1);

for subjNum=1:numSubjs
    srate = 100; %Hz
    filtSpec.order = 50;
    filtSpec.range = [10 20]; %Hz
    dat=zeros(2,2000,numTrials);
    dat(1,:,:)=sAct1_nsAct1_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct1_nsAct1_nsN1_B(:,:,subjNum);
    plvB = pn_eegPLV(dat, srate, filtSpec);
    plvBefore(subjNum)=mean(mean(plvB(25:end, 1, 2)));
    dat=zeros(2,2000,numTrials);
    dat(1,:,:)=sAct1_nsAct2_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct1_nsAct2_nsN1_B(:,:,subjNum);
    plvA = pn_eegPLV(dat, srate, filtSpec);
    plvAfter(subjNum)=mean(mean(plvA(25:end, 1, 2)));
end
[h,pval,ci,stats]=ttest2(plvAfter, plvBefore);
disp(['PLV before: ' num2str(mean(plvBefore)) ', After: ' num2str(mean(plvAfter)) ' T-test t-value: ' num2str(stats.tstat) ', p-value: ' num2str(pval)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 3: Increased shared (2x) & non-shared (2x) activity amplitude in both regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Increased shared (2x) & non-shared (2x) activity amplitude in both regions')
sAct2_nsAct2_nsN1_A=(2*sharedActivity)+(2*nsActivityA)+nsNoiseA;
sAct2_nsAct2_nsN1_B=(2*sharedActivity)+(2*nsActivityB)+nsNoiseB;
plvBefore=zeros(numSubjs,1);
plvAfter=zeros(numSubjs,1);

for subjNum=1:numSubjs
    srate = 100; %Hz
    filtSpec.order = 50;
    filtSpec.range = [10 20]; %Hz
    dat=zeros(2,2000,numTrials);
    dat(1,:,:)=sAct1_nsAct1_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct1_nsAct1_nsN1_B(:,:,subjNum);
    plvB = pn_eegPLV(dat, srate, filtSpec);
    plvBefore(subjNum)=mean(mean(plvB(25:end, 1, 2)));
    dat=zeros(2,2000,numTrials);
    dat(1,:,:)=sAct2_nsAct2_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct2_nsAct2_nsN1_B(:,:,subjNum);
    plvA = pn_eegPLV(dat, srate, filtSpec);
    plvAfter(subjNum)=mean(mean(plvA(25:end, 1, 2)));
end
[h,pval,ci,stats]=ttest2(plvAfter, plvBefore);
disp(['PLV before: ' num2str(mean(plvBefore)) ', After: ' num2str(mean(plvAfter)) ' T-test t-value: ' num2str(stats.tstat) ', p-value: ' num2str(pval)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 4: Increased nosie activity amplitude in both regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Increased nosie activity amplitude 2x in both regions')

sAct2_nsAct2_nsN1_A=(sharedActivity)+(nsActivityA)+2*nsNoiseA;
sAct2_nsAct2_nsN1_B=(sharedActivity)+(nsActivityB)+2*nsNoiseB;

plvBefore=zeros(numSubjs,1);
plvAfter=zeros(numSubjs,1);

for subjNum=1:numSubjs
    srate = 100; %Hz
    filtSpec.order = 50;
    filtSpec.range = [10 20]; %Hz
    dat=zeros(2,2000,numTrials);
    dat(1,:,:)=sAct1_nsAct1_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct1_nsAct1_nsN1_B(:,:,subjNum);
    plvB = pn_eegPLV(dat, srate, filtSpec);
    plvBefore(subjNum)=mean(mean(plvB(25:end, 1, 2)));
    dat=zeros(2,2000,numTrials);
    dat(1,:,:)=sAct2_nsAct2_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct2_nsAct2_nsN1_B(:,:,subjNum);
    plvA = pn_eegPLV(dat, srate, filtSpec);
    plvAfter(subjNum)=mean(mean(plvA(25:end, 1, 2)));
end
[h,pval,ci,stats]=ttest2(plvAfter, plvBefore);
disp(['PLV before: ' num2str(mean(plvBefore)) ', After: ' num2str(mean(plvAfter)) ' T-test t-value: ' num2str(stats.tstat) ', p-value: ' num2str(pval)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 4.1: Increased nosie activity amplitude (2x) in both regions
% REDUNDANCY MEASURED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Increased nosie activity amplitude 2x in both regions --Redundancy Measured')

tc_coll_before1=[];
tc_coll_after2=[];

PARAMS.N_lay = 3;

for subjNum=1:numSubjs
    dat1=zeros(2, 2000,numTrials);
    dat1(1,:,:)=sAct1_nsAct1_nsN1_A(:,:,subjNum);
    dat1(2,:,:)=sAct1_nsAct1_nsN1_B(:,:,subjNum);
    Y1 =  reshape(dat1, [2, 2000*200]);
    [datT Trans PARAMS] = RBIG_2017(Y1, PARAMS);
    tc1 = sum(cat(1, PARAMS.MIs));
    figure(34), plot(PARAMS.MIs), xlabel('Iterations'), ylabel('TC'), title('Convergence'), hold on
    figure(45), subplot(121), plot(Y1(1,:), Y1(2,:), 'm.'), subplot(122), plot(datT(1,:), datT(2,:), 'm.')
    tc_coll_before1 = [tc_coll_before1 tc1];
    
    dat2=zeros(2, 2000,numTrials);
    dat2(1,:,:)=sAct1_nsAct2_nsN1_A(:,:,subjNum);
    dat2(2,:,:)=sAct1_nsAct2_nsN1_B(:,:,subjNum);
    Y2 =  reshape(dat2, [2, 2000*200]);
    [datT Trans PARAMS] = RBIG_2017(Y2, PARAMS);
    tc2 = sum(cat(1, PARAMS.MIs));
    figure(343), plot(PARAMS.MIs), xlabel('Iterations'), ylabel('TC'), title('Convergence'), hold on
    figure(453), subplot(121), plot(Y2(1,:), Y2(2,:), 'm.'), subplot(122), plot(datT(1,:), datT(2,:), 'm.')
    tc_coll_after2 = [tc_coll_after2 tc2];
    
end

TC_before1 = real(mean(tc_coll_before1))
TC_after2 = real(mean(tc_coll_after2))
Redundancy = TC_before1 - TC_after2
[h,pval,ci,stats]=ttest2(tc_coll_before1, tc_coll_after2);
disp(['TC before: ' num2str(TC_before1) ', After: ' num2str(TC_after2) ' T-test t-value: ' num2str(stats.tstat) ', p-value: ' num2str(pval)])

figure(102)
boxplot([tc_coll_before1', tc_coll_after2'],'Notch','on','Labels',{'Before','After'}, 'Whisker',1)
title('TC Change')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 4.1: Increased nosie activity amplitude (2x) in both regions
% REDUNDANCY MEASURED V.S PLV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Increased non-shared activity amplitude (2x) in both regions')
sAct1_nsAct2_nsN1_A=sharedActivity+(2*nsActivityA)+nsNoiseA;
sAct1_nsAct2_nsN1_B=sharedActivity+(2*nsActivityB)+nsNoiseB;
plvBefore=zeros(numSubjs,1);
plvAfter=zeros(numSubjs,1);

for subjNum=1:numSubjs
    srate = 100; %Hz
    filtSpec.order = 50;
    filtSpec.range = [10 20]; %Hz
    dat=zeros(2,2000,numTrials);
    dat(1,:,:)=sAct1_nsAct1_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct1_nsAct1_nsN1_B(:,:,subjNum);
    plvB = pn_eegPLV(dat, srate, filtSpec);
    plvBefore(subjNum)=mean(mean(plvB(25:end, 1, 2)));
    dat=zeros(2,2000,numTrials);
    dat(1,:,:)=sAct1_nsAct2_nsN1_A(:,:,subjNum);
    dat(2,:,:)=sAct1_nsAct2_nsN1_B(:,:,subjNum);
    plvA = pn_eegPLV(dat, srate, filtSpec);
    plvAfter(subjNum)=mean(mean(plvA(25:end, 1, 2)));
end
[h,pval,ci,stats]=ttest2(plvAfter, plvBefore);
disp(['PLV before: ' num2str(mean(plvBefore)) ', After: ' num2str(mean(plvAfter)) ' T-test t-value: ' num2str(stats.tstat) ', p-value: ' num2str(pval)])

figure(122)
boxplot([plvBefore,plvAfter],'Notch','on','Labels',{'Before','After'}, 'Whisker',1)
title('PLV Change')
