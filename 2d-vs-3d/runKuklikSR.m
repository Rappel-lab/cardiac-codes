
% Feb 21, 2017
% Copyright David Vidmar (UCSD), Mahmood Alhusseini (Stanford U), and
% Wouter-Jan Rappel (UCSD)

clear all
clc
close all

side_size=64;

%Load Data
offnum = 1; %0 for A, 1 for B, 2 for C, 3 for D, 4 for E, 5 for F, 6 for G, 7 for H
[FileName,PathName] = uigetfile('*.txt','Select the Basket data file');

if (length(FileName)<=1)
    return;
end
FullPath=[PathName FileName];
%Ask for file specific input
prompt = {'Enter Start Time (ms):','How many windows (2000 ms/window, default = 2):','Sampling Frequency ','Enter Offset:','Left or Right Atrium:','Invert Traces (1 = no, 2 = yes)?'};
dlg_title = 'Raw Data Information';
num_lines = 1;
if strcmp(FileName(1:2),'S-')
    defaultans = {'1','2','977','A','Left','1'};
else
    defaultans = {'1','2','1000','A','Left','1'};
end
Range = inputdlg(prompt,dlg_title,num_lines,defaultans);
Fs = str2double(Range{3});
StartTime = round((Fs/1000)*str2double(Range{1}));
EndTime = min([StartTime + 2*Fs*str2double(Range{2}) - 1, 4000]);
%set offset of this data
offnum = double(upper(Range{4})) - 65; %0 for A, 1 for B, 2 for C, 3 for D, 4 for E, 5 for F, 6 for G, 7 for H
LeftOrRight = upper(Range{5});
InvertElectrodes = str2double(Range{6});

%read in 64xN QRS subtracted data (A1,A2,A3,...,H7,H8)
Lecg = dlmread(FullPath); 
Lecg = Lecg(:,StartTime:EndTime);
if InvertElectrodes == 2
    Lecg = -Lecg;
end

%Get Size and correct matrix
[myR myC] = size(Lecg);
%Lecg = Lecg(:,StartTime:EndTime);
Lecg = circshift(Lecg,-8*offnum);  %fix offset


%SINE RECOMPOSITION CODE --------------------------------------------------
%//////////////////////////////////////////////////////////////////////////
%--------------------------------------------------------------------------

%run all pre-processing of signal
siglen = size(Lecg,2) ; sig = zeros(64,siglen); rsig = zeros(64,siglen);
processSig;

%compute sine recomposition per Kuklik et al. 2015 paper
for i = 1:64
 for t = 1:siglen-1
    if ~isnan(pd(i)) %check that there is some data in this channel
         halfpd = round(pd(i)/2);
         dvdt = diffv(i,t); tt = (-halfpd:halfpd);
         swave = sin(2*pi*tt/pd(i))*abs(dvdt)*((1-sign(dvdt))/2);
         twave = t + (-halfpd:halfpd); swave(twave<=0 | twave>siglen) = []; 
         twave(twave<=0)=[]; twave(twave>siglen)=[];
         rsig(i,twave) = rsig(i,twave) + swave;
    else
        rsig(i,:)=nan; %no data, mark to be interpolated over
    end
 end
end

%interpolate signals to 65x65 grid
rsig_int = nan(side_size,side_size,size(rsig,2)); Lamap = reshape(1:64,[8,8]);
[X,Y] = meshgrid(linspace(1,side_size,8),linspace(1,side_size,8)); [Xq,Yq] = meshgrid(1:side_size,1:side_size); Xspl = 1:8;
for t = 1:size(rsig,2)
    %interpolate over saturated regions of signal
    for splnum = 1:8 %look through all splines for saturation
        splelectrodes = Lamap(:,splnum);
        rsig_spline=rsig(splelectrodes,t);
        if any(isnan(rsig_spline)) % if any electrodes on this spline are saturated, fill in
           goodinds = ~isnan(rsig_spline);
           rsig(splelectrodes,t) = ...
               interp1(Xspl(goodinds),rsig_spline(goodinds),Xspl,'linear','extrap');
        end
    end
    
    %interpolate full 8x8 signal to 65x65 grid
    sig8x8 = reshape(rsig(:,t),[8,8]);
    rsig_int(:,:,t)=interp2(X,Y,sig8x8,Xq,Yq,'linear');
end

%compute phase via Hilbert Transform
rsig_int = reshape(rsig_int,[side_size*side_size,siglen]);
hV = hilbert(rsig_int');
rphi = nan(side_size*side_size,size(hV,1));
for i = 1:side_size*side_size
rphi(i,:) = mod(angle(hV(:,i))+pi/2,2*pi); %shift pi/2 so phase inverts at activation
end
rphi = reshape(rphi,[side_size,side_size,size(rphi,2)]); 

    
%Animate Phase Videos
%default is play = 1, to stop, make plotfig = 0
% plotfig = 1; %0 = run analysis but don't animate results
% 
% while true
%     prompt = {'Play Movie? 1 for yes, 2 for no'};
%     dlg_title = 'Repeat Video';
%     num_lines = 1;
%     defaultans = {'1'};
%     Range = inputdlg(prompt,dlg_title,num_lines,defaultans);
%     ReplayVideo = str2double(Range{1});
%     if ReplayVideo ~= 1
%         break;
%     end
%     AnimatePhaseVids;
% end
patientnum='ID-00';
datapath='DATA/';
save(strcat(datapath,'rec64phase_', patientnum, '.mat'), 'rphi')
save(strcat(datapath,'phase_', patientnum, '.mat'), 'rsig')