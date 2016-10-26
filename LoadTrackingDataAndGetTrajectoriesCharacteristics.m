% load results (workspace) from DiaTrack 
clear all
[FileName,PathName] = uigetfile('*.mat','Select the DiaTrack output file');
load([PathName FileName]);
FrameTime=0.8; %frame time in seconds
mintrajsize=20; % the smallest (in frames) trajectory size to be considered later
[FileName1,PathName] = uigetfile('*.tif','Select the acceptor channel image');
dataApostThresh=rd_img16([PathName FileName1]);
[FileName2,PathName] = uigetfile('*.tif','Select the acceptor channel image');
dataDpostThresh=rd_img16([PathName FileName2]);
ratAtoD=double(dataApostThresh)./double(dataDpostThresh);
% extract all trajectories coordinates, intensity of features, time points
% at which trajectory appears and calcualte instantaneous diffusion
% coefficient...

%%
counterTrack=1;
for frame=firstFrameSub:lastFrameSub
for part=find(tracks{frame}.Predecessor==0)
framsav=frame;
partsav=part;
accumtrack(1,:)=[tracks{frame}.RefinedCooY(part),tracks{frame}.RefinedCooX(part),tracks{frame}.Intensity(part),frame];
while(tracks{frame}.Successor(part)&&frame<lastFrameSub)
part=tracks{frame}.Successor(part);
accumtrack(frame-framsav+2,:)=[tracks{frame+1}.RefinedCooY(part),tracks{frame+1}.RefinedCooX(part),tracks{frame+1}.Intensity(part),frame+1];
frame=frame+1;
end;

xtraj{counterTrack}=accumtrack(:,1);
ytraj{counterTrack}=accumtrack(:,2);
PartIntensity{counterTrack}=accumtrack(:,3);
timePart{counterTrack}=accumtrack(:,4);

v=accumtrack(1:frame-framsav+1,:);
if size(v,1)>1

vtot=zeros(1,size(v,1)-1);
for gap=1:(size(v,1)-1)
for posStart=1:(size(v,1)-gap)
vd=v(posStart+gap,:)-v(posStart,:);
vs=sum(vd.^2);
vtot(gap)=vtot(gap)+vs;
end;
vtot(gap)=vtot(gap)/(size(v,1)-gap);
end;
slopeavg=polyfit([0:size(v,1)-1],[0,vtot],1);
myTrajStat=slopeavg(1);
else
myTrajStat=Inf;
end;
diff(counterTrack)= myTrajStat;
 counterTrack=counterTrack+1;
frame=framsav;
 part=partsav;
accumtrack=[];
end;
end;


% set the minimum size trajectory to be considered (in frames) and extract
% the 'good' trajectories indexes...
%%
for i=1:length(xtraj)
 if length(xtraj{1,i})>mintrajsize
     goodtraj(i)=1;
 else
     goodtraj(i)=0;
 end
end

 indgoodtraj=find(goodtraj);
for i=1:length(indgoodtraj)
lengthGood(i)=length(timePart{indgoodtraj(i)});
end
%%%%
%%%%%%%%%%%% define the center of cell
%%
clear colormap
imagesc(mean(dataApostThresh+dataDpostThresh,3));colormap(gray)
title('Please select cell center or cSMAC by clicking on it','FontSize',16)
[ycell,xcell]=ginput(1);
close 

% for every 'good' trajectory find the average amplitude of the
% feature, at every time point, considering the 8 neighboring intensity of peak pixel.  
%%
set(gcbf,'pointer','watch');
h = waitbar(0,'Calculating trajectories stats...');
for j=1:length(indgoodtraj); 
times=timePart{indgoodtraj(j)};
intensities=PartIntensity{indgoodtraj(j)};
xinst=xtraj{indgoodtraj(j)};
yinst=ytraj{indgoodtraj(j)}; 
Xcm(j)=mean(xinst);
Ycm(j)=mean(yinst);
StartTraj(j,1)=xinst(1);
StartTraj(j,2)=yinst(1);
EndTraj(j,1)=xinst(length(xinst));
EndTraj(j,2)=yinst(length(yinst));
DistStartToCM(j)=sqrt((StartTraj(j,2)-ycell)^2+(StartTraj(j,1)-xcell)^2);
DistEndToCM(j)=sqrt((EndTraj(j,2)-ycell)^2+(EndTraj(j,1)-xcell)^2);
DeltaRtoCM(j)=DistStartToCM(j)-DistEndToCM(j);
DistTrajCMtoCellCM(j)=sqrt((Ycm(j)-ycell)^2+(Xcm(j)-xcell)^2);
[xc,yc,theta,minor,major,mp20,mp02,mp11]=trajectoryPropsByMoments(xinst,yinst);
MajorToMinorRatioTrajectory(j)=major/minor;

masks=zeros(size(FiltImg,1),size(FiltImg,2),length(times));
clear intfeat intfeatRplusG
for i=1:length(xinst)
    centx=round(yinst(i));
    centy=round(xinst(i));
    if centx>=size(FiltImg,2)
        centx=size(FiltImg,2)-1;
    elseif centx<=1
        centx=2;
    end
    if centy>=size(FiltImg,1)
        centy=size(FiltImg,1)-1;
    elseif centy<=1
        centy=2;
    end
    masks(centy,centx,i)=1;
    masks(centy-1,centx-1,i)=1;
    masks(centy-1,centx,i)=1;
    masks(centy,centx-1,i)=1;
    masks(centy+1,centx,i)=1;
    masks(centy,centx+1,i)=1;
    masks(centy+1,centx+1,i)=1;
    masks(centy+1,centx-1,i)=1;
    masks(centy-1,centx+1,i)=1;
    intfeat(i)=squeeze(nansum(nansum(masks(:,:,i).*double(ratAtoD(:,:,times(i))))))/9;
    intfeatRplusG(i)=squeeze(nansum(nansum(masks(:,:,i).*(double(dataApostThresh(:,:,times(i)))+double(dataDpostThresh(:,:,times(i))))/2)))/9;

end
IntPeakNeighbour{j}=intfeat;
IntPeakRplusGNeighbour{j}=intfeatRplusG;
%%%%%%%%

fited=polyfit(times*FrameTime,intfeat',1);
LengthTrajectory(j)=length(times);
SlopeIntensity(j)=fited(1);

fited=polyfit(times*FrameTime,intfeatRplusG',1);
SlopeIntensityRplusG(j)=fited(1);


           if ishandle(h)
               waitbar((j+1)/length(indgoodtraj),h)
           else 
               break
           end

end
if ishandle(h)
close(h)
end
set(gcbf,'pointer','arrow');

clear dataApostThresh dataDpostThresh ratAtoD
save([PathName FileName1(1:length(FileName1)-4) '_' FileName2(1:length(FileName2)-4) '_stats.mat'])
