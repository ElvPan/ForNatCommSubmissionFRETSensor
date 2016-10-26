% Load images in both acceptor and donor channels
[FileName1,PathName] = uigetfile('*.tif','Select the acceptor channel image');
dataA=rd_img16([PathName FileName1]);
[FileName2,PathName] = uigetfile('*.tif','Select the acceptor channel image');
dataD=rd_img16([PathName FileName2]);
%%%%% Prompt user for some parameters
prompt = {'Enter pixel size (microns):'};
dlg_title = 'Enter Analysis Parameters';
num_lines = 1;
defaultans = {num2str(0.1)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
pixelsize=str2double(answer{1}); % in um

% combine data from Donor and Acceptor
for i=1:size(dataA,3)
combDandA(:,:,i)=double(dataA(:,:,i))+double(dataD(:,:,i));
end
% perform ICS on combined images
[G] = corrfunc(combDandA);
% crop correlation function around central peak
[cropG,rectG] = autocrop(G,20);
% fit croped corr. fn.
[a] = gaussfitold(cropG,'2d',pixelsize,'y');
% calcualte Degree of aggregation
for i=1:size(a,1)
    aveim(i)=squeeze(mean(mean(combDandA(:,:,i),1),2));
    DA(i)=aveim(i)*(pi*a(i,1)*a(i,2)*a(i,3));
end