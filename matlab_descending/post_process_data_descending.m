%%Convention is Master_slave.unw, and crossmul takes master*conj(slave),
%%such that a negative values in the unwrapped phase corresponds to a
%%subsidence (since phase=phi_1-phi_2)

%Fire began 8/13/2020, was 100% contained 12/2/20
%greatest extent by 10/18/20 ?
%Data processed was Path 56 Frame 455
%pixel spacing in dem.rsc is 15m

%% Load all of the InSAR files you are considering
clear all; close all; clc;

[A geoR]=geotiffread('20200709_20200919.tif');

nr=8222; naz=3420;% image size:  nr = # of x (range) pixels; 
                             % naz = # of y (azimuth) pixels
n=17; % number of slcs 

%Just read in coherent inteferograms with a temporal
%baseline of 48 days or less
Bperp0=load('Bperp36_descending.out');
Tm0=load('Tm36_descending.out');
deltime0=load('deltime36_descending.out'); 
timedeltas0=load('timedeltas36_descending.out');


Bperp1=zeros(size(Bperp0));
Tm1=zeros(size(Tm0));
deltime1=zeros(size(deltime0));


cells2=importdata('sbas_list_active2_descending');
N2=length(cells2);

sbas_dates0=string(cells2.textdata);
sbas_data0=string(cells2.data);
cells=importdata('intlist_active2_descending');
N=length(cells);
lambda=5.6; %wavelength

unw_phase=zeros(nr,naz,N);
amps=zeros(nr,naz,N);
ints=zeros(nr,naz,N);
coh=zeros(nr,naz,N);
date_pair=cell(2,N);
doy_pair=cell(2,N);

%Read in the unwrapped phase (unw), coherence (coh), amplitude (amp) and
%unimodally-corrected unwrapped phase (uni)
for i=1:N
    disp(i)
    strint=cells{i};
    strint1=strcat('ints/',strint);
    strunw1=strrep(strint,'.int','.unw');
    strunw=strcat('unws/',strunw1);
    stramp1=strrep(strint,'.int','.amp');
    stramp=strcat('amps/',stramp1);
    strcc1=strrep(strint,'.int','.cc'); 
    strcc=strcat('ccs/',strcc1);
    struni1=strrep(strint,'.int','.int');
    struni=strcat('ints/',struni1);
% correlations
    filename_c=sprintf('%s',strcc); 
    fid=fopen(filename_c);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat((nr+1):end,:);
    coh(:,:,i)=temp;
    fclose(fid);
% unwrapped phase
    filename=sprintf('%s',strunw);  
    fid=fopen(filename);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat(nr+1:end,:);
    unw_phase(:,:,i)=temp;
    fclose(fid);
% Interferograms
   filename=sprintf('%s',struni);  
   fid=fopen(filename);
   dat=fread(fid,[2*nr,inf],'float','ieee-le');
   temp=dat(1:2:end,1:naz)+1i*dat(2:2:end,1:naz);
   phase(:,:,i)=temp;
   fclose(fid);
% Amplitude
    filename=sprintf('%s',strint1);  
    fid=fopen(filename);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat(1:2:2*nr-1,:)+1i*dat(2:2:2*nr,:);
    ints(:,:,i)=temp;
    fclose(fid);
% Amplitude
    filename=sprintf('%s',stramp);  
    fid=fopen(filename);
    dat=fread(fid,[2*nr,inf],'float','ieee-le');
    temp=dat(1:2:2*nr-1,:)+1i*dat(2:2:2*nr,:);
    amps(:,:,i)=temp;
    fclose(fid);
% date information
split1=strsplit(strint,'_');
strint2=split1{2};
split2=strsplit(strint2,'.');
d1=split1{1};
d2=split2{1};

date1=strcat(d1(5:6),'/',d1(7:8),'/',d1(1:4));
date2=strcat(d2(5:6),'/',d2(7:8),'/',d2(1:4));

date1_vec=datetime(date1,'InputFormat','MM/dd/yyyy');
date2_vec=datetime(date2,'InputFormat','MM/dd/yyyy');

doy1=day(date1_vec,'dayofyear');
doy2=day(date2_vec,'dayofyear');

date_pair{1,i}=date1;
date_pair{2,i}=date2;
doy_pair{1,i}=doy1;
doy_pair{2,i}=doy2;
end


% Get an average coherence file 
avecc = mean(coh,3);
cor_mask=avecc;
cor_mask(cor_mask<0.5)=nan;
cor_mask(cor_mask>0.5)=1;

%generated coherence mask
mask=cor_mask;


%% Now we will do the atmospheric correction
%remove topographically-correlated atmospheric noise from interferograms
% Depending on how large the scene is, you may want to base the correction
% solely on a few areas with the topographic relief changes quite a bit.
% But if the scene is small, you can indeed use the whole scene.
% This bit also requires a calibration pixel (or set of pixels) to ensure
% that all of the scenes are set to the same "datum"

% Load the dem
nr0=24667;
naz0=10261;
fid=fopen('elevation.dem','r');
dem0=fread(fid,[nr0,naz0],'int16'); % x length first, y length second
fclose(fid);

dem=imresize(dem0,[nr,naz]);
masked_unw_phase = unw_phase;%.*mask;
%masked_unw_phase1 = uni.*mask;

% Pick some pixels for calibration
pixels = [1140 1070]; % [range azimuth]
[sz2,~] = size(pixels);

phase = zeros(nr,naz,N);
phase1 = zeros(nr,naz,N);
corrections=zeros(nr,naz,N);
phase_rshp = [];
dem_rshp = [];
for int = 1:N
     block_phase = masked_unw_phase(:,:,int);
     indx = isnan(block_phase);
     line = polyfit(dem(~indx),block_phase(~indx),1);
     correction = (line(1)*dem + line(2)); 
     corrections(:,:,int)=correction;
     %phase(:,:,int) = masked_unw_phase(:,:,int);
     phase(:,:,int) = masked_unw_phase(:,:,int) - correction;
     
end


%% Pluck out Scenes you want to use, only use scenes that satisfy a skewness criteria

%scene_inds=find(skew_val==1);
scene_inds=linspace(1,N,N);
Ngood=length(scene_inds);
phase_good=phase(:,:,scene_inds);
coh_good = coh(:,:,scene_inds);
doy_pair_good=doy_pair(:,scene_inds);

%% Phase Standard Deviation

h = 1/11*ones(11,1);
H = h*h';

for i=1:Ngood
    phase_good_std(:,:,i) = stdfilt(phase_good(:,:,i),true(11));
    coh_good_std(:,:,i) = stdfilt(coh_good(:,:,i),true(11));
    phase_good_mean(:,:,i) = filter2(H,phase_good(:,:,i));
    coh_good_mean(:,:,i) = filter2(H,coh_good(:,:,i));
end




phase_good_std2=phase_good_std;
val=1.75;
val=pi/2;
val=pi;
val=sqrt(pi);
val=pi/sqrt(3);
%val=1.5;
phase_good_std2(phase_good_std2<val)=0;
phase_good_std2(phase_good_std2>=val)=1;
se = strel('disk',10);


%4489 corresponds to ~ 1km square area
for i=1:Ngood
    opened(:,:,i)=bwareaopen(phase_good_std2(:,:,i),4489);
end


%std filter mask
test=nanvar(opened,0,3);
test_mask=test;
test_mask(test_mask>0.12)=1;
test_mask(test_mask<0.12)=nan;


bounding_mask=ones(nr,naz);
bounding_mask(1:915,:)=nan;
bounding_mask(5906:end,:)=nan;

bounding_mask(1500:2250,1:400)=nan;
bounding_mask(1:1600,1:900)=nan;
bounding_mask(1:500,1600:2400)=nan;
bounding_mask(600:1950,2660:end)=nan;
bounding_mask(1930:2300,2875:end)=nan;


masked_opened=(bounding_mask.*test_mask.*opened);
masked_opened1=masked_opened;
masked_opened1(isnan(masked_opened1))=0;
masked_opened(masked_opened==0)=nan;
masked_phase=(bounding_mask.*test_mask.*opened.*phase_good);
masked_phase(masked_phase==0)=nan;
masked_std=(bounding_mask.*test_mask.*opened.*phase_good_std);
masked_std(masked_std==0)=nan;
masked_coh=(bounding_mask.*test_mask.*opened.*coh_good);
masked_coh(masked_coh==0)=nan;


for i=1:N
    masked_phase_mean(i)=nanmean(nanmean(masked_phase(:,:,i)));
    masked_std_mean(i)=nanmean(nanmean(masked_std(:,:,i)));
end

%%

%get converse map, all non fire pixels
total_fire_map=bounding_mask.*test_mask.*opened;
total_fire_map(isnan(total_fire_map))=0;
nonfire_map=1-total_fire_map;
nonfire_map(isnan(nonfire_map))=0;

masked_nonfire_coh=nonfire_map.*coh_good;

%view time series
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:Ngood
    subplot(3,1,1),imagesc(masked_phase(:,:,i).');
    title(strcat(string(i),'th interferogram_',string(date_pair(1,i)),'-',string(date_pair(2,i))));
    subplot(3,1,2),imagesc(phase_good(:,:,i).');
    title(strcat(string(i),'th interferogram_',string(date_pair(1,i)),'-',string(date_pair(2,i))));
    subplot(3,1,3),imagesc(coh_good(:,:,i).');
    title(strcat(string(i),'th interferogram_',string(date_pair(1,i)),'-',string(date_pair(2,i))));
    pause(0.15);
    new1 = strrep(string(date_pair(1,i)),'/','_');
    new2 = strrep(string(date_pair(2,i)),'/','_');
    st=strcat('view_descending_',string(i),'th interferogram_',new1,'-',new2);
    geotiffwrite(st,masked_phase(:,:,i).',geoR);
end
 
%% 'SBAS' for Fire Burn Area new
save descending_masked_opened masked_opened -v7.3

I=masked_opened;
I(isnan(I))=0;

I2=I(:,:,1:13);

%a is 8 x 1
a=zeros(size(I2,1),size(I2,2),8);
%B is 14 x 8, #ints x #scenes
B=zeros(13,8);
B(1,1)=1;
B(2,1)=1;
B(3,2)=1;
B(4,2)=1;
B(5,2)=1;
B(6,3)=1;
B(7,4)=1;
B(8,2)=1;
B(9,3)=1;
B(10,4)=1;
B(11,5)=1;
B(12,5)=1;
B(13,6)=1;
%B(14,7)=1;

B(1,2)=-1;
B(2,3)=-1;
B(3,3)=-1;
B(4,4)=-1;
B(5,5)=-1;
B(6,5)=-1;
B(7,5)=-1;
B(8,6)=-1;
B(9,6)=-1;
B(10,6)=-1;
B(11,6)=-1;
B(12,8)=-1;
B(13,8)=-1;
%B(14,8)=-1;


%B(14,:)=0;
%B(:,7)=0;

Bi=pinv(-1.*B);

for i=1:size(I,1)
    for j=1:size(I,2)
        pluck=squeeze(I2(i,j,:));
        out=Bi*pluck;
        a(i,j,:)=out;
    end
end

% a_filt2=a;
% a_filt2(a_filt2<0)=0;
% a_2=a_filt2;
% a_2(a_2>0)=1;


a_filt=a;
a_filt(a_filt<0)=0;
filt_val=sum(abs(B));
a_1=zeros(size(a_filt));

for i=1:length(filt_val)
    pluck=a_filt(:,:,i);
    val=0.25./filt_val(i);
    pluck(pluck<val)=0;
    pluck(pluck>=val)=1;
    a_1(:,:,i)=pluck;
end


%% ascending 'SBAS' for Fire Burn Area new
load ascending_masked_opened.mat

I=ascending_masked_opened;
I(isnan(I))=0;
I2=I;
%a is 8 x 1
a=zeros(size(I2,1),size(I2,2),9);
%B is 14 x 8, #ints x #scenes
B=zeros(18,9);
B(1,1)=1;
B(2,3)=1;
B(3,1)=1;
B(4,2)=1;
B(5,3)=1;
B(6,4)=1;
B(7,5)=1;
B(8,3)=1;
B(9,4)=1;
B(10,5)=1;
B(11,6)=1;
B(12,4)=1;
B(13,5)=1;
B(14,6)=1;
B(15,5)=1;
B(16,6)=1;
B(17,7)=1;
B(18,8)=1;

h=-1;
B(1,4)=h;
B(2,4)=h;
B(3,5)=h;
B(4,6)=h;
B(5,6)=h;
B(6,6)=h;
B(7,6)=h;
B(8,7)=h;
B(9,7)=h;
B(10,7)=h;
B(11,7)=h;
B(12,8)=h;
B(13,8)=h;
B(14,8)=h;
B(15,9)=h;
B(16,9)=h;
B(17,9)=h;
B(18,9)=h;


Bi=pinv(-1.*B);


for i=1:size(I,1)
    for j=1:size(I,2)
        pluck=squeeze(I2(i,j,:));
        out=Bi*pluck;
        a(i,j,:)=out;
    end
end

a_filt=a;
a_filt(a_filt<0)=0;
filt_val=sum(abs(B));
a_2=zeros(size(a_filt));

for i=1:length(filt_val)
    pluck=a_filt(:,:,i);
    val=1.5./filt_val(i);
    pluck(pluck<val)=0;
    pluck(pluck>=val)=1;
    a_2(:,:,i)=pluck;
end

ascending_map = zeros(nr,naz,n);
ascending_map(:,:,5)=a_2(:,:,4);
ascending_map(:,:,6)=a_2(:,:,5);
ascending_map(:,:,7)=a_2(:,:,6);
ascending_map(:,:,8)=a_2(:,:,7);
ascending_map(:,:,9)=a_2(:,:,8);
ascending_map(:,:,10)=a_2(:,:,9);

total_map=zeros(nr,naz,33);
total_map(:,:,4)=a_1(:,:,1);                  %07/21
total_map(:,:,8)=a_1(:,:,2);                  %08/14
total_map(:,:,9)=ascending_map(:,:,5);        %08/22
total_map(:,:,10)=a_1(:,:,3);                 %08/26
total_map(:,:,11)=ascending_map(:,:,6);       %09/03
total_map(:,:,12)=a_1(:,:,4);                 %09/07
total_map(:,:,13)=ascending_map(:,:,7);       %09/15
total_map(:,:,14)=a_1(:,:,5);                 %09/19
total_map(:,:,15)=ascending_map(:,:,8);       %09/27
total_map(:,:,16)=a_1(:,:,6);                 %10/01
total_map(:,:,17)=ascending_map(:,:,9);       %10/09
total_map(:,:,18)=a_1(:,:,7);                 %10/13
total_map(:,:,19)=ascending_map(:,:,10);      %10/21
total_map(:,:,21)=a_1(:,:,8);                 %11/06

total_map(total_map==0)=nan;

    
cellsNN=importdata('total_dates');

%%

%view time series
%figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:33
    imagesc((total_map(:,:,i)).');
    title(strcat(string(i),'th time step ',string(filename)));
    pause(0.1);
    st=strcat('view_time_step_descending_',string(i),'_date_',string(cellsNN(i)));
    geotiffwrite(st,total_map(:,:,i).',geoR);
end


i=4;
geotiffwrite(strcat('coherence_example_',string(i),'.tif'),coh_good(:,:,i).',geoR);
geotiffwrite(strcat('phase_example_',string(i),'.tif'),phase_good(:,:,i).',geoR);
geotiffwrite(strcat('phase_std_example_',string(i),'.tif'),phase_good_std(:,:,i).',geoR);

geotiffwrite(strcat('test.tif'),phase_good_std(:,:,i).',geoR);



%% Let's estimate the first order interferometric phase variance from Bamler and Hartl
coh_vec=0:0.1:1;

sig_var = (pi*pi/3) - pi*asin(abs(coh_vec)) + (asin(abs(coh_vec))).^2 - 0.5*dilog((abs(coh_vec)).^2);

sig_var = (pi*pi/3) - pi*asin(abs(coh_vec)) + (asin(abs(coh_vec))).^2 - 0.5*polylog(2,(abs(coh_vec)).^2);


plot((pi*pi/3)*ones(length(coh_vec),1))
hold on
plot(-1.*pi*asin(abs(coh_vec)))
hold on
plot((asin(abs(coh_vec))).^2)
hold on
plot(-0.5*dilog((abs(coh_vec)).^2))

