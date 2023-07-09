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
Bperp0=load('Bperp36.out');
Tm0=load('Tm36.out');
deltime0=load('deltime36.out'); 
timedeltas0=load('timedeltas36.out');


Bperp1=zeros(size(Bperp0));
Tm1=zeros(size(Tm0));
deltime1=zeros(size(deltime0));


cells2=importdata('sbas_list_active2');
%cells2=importdata('sbas_list36');
N2=length(cells2);

sbas_dates0=string(cells2.textdata);
sbas_data0=string(cells2.data);
cells=importdata('intlist_active2');
%cells=importdata('intlist36');
N=length(cells);
lambda=5.6; %wavelength

unw_phase=zeros(nr,naz,N);
%uni=zeros(nr,naz,N);
amps=zeros(nr,naz,N);
ints=zeros(nr,naz,N);
coh=zeros(nr,naz,N);
%uni=zeros(nr,naz,N);
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

% Get an average amplitude file 
% amp_mean = mean(amp,3);
% amp_mask=amp_mean;
% amp_mask(amp_mask<0.1)=nan;
% amp_mask(amp_mask>0.1)=1;

%generated combined amplitude and coherence mask
mask=cor_mask;
%mask(mask<2)=nan;
%mask(mask==2)=1;

%clear amp
%clear coh

%phase_cal=unw_phase;
% 
% 
% %view time series
% figure('units','normalized','outerposition',[0 0 1 1]);
% for i=1:N
%     %strint=cells{i};   
%     %filename=sprintf('%s',strint); 
%     %figure('units','normalized','outerposition',[0 0 1 1]);
%     %subplot(2,1,1),imagesc(angle(ints(:,:,i)).');
%     subplot(2,1,1),imagesc(unw_phase(:,:,i).');
%     title(strcat(string(i),'th interferogram ',string(date_pair(1,i)),'-',string(date_pair(2,i))));
%     %ylim([460 640]);
%     %xlim([1500 1900]);
%     %caxis([-pi pi]);
%     subplot(2,1,2),imagesc(coh(:,:,i).');
%     title(strcat(string(i),'th interferogram ',string(date_pair(1,i)),'-',string(date_pair(2,i))));
%     caxis([0 1]);
%     %ylim([460 640]);
%     %xlim([1500 1900]);
%     pause(0.1);
%    % close all
%  end

%% Closure Matrices

%Set up SBAS
dates=load('dates');
dtk_old=timedeltas0;
Tm_old=Tm0;
Ndate=length(dates);
cellssbas=cells2;
Tm=zeros(N,Ndate-1);
A=zeros(N,Ndate);
%Tm=zeros(size(Tm));
%A=zeros(size(Tm));
dtk=dtk_old;

for i=1:N
    disp(i)
    strint=cells{i};   
    filename=sprintf('%s',strint); 
    s1=split(filename,'_');
    s2=string(s1(2));
    s3=split(s2,'.');
    date1=string(s1(1));
    date2=string(string(s3(1)));
    dates_vec(i,1)=double(date1);
    dates_vec(i,2)=double(date2);
    
    sbas_dates=cellssbas.textdata;
    sbas_dates_1=string(sbas_dates(:,1));
    sbas_dates_2=string(sbas_dates(:,2));
    
    return1a=find(contains(sbas_dates_1,date1));
    return1b=find(contains(sbas_dates_2,date2));
    ind1=intersect(return1a,return1b);
    
    sbas_data=cellssbas.data;
    temps=sbas_data(:,1);
    perps=sbas_data(:,2);
    Btemp(i)=temps(ind1);
    Bperp(i)=perps(ind1);
    
    Tm(i,:)=Tm_old(ind1,:);
    
    date_ind1=find(contains(string(dates),date1));
    date_ind2=find(contains(string(dates),date2));
    A(i,date_ind1)=-1;
    A(i,date_ind2)=1;
    
    %dtk(i,:)=dtk_old(ind1,:);
end



J=nchoosek(1:1:n,2);
K=nchoosek(1:1:n,3);

J1=J(:,1);
J2=J(:,2);

%Chi=zeros(nr,naz,length(K));
B1=zeros(length(K),n);
bads=[];

for i=1:length(K)
    inds=K(i,:);
    
    return1a=find(dates_vec(:,1)==dates(inds(1)));
    return1b=find(dates_vec(:,2)==dates(inds(2)));
    ind1=intersect(return1a,return1b);
    
    return2b=find(dates_vec(:,1)==dates(inds(2)));
    return2c=find(dates_vec(:,2)==dates(inds(3)));
    ind2=intersect(return2b,return2c);
    
    return3a=find(dates_vec(:,1)==dates(inds(1)));
    return3c=find(dates_vec(:,2)==dates(inds(3)));
    ind3=intersect(return3a,return3c);
    
    size_check=[size(ind1)==[1,1],size(ind2)==[1,1],size(ind3)==[1,1]];
    %if sum(size_check)==6
    %    Chi(:,:,i) = angle((phase(:,:,ind1).*phase(:,:,ind2).*conj(phase(:,:,ind3))));
    %elseif sum(size_check)<6
    %    bads=[bads i];
    %end   
    B1(i,ind1)=1;
    B1(i,ind2)=1;
    B1(i,ind3)=-1;
end

B=B1;
B(bads,:)=[];
%Chi_cull=Chi;
%Chi_cull(:,:,bads)=[];

%%
% %compute svd of G
% G=B;
% [U S V]=svd(G);
% k=n-1;
% Si=zeros(rank(G)-k,rank(G)-k);
% for i=1:rank(G)-k
% %keep all the singular values
%  if abs(S(i,i))<1.e-8
%    Si(i,i)=0;
%  else
%    Si(i,i)=1./S(i,i);
%  end
% end
% 
% Vp=V(:,1:(rank(G)-k));
% Up=U(:,1:(rank(G)-k));
% Bi=Vp*Si*Up';
% 
% 
% for i=1:nr
%     for j=1:naz
%         pChi(i,j,:)=Bi*squeeze((Chi_cull(i,j,:)));            
%     end
% end
% 
% 
% 
% %Correct igrams with corrected stack pChi
% igrams=angle(phase);
% igrams_coru=wrapToPi(igrams-pChi);
% igrams_cor=(igrams-pChi);


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
     
%     block_phase = masked_unw_phase1(:,:,int);
%     indx = isnan(block_phase);
%     line = polyfit(dem(~indx),block_phase(~indx),1);
%     correction = (line(1)*dem + line(2)); 
%     corrections1(:,:,int)=correction;
%     %phase1(:,:,int) = masked_unw_phase1(:,:,int);
%     phase1(:,:,int) = masked_unw_phase1(:,:,int) - correction;
    
end
% 
% 
% %view time series
% figure('units','normalized','outerposition',[0 0 1 1]);
% for i=1:N
%     %strint=cells{i};   
%     %filename=sprintf('%s',strint); 
%     %figure('units','normalized','outerposition',[0 0 1 1]);
%     %subplot(2,1,1),imagesc(angle(ints(:,:,i)).');
%     subplot(2,1,1),imagesc(phase(:,:,i).');
%     title(strcat(string(i),'th interferogram ',string(date_pair(1,i)),'-',string(date_pair(2,i))));
%     %ylim([460 640]);
%     %xlim([1500 1900]);
%     %caxis([-pi pi]);
%     subplot(2,1,2),imagesc(coh(:,:,i).');
%     title(strcat(string(i),'th interferogram ',string(date_pair(1,i)),'-',string(date_pair(2,i))));
%     caxis([0 1]);
%     %ylim([460 640]);
%     %xlim([1500 1900]);
%     pause(0.1);
%    % close all
%  end

%clear masked_unw_phase
%clear unw_phase
%clear ints

%% Pluck out Scenes you want to use, only use scenes that satisfy a skewness criteria
% val=0.0375;
% for i=1:N
% ff=coh(:,:,i);
% skew(i)=nanvar(ff(:));
% end
% skew_val=skew;
% skew_val(abs(skew_val)>val)=nan;
% skew_val(abs(skew_val)<=val)=1;

%scene_inds=find(skew_val==1);
scene_inds=linspace(1,N,N);
Ngood=length(scene_inds);
phase_good=phase(:,:,scene_inds);
coh_good = coh(:,:,scene_inds);
doy_pair_good=doy_pair(:,scene_inds);

%clear coh
%clear phase
%clear block_phase
%clear corrections
%clear dem0
%clear masked_unw_phase
%clear unw_phase
%clear phase1

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




% 
% for i=1:Ngood
%     phase_good_std_closed(:,:,i) = imopen(phase_good_std2(:,:,i),se);
% end

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

% i=1;
% figure,subplot(2,1,1),imagesc(a_1(:,:,i).');
% hold on
% subplot(2,1,2),imagesc(a_2(:,:,i).');

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


%B(14,:)=0;
%B(:,7)=0;
%B=B(4:15,2:end);
%Bi=pinv(B);
Bi=pinv(-1.*B);
%I2=I(:,:,4:15);
%a=a(:,:,2:end);


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


%% 'SBAS' for Fire Burn Area original

% I=masked_opened;
% I(isnan(I))=0;
% 
% A4   = I(:,:,1);
% A4_mean=A4;
% A4_mean(A4_mean<1)=0;
% A5_1 = I(:,:,2);
% A5_2 = I(:,:,3)-I(:,:,1);
% A5_sum=A5_1+A5_2;
% A5_mean=A5_sum;
% A5_mean(A5_mean<2)=0;
% A6   = I(:,:,4)-I(:,:,1);
% A6_mean=A6;
% A6_mean(A6_mean<1)=0;
% A7_1 = I(:,:,5)-I(:,:,1);
% A7_2 = I(:,:,6)-I(:,:,2);
% A7_3 = I(:,:,7)-I(:,:,4)-I(:,:,1);
% A7_sum=A7_1+A7_2+A7_3;
% A7_mean=A7_sum;
% A7_mean(A7_mean<3)=0;
% A8_1 = I(:,:,8)-I(:,:,1);
% A8_2 = I(:,:,9)-I(:,:,2);
% A8_3 = I(:,:,10)-I(:,:,4)-I(:,:,1);
% A8_sum=A8_1+A8_2+A8_3;
% A8_mean=A8_sum;
% A8_mean(A8_mean<3)=0;
% 
% 
% A11_1_1= I(:,:,12)-A7_1;
% A11_1_2= I(:,:,12)-A7_2;
% A11_1_3= I(:,:,12)-A7_3;
% A11_1_sum = A11_1_1+A11_1_2+A11_1_3;
% 
% A11_2_1= I(:,:,13)-A8_1;
% A11_2_2= I(:,:,13)-A8_2;
% A11_2_3= I(:,:,13)-A8_3;
% A11_2_sum = A11_2_1+A11_2_2+A11_2_3;
% 
% A11_sum = A11_1_sum+A11_2_sum;
% A11_mean=A11_sum./6;
% A11_mean_mask=A11_mean;
% A11_mean_mask(A11_mean_mask<1)=0;
% A11_mean_mask(A11_mean_mask>1)=0;
% 
% 
% A9   = -(-I(:,:,14) - A11_mean_mask);
% A9_mean=A9;
% A9_mean(A9_mean>1)=0;
% A9_mean(A9_mean<1)=0;
% 
% A11_3= -(I(:,:,14)-A9);




% descending_map = zeros(nr,naz,n);
% descending_map(:,:,4)=A4_mean;
% descending_map(:,:,5)=A5_mean;
% descending_map(:,:,6)=A6_mean;
% descending_map(:,:,7)=A7_mean;
% descending_map(:,:,8)=A8_mean;
% descending_map(:,:,9)=A9_mean;      %10/13
% descending_map(:,:,11)=A11_mean_mask; %1106
% 
% save descending_map descending_map -v7.3

%%  Do Map Progression

%ascending_map=load('../subdir_ascending/ascending_map.mat');
%load ascending_map.mat

%Descending Dates
% 1   20200709
% 2   20200721
% 3   20200802
% 4   20200814      x
% 5   20200826      x
% 6   20200907      x
% 7   20200919      x
% 8   20201001      x
% 9   20201013      x
% 10  20201025
% 11  20201106      x
% 12  20201118
% 13  20201130
% 14  20201212
% 15  20201224
% 16  20210105
% 17  20210117

% Ascending Dates
% 1   20200705
% 2   20200717
% 3   20200729
% 4   20200810
% 5   20200822      x
% 6   20200903      x
% 7   20200915      x
% 8   20200927      x
% 9   20201009      x
% 10  20201021      x
% 11  20201114
% 12  20201126
% 13  20201208
% 14  20210101
% 15  20210113
% 16  20210125


%Total Dates
%1   1   20200705      a
%2   1   20200709      d
%3   2   20200717      a
%4   2   20200721      d
%5   3   20200729      a
%6   3   20200802      d
%7   4   20200810      a
%8   4   20200814      dx
%9   5   20200822      ax
%10  5   20200826      dx
%11  6   20200903      ax
%12  6   20200907      dx
%13  7   20200915      ax
%14  7   20200919      dx
%15  8   20200927      ax
%16  8   20201001      dx
%17  9   20201009      ax
%18  9   20201013      dx
%19  10  20201021      ax
%20  10  20201025      d
%21  11  20201106      dx
%22  11  20201114      a
%23  12  20201118      d
%24  12  20201126      a
%25  13  20201130      d
%26  13  20201208      a
%27  14  20201212      d
%28  15  20201224      d
%29  14  20210101      a
%30  16  20210105      d
%31  15  20210113      a
%32  17  20210117      d
%33  16  20210125      a


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

