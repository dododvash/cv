%%
%ADVANCED IMAGE PROCESSING LAB
%STUDENT : DAVID DVASH
%DATE : 22.05.2021

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%`````````````Code Project``````````````%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Parameters:
r=1.18; %Power Parameter of Masi Entropy
R=1; G=2; B=3; % R G B Valuses
L=256; % Num of Levels Colors
cnt=1; % Aling Counter for images
NUMofimage=5; % number of image we want create

while (cnt<NUMofimage+1)
% states:
if cnt==1
Image=imread('Capture8.PNG'); % Read Image to RGB format
elseif cnt==2
Image=imread('Capture5.PNG'); % Read Image to RGB format
elseif cnt==3
Image=imread('Capture4.PNG'); % Read Image to RGB format
elseif cnt==4
Image=imread('Capture6.PNG'); % Read Image to RGB format
elseif cnt==5
Image=imread('Capture7.PNG'); % Read Image to RGB format
end
tic
%Step1 - Create Entropy and Histograma
[ROW,COM]=size(Image(:,:,1));
SegmentationImage=ones(ROW,COM,3);

%RED Color:
ImageHistR=imhist(Image(:,:,R)); % Create Histogram
SumImageHistR=sum(ImageHistR); % Sumation Of Histogram 
%Create Probabilities Of Image RGB:
for i=1:L
    PR(i)=ImageHistR(i)/SumImageHistR; % P is Arry with Probabilities
end
%GREEN Color:
ImageHistG=imhist(Image(:,:,G)); % Create Histogram
SumImageHistG=sum(ImageHistG); % Sumation Of Histogram 
%Create Probabilities Of Image RGB:
for i=1:L
    PG(i)=ImageHistG(i)/SumImageHistG; % P is Arry with Probabilities
end
%BLUE Color:
ImageHistB=imhist(Image(:,:,B)); % Create Histogram
SumImageHistB=sum(ImageHistB); % Sumation Of Histogram 
%Create Probabilities Of Image RGB:
for i=1:L
    PB(i)=ImageHistB(i)/SumImageHistB; % P is Arry with Probabilities
end

figure (cnt+NUMofimage);
subplot(4,1,1);
imshow(uint8(Image));
title('Original Image');
subplot(4,1,2);
p=stem(PR);
title('Probabilities Of Image RED Color');
p.Color = 'r';
subplot(4,1,3);
p=stem(PG);
title('Probabilities Of Image GREEN Color');
p.Color = 'g';
subplot(4,1,4);
p=stem(PB);
title('Probabilities Of Image BLUE Color');
p.Color = 'b';

%Step 2 - Create Masi Entropy

[EjR , EjG , EjB] = Step2_MasiEntropy(Image,r);

%Step 3 - Create Optimal Threshold
[thR,thG,thB] = Step3_OptimalThreshold(Image, EjR , EjG , EjB);
tempR= MY_FUNC(Image,thR,1);
tempG= MY_FUNC(Image,thG,2);
tempB= MY_FUNC(Image,thB,3);

%Step 4 -
%Marge Red, Green, Blue value:
SegmentationImage = Step4_MargeRGBvalue(SegmentationImage, tempR , tempG , tempB);
%Create Histogarm in RGB
% To Original Image:
hist_Image = Step4_HistogarmRGB(Image);
% To Segmentation Image:
hist_Sigment = Step4_HistogarmRGB(SegmentationImage);
toc
%step 5 plots
figure (cnt);
subplot(2,2,1);
imshow(uint8(SegmentationImage));
title(['Segmentation Image, r=' num2str(r)]);

subplot(2,2,2);
imshow(uint8(Image));
title('Original Image');

subplot(2,2,3);
p=plot(0:255,hist_Sigment(:,1),0:255,hist_Sigment(:,2),0:255,hist_Sigment(:,3));
p(1).Color = 'r'; p(2).Color = 'g'; p(3).Color = 'b';
title(['Histogram Segmentation, r=' num2str(r)]);

subplot(2,2,4);
p=plot(0:255,hist_Image(:,1),0:255,hist_Image(:,2),0:255,hist_Image(:,3));
p(1).Color = 'r'; p(2).Color = 'g'; p(3).Color = 'b';
title('Histogram Original');

cnt=cnt+1;
end
%%
clear all;
close all;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%`````````````FUNCTION``````````````%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function temp = MY_FUNC(Image,th,n)
th1=th;
[ROW,COM]=size(Image(:,:,n));
 temp=ones(ROW,COM,1);
Ltrh=length(th1);
for i=1:1:ROW
    for j=1:1:COM
        for t=0:Ltrh-1
            if (Image(i,j,n) < th1(Ltrh-t))
                temp(i,j)=th1(Ltrh-t);
            end
            if (Image(i,j,n) > th1(Ltrh))
                temp(i,j)=255;  
            end
        end
    end
end
end


function [EjR , EjG , EjB] = Step2_MasiEntropy(Image,r)
R=1; G=2; B=3; % R G B Valuses
%RED Color:
TEMPR=entropy(Image(:,:,R)); % Create regular Entropy
EjR=(-1)*((1/(1-r))*log(1-(1-r)*TEMPR)); % Create Masi Entropy
%GREEN Color:
TEMPG=entropy(Image(:,:,G)); % Create regular Entropy
EjG=(-1)*((1/(1-r))*log(1-(1-r)*TEMPG)); % Create Masi Entropy
%BLUE Color:
TEMPB=entropy(Image(:,:,B)); % Create regular Entropy
EjB=(-1)*((1/(1-r))*log(1-(1-r)*TEMPB)); % Create Masi Entropy

end


function [thR,thG,thB] = Step3_OptimalThreshold(Image, EjR , EjG , EjB)
R=1; G=2; B=3; % R G B Valuses
%RED Color:
EjR=round(EjR);
thR=multithresh(Image(:,:,R),EjR);% Multilevel image thresholds using Otsu’s method
ImageR=Image(:,:,R);

%GREEN Color:
EjG=round(EjG);
thG=multithresh(Image(:,:,G),EjG);% Multilevel image thresholds using Otsu’s method

%BLUE Color:
EjB=round(EjB);
thB=multithresh(Image(:,:,B),EjB);% Multilevel image thresholds using Otsu’s method

%th=[thR,thG,thB];
end


function SegmentationImage = Step4_MargeRGBvalue(SegmentationImage, tempR , tempG , tempB)
R=1; G=2; B=3; % R G B Valuses
SegmentationImage(:,:,R)=tempR;
SegmentationImage(:,:,G)=tempG;
SegmentationImage(:,:,B)=tempB;

end


function hist_Image = Step4_HistogarmRGB(Image)
R=1; G=2; B=3; % R G B Valuses
%Step 4 - Create Histogarm in RGB
% To Original Image:
hist_Image_R=imhist(uint8(Image(:,:,R)));
%hist_Image_R=cumsum(hist_Image_R);

hist_Image_G=imhist(uint8(Image(:,:,G)));
%hist_Image_G=cumsum(hist_Image_G);

hist_Image_B=imhist(uint8(Image(:,:,B)));
%hist_Image_B=cumsum(hist_Image_B);

hist_Image=[hist_Image_R ,hist_Image_G ,hist_Image_B ];
end
