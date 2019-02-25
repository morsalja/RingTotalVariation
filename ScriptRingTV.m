clc
clear

load('FBPNegative40.mat')

x=FBPNegative40;
lambdaIn =0.01; 
alpIn = 2; 
Center=[256,256];
itterRTV_IN=30;   

y= RingTotalVariation(x,lambdaIn,alpIn,Center,itterRTV_IN);

figure(1)
imshow([x y],[0 0.7])

figure(2)
imshow(abs(x-y),[])
