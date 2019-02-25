clc
clear

load('FBPNegative40.mat')

x=FBPNegative40;
lambdaIn =0.01; %weight of the TV term in the cost function (see eq(7) in the manuscript)
alpIn = 2; % 3 and 4 is the best then 2 then 1,2 is the best, 1 is the best
Center=[256,256];% CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE CHANGEABLE
itterRTV_IN=30;   

y= RingTotalVariation(x,lambdaIn,alpIn,Center,itterRTV_IN);

figure(1)
imshow([x y],[0 0.7])

figure(2)
imshow(abs(x-y),[])