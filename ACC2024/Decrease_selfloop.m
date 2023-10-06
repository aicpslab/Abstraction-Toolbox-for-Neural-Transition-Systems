%% Initialization
clear all 
clc
close all
Num=100;
%% Load Data and Reverse mapminmax
names = {'Angle.mat','CShape.mat','LShape.mat','PShape.mat','Spoon.mat'};
for n=1:5
      load(['Results/' names{n}]) %loading the model
    for i=1:size(intervals,2)
        intvl_lower = mapminmax('reverse',intervals{i}(:,1),ps_input);
        intvl_upper = mapminmax('reverse',intervals{i}(:,2),ps_input);
        intervals{i}=[intvl_lower,intvl_upper];        
    end

%% Data Processing
for i= 1:size(IntersectSet,1)
   if IntersectSet(i,i)==1
      IntersectSet(i,i)=Selfless(demos,intervals{i},Num);
   end
end

%% Let's see how it works 
ModelGraph=digraph(IntersectSet);
figure
plot(ModelGraph,'r')
end