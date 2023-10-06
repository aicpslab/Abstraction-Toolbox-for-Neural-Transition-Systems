clear
close all
clc
%%

duration = 1000;
tol = 0.01;
maximum_entropy=300;
Dimension=2;
% Initial Input Set

names = {'Angle','CShape','LShape','PShape','Spoon'};
figName = [];
%% Load Data
for n=1:5
   if (n<0 || n>size(names,2))
        disp('Wrong model number!')
        disp('Please try again and type a number between 1-30.')
    elseif n == 0
        return
    else
      load(['DataSet/' names{n}],'demos','dt') %loading the model
   end
%load(['DataSet/Angle'],'demos','dt')   

%% Training of the initial NN

% 1.Structrues of the NN
NeuronNum_switch=40;
NeuronNum_single=200;
% 2.Extract Data
for i = 1:size(demos,2)
    TrajData{i} = demos{i}.pos;
    plot(TrajData{i}(1,:),TrajData{i}(2,:))
    hold on
end
%
TrajNum=size(demos,2);
xs=zeros((size(TrajData{1},1))*TrajNum,2);
t=zeros((size(TrajData{1},1))*TrajNum,2);
for i = 1:TrajNum
    Begin=(size(TrajData{1},2)-1)*(i-1)+1;
    End = (size(TrajData{1},2)-1)*(i-1)+size(TrajData{1},2);
    xs(Begin:End-1,1) = TrajData{i}(1,1:end-1);
    xs(Begin:End-1,2) = TrajData{i}(2,1:end-1);
    t(Begin:End-1,1) = TrajData{i}(1,2:end);
    t(Begin:End-1,2) = TrajData{i}(2,2:end);
end


        % 2. Maybe try BP this time
         % apply mapminmax process
        [xsn, ps_input] = mapminmax(xs',-1,1);
        tn = mapminmax('apply',t',ps_input);
      %  net = newff(xsn,tn,[20,20]);
        % net = newff(xsn,tn,10);
        % view(net)
        % net.trainParam.epochs = 2000;
        % net.trainParam.goal = 1e-8;
        % net.trainParam.lr = 0.1;
        % net = train(net,xsn,tn);
     %   ffnn=ffnetwork({net.IW{1},net.LW{2,1},net.LW{3,2}},net.b',{'ReLu','ReLu','purelin'});
    
     %ffnn=ffnetwork({net.IW{1},net.LW{2,1}},net.b',{'ReLu','purelin'});

%% Maybe we try ELM
ELMs1=ELM.GenerateELM(size(xs,2),40,'ReLu',size(t,2));
%ELMs1=trainELMLipridge(ELMs1,xsn,tn);
ELMs1=trainELM(ELMs1,xsn,tn);
ffnn=ffnetwork(ELMs1.weight,ELMs1.bias,{'ReLu','purelin'});
%% Initialization
% How many trajectories we want to obtain
% Sampling data from working zone

switch(n)
    case 1
        inputbound =[-44.35;-44.3;-3.12;-3.1];
        thetaInterval=[-46 -44;-3.1 -2];
    case 2
         thetaInterval=[-5 5;25 42];
         inputbound = [-1;0;28;45];
    case 3
          thetaInterval=[-31.5 -31;41.5 43];
          inputbound =[-31.2;-31;42.8;42.9];
    case 4
        inputbound =[-17;-19;-20;-16];
        thetaInterval=[-23 -18;-20 -16];
    case 5
        inputbound =[-47.02;-47;0;0.2];
        thetaInterval=[-48 -40;0 1];
end
%% Load Data Set

% Find the bounderies
%lowerbound= min(tn')-0.1 ;
%upperbound= max(tn')+0.1;
lowerbound= min(xsn')-0.1;
upperbound= max(xsn')+0.1;
init_interval{1}=[lowerbound',upperbound'];

%% Data-driven Partitioning

P=partitions(init_interval,xsn',tn');
intervals=ME(P,tol,maximum_entropy,Dimension);
P.intervals=intervals;
figure
    partitions.intervalplot(intervals,'empty','b')
grid on

% %% Training of the neural networks
% ELMs1=ELM.GenerateELM(size(xs,2),NeuronNum_switch,'ReLu',size(t,2));
% [P1,ELMs]=MergePatitions(P,ELMs1,e);
% 
% % Plot intervals
% figure
% partitions.intervalplot(P.intervals,'empty','black')
% figure
% partitions.intervalplot(P.intervals,'empty','black')
% partitions.intervalplot(P1.intervals,'full','red')
% %title('Invariantspace using Bisection method')

%% In this run we only apply BP as the reference neural network.
%% Verfify whether can it approximate the dynamics well
% 1. Reverse input bound
thetaInterval_lowerboundn = mapminmax('apply',thetaInterval(:,1),ps_input);
thetaInterval_upperboundn = mapminmax('apply',thetaInterval(:,2),ps_input);
thetaIntervaln=[thetaInterval_lowerboundn,thetaInterval_upperboundn];
TestTraj = 10;%number of test trajectories
%Randomly generate staring state
%thetaInterval=[floor(min(xs(:,1))) ceil(max(xs(:,1)));floor(min(xs(:,2))) ceil(max(xs(:,2)))];
RandomStateInput = intervalCompute.randomPoint(thetaIntervaln,TestTraj);
%Random start state should be within the select input space
j=1;
segmentIndex=P.intervals;
inputspace1=P.intervals;
for i=1:TestTraj
    for k = 1:size(segmentIndex,2)
       % for z = 1:size(segmentIndex{k},1)/2
            %  if(RandomStateInput(1,i)>segmentIndex{k}(2*z-1,1))&&(RandomStateInput(1,i)<segmentIndex{k}(2*z-1,2))
                    %if(RandomStateInput(2,i)>segmentIndex{k}(2*z,1))&&(RandomStateInput(2,i)<segmentIndex{k}(2*z,2))
                       Traj{1,j}=RandomStateInput(:,i);
                   % end
              %end
        %end
    end
    j=j+1;
end
validNum=size(Traj,2);
%Generate the system response and break while it's out of input space

% 1.Input uk from uniform distribution

for j=1:size(Traj,2)
  %  uk = zeros(size(2:1:duration,2)+1,1);
  % uk = u*ones(size(uk,1),1)-1+(1-(-1))*rand(size(uk,1),1);
    for i = 2:1:duration
         flag=0;
         %Traj{1,j}(1:2,i)= net(Traj{1,j}(:,i-1));
         Traj{1,j}(1:2,i)= ELMpredict(ELMs1,Traj{1,j}(:,i-1));
         TrajflagDDM{j}(i-1)=k; 
         flag=1;

             if(k==size(inputspace1,2))&&(flag==0) 
                  TrajEndflag(j)=1; 
                  j=j+1;
                  i=2;
              break;
             end
    end
end

for i =1:size(intervals,2)
    for j=1:size(intervals{i},2)
        intervals_reverse_mapminmax{i}(:,j)=mapminmax('reverse',intervals{i}(:,j),ps_input);
    end
end

figure
for i=1:validNum
    Trajn= mapminmax('reverse',Traj{1,i},ps_input);
    plot(Trajn(1,:),Trajn(2,:))
    hold on
end
hold on
partitions.intervalplot(intervals_reverse_mapminmax,'empty','k')

%title('DataDriven State Denpendent Swtiching Model')
%xlabel('\theta_1')
%ylabel('\theta_2')
grid on

%% Relationship between different partitions

% test
CalSetOutput(ffnn,[-100,100;-100,100])

IntersectSet= zeros(size(intervals,2));
for i = 1: size(intervals,2)
    for j = 1: size(intervals,2)
        I=intervals{j};
        IntersectSet(i,j)=ifNextIntersect(ffnn,intervals{i},intervals{j});
    end
    %  hold on
    %NextSetInput=[];
    %NextSetOutput=[];
end

save(strcat(['Results/',names{n}],'.mat'))
%% Generate Abstraction Graph
end
% 1.Partition graph
figure
    partitions.intervalplot(intervals,'empty','b')
grid on

%%
% figure
%     partitions.intervalplot(NextSetInput,'empty','b')
% grid on

Num=ceil(sqrt(size(IntersectSet,1)));
% 1.Compute Each position in Graph
relativeRate=2400/size(IntersectSet,1);
k=1;
for i = 1:Num
    for j =1:Num
        Location{k}=strcat(['x="',num2str(i*relativeRate+5)],['" y="',num2str(j*relativeRate+5),'"']);
        point{k}{1}=[num2str(i*relativeRate+5)];
        point{k}{2}=[num2str(j*relativeRate+5)];
        k=k+1;
    end
end
% 2.Add node name
for i = 1:size(IntersectSet,1)
    b=strcat(['location id= id',num2str(i-1)],[' ',Location{i}]);
    LocationInfo{i}= replace(b,"''",'"');
    points{i}=point{i};
    %nodename{i}=strcat(['id',num2str(i-1)]);
    nodename{i}=strcat('[',num2str(P.intervals{i}(1,1)),',',num2str(P.intervals{i}(1,2)),']');
    %nodename{i}=LocationInfo{i};
end
% 2.
%ModelGraph=digraph(IntersectSet,nodename);
ModelGraph=digraph(IntersectSet);
figure
plot(ModelGraph,'r')
