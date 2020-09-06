%% MONTE CARLO INVERSION OF R WAVE DISPERSION (JOINT PHASE AND GROUP)
% REQUIRES HERRMANN'S CODES TO BE INSTALLED IN COMPUTER
%% Clean Up and Begin!
%tic %keep track of run time
clc; close all; clearvars; format long g; %warning off
cd bin; path_bin=pwd; cd .. %add bin to path
addpath(path_bin); %add bin to path
disp('Surface Wave Dispersion Inversion')
disp(' ')

%tic
%% ITERATIONS: Number of iterations to run.
IT=100000;

%% Models
%Load initial model
Model0=load('model0.5');
sM0=size(Model0,1);
Depths=cumsum(Model0(:,1));

%% LOAD OBSERVED CURVES
% Load the observed dispersion curves for all the waves.
R_Cobs= load('R_C.dat'); % Rayleigh Phase Curve
R_Uobs= load('R_U.dat'); %Rayleigh Group Curve

%% INITIATE VARIABLES
cc=1;
Total_errors=zeros(1,IT);
Accepted=zeros(1,IT);
P=zeros(1,IT);
LM_Total=zeros(1,IT);
CEROS=zeros(size(Model0,1));


%% Periods
% Periods for the computation of the Dispersion curves. Right now, is all
% the range from 1 s to 50 s. These are read from the file
T=(R_Cobs(:,1))';

%% Vp/Vs
% Vp/Vs ratio for the P wave velocity needed for the calculation. The
% inversion is based mainly in the Vs model
VpVs=1.7531;

%% MAX and MIN Vs

maxVstop=4.0; % Max Vs to accept in the cadidate models
minVstop=2.0; % Max Vs to accept in the cadidate models

maxVsbottom=5; % Max Vs to accept in the cadidate models
minVsbottom=3; % Max Vs to accept in the cadidate models

%% Compute Dispersion for Initial Model
% Compute the dispersion curve of the initial model to estimate initial
% RMSE
Disp0 = dispR_surf96(T,Model0,VpVs);
% Order the information 
R_C= [T' Disp0(1:size(T,2)) ]; % Rayleigh Phase Curve
R_U= [T' Disp0((size(T,2)+1):end)]; %Rayleigh Group Curve

%Size of the vector
tam=size(R_U(:,2),1);

%% INITIAL ROOT MEAN SQUARE ERROR

R_C_error0=sqrt(sum((R_Cobs(:,2)-R_C(:,2)).^2)/tam);
R_U_error0=sqrt(sum((R_Uobs(:,2)-R_U(:,2)).^2)/tam);

% Initial RMSE
Total_errors(1)=R_C_error0+R_U_error0;
% FIRST MODEL = Initial Model
MODELS(:,:,1)=Model0;

%% INITIAL VEROSIMILITY & PROBABILITY

fv=0.022; %Variance for the inversion (Sigma square)
i2fv=1/(2*fv); %to save time
fv2=0.01; %Variance for the inversion (Sigma square)
i2fv2=1/(2*fv2); %to save time

% Vesosimilitud de cara parametro
LMi_R_C=exp( -sum( ((R_C(:,2)-R_Cobs(:,2)).^2) *i2fv ) );
LMi_R_U=exp( -sum( ((R_U(:,2)-R_Uobs(:,2)).^2) *i2fv ) );

LM_Total(1)=LMi_R_C*LMi_R_U;

% Primera probabilidad es 1
P(1)=1;

%% CONSTANTS
% Choose a random central layer to perturb.
% FIRST LAYER DOESNT COUNT
a = 2; %First layer in the range to perturb
b = 11; %Last layer to in the range to perturb
sl=5; %Step layer
%Fl=21; %Final layer. This layer and the 2 subsecuent layers move to change keep the rest of the model consistet
% Random amplitude of the assymetric perturbation
c = 0; %Pertub at least 1 layer: Could be 0 to allow for stronger variations
d = 2; %Perturb Max 
smoothfactor=0.25; %smooth factor for the code. Depends on the case %0.25 funciona bien
%d=5
for i=2:IT
%% Perturbab the Model

%% Choose what to do
Choose=rand(1);

%% If Choose is >= 0.5 perturb Vs model
if Choose >=0.5 
%% Test Gaussianas
% Choose a random central layer to perturb.
r1 = round((b-a).*rand(1) + a); %Number of the random layer to be perturb
% Random amplitude of the assymetric perturbation
e1=(d-c).*rand(1) + c; %How large is the perturbation on one side
e2=(d-c).*rand(1) + c; %How large is the perturbation of the other
%Random Amplitud of the perturbation 
h = (1--1).*rand(1) + -1;
%Gaussuan Perturnation
G = gauss2mf(1:sM0,[e1 r1 e2 r1])*h;
%Only perturb from the top to the last layer allowed
G(b+1:end)=0;
% NEW PERTURB MODEL!
v=MODELS(:,2,i-1)+G';
%% Check the model

if min(v(a:sl))<=minVstop+0.05 || max(v(a:sl))>=maxVstop-0.05
    G=G*0.25;
    v=MODELS(:,2,i-1)+G';
end

if min(v(sl:b))<=minVsbottom+0.05 || max(v(sl:b))>=maxVsbottom-0.05
    G=G*0.25;
    v=MODELS(:,2,i-1)+G';
end

%v1=smooth(v(a:b),smoothfactor,'sgolay');
v1=smooth(v(a:b),smoothfactor);
v2=v(b+1:end);
v=[v(2); v1 ; v2];
Thickness=MODELS(:,1,i-1);
%%
% If we are going to perturb a depth

else
r1 = round((b-a).*rand(1) + a); %Number of the random layer to be perturb
%Get the thickness of the layer
Thickness=MODELS(:,1,i-1);
v=MODELS(:,2,i-1);
%Get the random value to perturb the layer
dz = (2.5--2.5).*rand(1) + -2.5;
%Add the value to the layer
Thickness(r1)=Thickness(r1)+dz;
% If too much or too little go back
if min(Thickness(a:b))<1 || max(Thickness(a:b))>10
    Thickness=MODELS(:,1,i-1);
else
% If all is OK The last 3 layers after the value move to accomodate the model
Thickness(b+1)=Thickness(b+1)-(dz/3);
Thickness(b+2)=Thickness(b+2)-(dz/3);
Thickness(b+3)=Thickness(b+3)-(dz/3);
end
end

% Final fix
v(1)=v(2);
Thickness(1)=0;

ModelTest=[Thickness v];

%% Compute Dispersion for Perturb Model
try
DispTest = dispR_surf96(T,ModelTest,VpVs);

R_C_test= [T' DispTest(1:size(T,2)) ]; % Love Phase Curve
R_U_test= [T' DispTest(size(T,2)+1:end)]; %Love Group Curve


%% ROOT MEAN SQUARE ERROR

R_C_error=sqrt(sum((R_C_test(:,2)-R_Cobs(:,2)).^2)/tam);
R_U_error=sqrt(sum((R_U_test(:,2)-R_Uobs(:,2)).^2)/tam);

% Total errors
errori=R_C_error+R_U_error;

%% VEROSIMILITY (L)

if IT<IT/2 % Use origina Sigma
    LM_R_C=exp( -sum( ((R_C_test(:,2)-R_Cobs(:,2)).^2) *i2fv ) );
    LM_R_U=exp( -sum( ((R_U_test(:,2)-R_Uobs(:,2)).^2) *i2fv ) );
else %Reduce Sigma half way 
    LM_R_C=exp( -sum( ((R_C_test(:,2)-R_Cobs(:,2)).^2) *i2fv2 ) );  
    LM_R_U=exp( -sum( ((R_U_test(:,2)-R_Uobs(:,2)).^2) *i2fv2 ) );
end
 
%Total Likelyhood
LM_Totaltest=LM_R_C*LM_R_U;

%% PROBABILITY (P)

Ptest=min(1,LM_Totaltest/LM_Total(i-1));

%% RANDOM VALUE

RR=rand(1);

%% TEST the Model
if RR <= Ptest % ACCEPT THE MODEL
    % SAVE MODEL
    MODELS(:,:,i)=ModelTest;
    % SAVE PROBABILITY
    P(i)=Ptest;
    % SAVE VEROSIMILITUD
    LM_Total(i)= LM_Totaltest;
    % SAVE ERROR
    Total_errors(i)=errori;
    
    %Keep Track of Acepted models
    Accepted(cc)=i;
    cc=cc+1;
    
    if Choose >=0.5
    X = ['YES! A Perturbation of Vs was accepted at iteration ',num2str(i), ' with RMSE ', num2str(errori)];
    else
    X = ['YES! A Perturbation of H was accepted at iteration ',num2str(i), ' with RMSE ', num2str(errori)];  
    end
    disp(X)
else % REJECT MODEL
    % SAVE PREVIOUS MODEL
    MODELS(:,:,i)=MODELS(:,:,i-1);
    % SAVE PREVIOUS  PROBABILITY
    P(i)=P(i-1);
    % SAVE PREVIOUS VEROSIMILITUD
    LM_Total(i)=LM_Total(i-1);
    % SAVE ERROR
    Total_errors(i)=Total_errors(i-1);
end

catch
    disp('SHITTY MODEL GAVE A PROBLEM! FUCK IT!!!')
     % SAVE 2 MODELS AGO, TO SKIP A PROBLEMATIC PREVIOUS MODEL
    MODELS(:,:,i)=MODELS(:,:,i-2);
    % SAVE 2 MODELS AGO  PROBABILITY
    P(i)=P(i-2);
    % SAVE 2 MODELS AGO VEROSIMILITUD
    LM_Total(i)=LM_Total(i-2);
    % SAVE 2 MODELS AGO ERROR
    Total_errors(i)=Total_errors(i-2);  
    % REMOVE CRAP
    system('rm tmp* sobs.d lu lc rc ru disp_obs.dsp NUL temp.dsp ')    
    
end

end

%toc
%% Accepted
Accepted=Accepted(1:cc-1);

%%
load model0.Crust2
%% PLOT RMSE for selection
%%

figure
plot(Total_errors)
xlabel('Iterations')
ylabel('RMSE')
title('Root Mean Square Error')
CUT=ginput(1);
CUT=round(CUT(1));

%% PLOT

figure
subplot(3,2,1)
plot(1:IT, Total_errors)
xlabel('Iterations')
ylabel('RMSE')
subplot(3,2,3)
plot(1:IT, P)
xlabel('Iterations')
ylabel('Probability')
subplot(3,2,5)
plot(1:IT, LM_Total)
xlabel('Iterations')
ylabel('Likelyhood')
subplot(3,2,[2 4 6])
SAMPLES=round(1:(IT-1)/100:IT);
Vector=[SAMPLES' Total_errors(SAMPLES)'];
colo=jet(size(SAMPLES,2));
for samp=1:size(SAMPLES,2)
   stairs(MODELS(:,2,Vector(samp)), cumsum(MODELS(:,1,Vector(samp))),'color',colo(samp,:))
   hold on
end
stairs(model0(:,2), cumsum(model0(:,1)),'k', 'LineWidth',2)
set(gca,'XDir','normal','YDir','reverse');
ylabel('Depth (km)')
xlabel('Wave Speed (km/s)')
xlim([2 5])
ylim([Depths(1) Depths(b)])
grid on
title('Velocity Model')

%% CALCULATE STATISTICS MODELS 
%%
Mean_Model=[mean(MODELS(:,1,CUT:IT),3) mean(MODELS(:,2,CUT:IT),3)];
Mode_Model=[mode(MODELS(:,1,CUT:IT),3) mode(MODELS(:,2,CUT:IT),3)];
Median_Model=[median(MODELS(:,1,CUT:IT),3) median(MODELS(:,2,CUT:IT),3)];

%% COMPUTE DISPERSIONS OF STATISTICS MODELS
%% I

Disp_Mean = dispR_surf96(T,Mean_Model,VpVs);

R_C_Mean= [T' Disp_Mean(1:size(T,2)) ]; % Love Phase Curve
R_U_Mean= [T' Disp_Mean((size(T,2)+1):end)]; %Love Group Curve

Disp_Mode = dispR_surf96(T,Mode_Model,VpVs);

R_C_Mode= [T' Disp_Mode(1:size(T,2)) ]; % Love Phase Curve
R_U_Mode= [T' Disp_Mode((size(T,2)+1):end)]; %Love Group Curve

Disp_Median = dispR_surf96(T,Median_Model,VpVs);

R_C_Median= [T' Disp_Median(1:size(T,2)) ]; % Love Phase Curve
R_U_Median= [T' Disp_Median((size(T,2)+1):end)]; %Love Group Curve

%% PROBABILITY MATRIX

A=MODELS(1:b,2,CUT:IT);
B=cumsum(MODELS(1:b,1,CUT:IT));
Xedges=1.5:0.05:5;
Yedges=0:0.5:50;
N = histcounts2(A(:),B(:),Xedges,Yedges,'Normalization', 'probability');
[YY,III]=max(N,[],2); 
HPDEPTHS=Yedges(III);
BBB=N(:,III);
[YYY2,III2]=max(BBB,[],1);
HPVs=Xedges(III2);

HPM1=[HPVs',HPDEPTHS'];
sHPM1=sortrows(HPM1,2);
[C,IA,IC]=unique(sHPM1(:,2));
usHPM1=[sHPM1(IA,1) sHPM1(IA,2)];


%% FIGURES OF MODELS
%% I: Velocity models from statistics to the Chain

figure

%Rayleigh
subplot(1,2,1)
plot(R_Cobs(:,1),R_Cobs(:,2),'ks')
hold on
plot(R_Uobs(:,1),R_Uobs(:,2),'ko')
plot(R_C_Mean(:,1),R_C_Mean(:,2),'Color','r')
plot(R_U_Mean(:,1),R_U_Mean(:,2),'--','Color','r')
plot(R_C_Mode(:,1),R_C_Mode(:,2),'Color','b')
plot(R_U_Mode(:,1),R_U_Mode(:,2),'--','Color','b')
plot(R_C_Median(:,1),R_C_Median(:,2),'Color','g')
plot(R_U_Median(:,1),R_U_Median(:,2),'--','Color','g')
plot(R_C(:,1),R_C(:,2),'Color',[125/255 125/255 125/255])
plot(R_U(:,1),R_U(:,2),'--','Color',[125/255 125/255 125/255])
legend('R C Obs','R U Obs','R C Mean','R U Mean','R C Mode','R U Mode','R C Median','R U Median','R C Initial','R U Initial','Location','southeast')
xlabel('Period (s)')
ylabel('Wave Speed (km/s)')
grid on
title('Rayleigh Dispersion Curves')

subplot(1,2,2)

mmm=max(log10(N(:)'));
AAA=mmm./log10(N');
imagesc(Xedges,Yedges,AAA)
title('Probability of Models')
hold on
stairs(Model0(:,2), cumsum(Model0(:,1)),'k', 'LineWidth',1.5)
stairs(model0(:,2), cumsum(model0(:,1)),'Color',[200/255 200/255 200/255], 'LineWidth',1.5)
stairs(Mean_Model(:,2),cumsum(Mean_Model(:,1)),'Color','r','LineWidth',1.5)
stairs(Mode_Model(:,2),cumsum(Mode_Model(:,1)),'Color','b','LineWidth',1.5)
stairs(Median_Model(:,2),cumsum(Median_Model(:,1)),'Color','g','LineWidth',1.5)
stairs(usHPM1(:,1),usHPM1(:,2),'Color','y','LineWidth',1.5)

ylabel('Depth (km)')
xlabel('Wave Speed (km/s)')
set(gca,'XDir','normal','YDir','reverse');
ccc = colorbar;
ccc.Label.String = 'log10(Probability)';
colormap(viridis(10))
grid on
legend('Inital Model','Solution','Mean Model','Mode Model','Median','Location','southwest')



% stairs(Mean_Model(:,2),cumsum(Mean_Model(:,1)),'Color',[210/255 0/255 0/255])
% stairs(Mode_Model(:,2),cumsum(Mode_Model(:,1)),'Color',[0/255 0/255 170/255])
% stairs(Median_Model(:,2),cumsum(Median_Model(:,1)),'Color',[0/255 150/255 100/255])


%  Candidate1=(Mean_Model1(1:b,2) + Mode_Model1(1:b,2) +Median_Model1(1:b,2))./3;
%  Candidate2=(Mean_Model(1:b,2) + Mode_Model(1:b,2) +Median_Model(1:b,2))./3;
%  Candidate3=(Mean_Model1(1:b,2) + Mode_Model1(1:b,2) + Median_Model1(1:b,2) + Mean_Model(1:b,2) + Mode_Model(1:b,2) +Median_Model(1:b,2))./6;
%  Candidate4= smooth(XX(I),0.35,'sgolay');
%  Candidate5=mean([Candidate1 Candidate2 Candidate3 Candidate4],2);

%% OUTPUT CANDIDATE MODELS

%dlmwrite('MC_Inversion_Output.txt',[Depths(a:b) Candidate1 Candidate2 Candidate3 Candidate4 Candidate5],'delimiter','\t','precision',4)

%% Final Clean Up

system('rm disp_obs.dsp sobs.d start.mod temp.dsp');

