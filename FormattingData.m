%% 
Data = MacroPhageData; % Get data from Excel

[ L , W ] = size(Data); 

D = [];

y = 1;
for x = 1:L
    
    if( isnan(Data(x,1)) == 0)
        
        D(y,:) =   Data(x,:);
        y = y+1;
        
    end
    
end

q = permute(reshape(D, 26, 27, 14), [1, 3, 2]); % Leaves data as 26 time points and avg of 27 Cytokines with width experiment 

D = q;
%
%only need 3,6,9,12,13 columns

D = D(:,[3,6,9,12,13],:);
D = D(1:25,:,:);

%%%% Devide into the 4 types
Mo = D(1:7,:,:);
IL4 = D([1,8,9,10,11,12,13],:,:);
IL10 = D([1,14,15,16,17,18,19],:,:);
LPS = D([1,20,21,22,23,24,25],:,:);

%%% Convert into columns for fitting 
Mo = reshape(Mo,[35,27]);
IL4 = reshape(IL4,[35,27]);
IL10 = reshape(IL10,[35,27]);
LPS = reshape(LPS,[35,27]);

t = [ 0 .5 1 3 6 12 24 0 .5 1 3 6 12 24 0 .5 1 3 6 12 24 0 .5 1 3 6 12 24];

for x = 1:27 %% Fit for each cytokine
    %%%%% holding values 
    m = Mo(1:28,x);
    i4 = IL4(1:28,x);
    i10 = IL10(1:28,x);
    lps = LPS(1:28,x);
    
    %%%%% each get their own time
    tm = t;
    t4 = t;
    t10 = t;
    tlps = t;
    for y=1:28 %% have to process in case reading doesnt exist 
        if (isnan(m(y)) == 1)
            tm(y) = []; 
            m(isnan(m)) = [];
        end
        
        if(isnan(i4(y)) == 1) 
            t4(y) = []; 
            i4(isnan(i4)) = [];
        end
        
        if(isnan(i10(y)) == 1) 
            t10(y) = []; 
            i10(isnan(i10)) = [];
        end
        
        if(isnan(lps(y)) == 1) 
            tlps(y) = [];   
            lps(isnan(lps)) = [];
        end
    end
    [A,B,C] = fit(tm',m,'smoothingspline');
     M(x).func = A;
     M(x).err = B;
     M(x).par = C;
    [A,B,C] = fit(t4',i4,'smoothingspline');
     I4(x).func = A;
     I4(x).err = B;
     I4(x).par = C;
    [A,B,C] = fit(t10',i10,'smoothingspline');
     I10(x).func = A;
     I10(x).err = B;
     I10(x).par = C;
    [A,B,C] = fit(tlps',lps,'smoothingspline');
     LP(x).func = A;
     LP(x).err = B;
     LP(x).par = C;
end

%%
%%%% getting errors
MavgErr = zeros(5,1);
I4avgErr = zeros(5,1);
I10avgErr = zeros(5,1);
LPSavgErr = zeros(5,1);
for x = 1:27
    
    MavgErr(1) =  MavgErr(1) + M(x).err.sse;
    MavgErr(2) =  MavgErr(2) + M(x).err.rsquare;
    MavgErr(3) =  MavgErr(3) + M(x).err.dfe;
    MavgErr(4) =  MavgErr(4) + M(x).err.adjrsquare;
    MavgErr(5) =  MavgErr(5) + M(x).err.rmse;
    
    I4avgErr(1) =  I4avgErr(1) + I4(x).err.sse;
    I4avgErr(2) =  I4avgErr(2) + I4(x).err.rsquare;
    I4avgErr(3) =  I4avgErr(3) + I4(x).err.dfe;
    I4avgErr(4) =  I4avgErr(4) + I4(x).err.adjrsquare;
    I4avgErr(5) =  I4avgErr(5) + I4(x).err.rmse;
    
    I10avgErr(1) =  I10avgErr(1) + I10(x).err.sse;
    I10avgErr(2) =  I10avgErr(2) + I10(x).err.rsquare;
    I10avgErr(3) =  I10avgErr(3) + I10(x).err.dfe;
    I10avgErr(4) =  I10avgErr(4) + I10(x).err.adjrsquare;
    I10avgErr(5) =  I10avgErr(5) + I10(x).err.rmse;
    
    LPSavgErr(1) =  LPSavgErr(1) + LP(x).err.sse;
    LPSavgErr(2) =  LPSavgErr(2) + LP(x).err.rsquare;
    LPSavgErr(3) =  LPSavgErr(3) + LP(x).err.dfe;
    LPSavgErr(4) =  LPSavgErr(4) + LP(x).err.adjrsquare;
    LPSavgErr(5) =  LPSavgErr(5) + LP(x).err.rmse;
end

MavgErr  = MavgErr./27;
I4avgErr =  I4avgErr./27;
I10avgErr =  I10avgErr./27;
LPSavgErr =  LPSavgErr./27;

%%
%%%% fiiting for each sample 


D = D(:,[3,6,9,12,13],:);
D = D(1:25,:,:);

%%%% Devide into the 4 types
Mo = D(1:7,:,:);
IL4 = D([1,8,9,10,11,12,13],:,:);
IL10 = D([1,14,15,16,17,18,19],:,:);
LPS = D([1,20,21,22,23,24,25],:,:);
t = [0 .5 1 3 6 12 24];

for x = 1:27
    % Each cytokine with all 4 experiments and no average 
    M = Mo(:,1:4,x);
    I4 = IL4(:,1:4,x);
    I10 = IL10(:,1:4,x);
    L = LPS(:,1:4,x);
    for y=1:4
    m = M(:,y);
    i4 = I4(:,y);
    i10 = I10(:,y);
    l = L(:,y);
    if( y == 3)|| (y==1)
    tm = t;
    t4 = t;
    t10 = t;
    tlps = t;   
    m = M(:,y);
    i4 = I4(:,y);
    i10 = I10(:,y);
    l = L(:,y);
    else
    tm = t;
    t4 = t;
    t10 = [0 .5 3 6 12 24];
    tlps = t;    m = M(:,y);
    i4 = I4(:,y);
    i10 = I10([1 2 4 5 6 7],y);
    l = L(:,y);
    end
   
    
            
    m(isnan(i10)) = [];
    i4(isnan(i10)) = [];
    i10(isnan(i10)) = [];
    l(isnan(i10)) = [];
    
    
    [A,B,C] = fit(tm',m,'smoothingspline')  ;      
    Mof(x,y).func = A;
    Mof(x,y).err = B;
    Mof(x,y).para = C;
    
    [A,B,C] = fit(t4',i4,'smoothingspline') ;  
    IL4f(x,y).func = A;
    IL4f(x,y).err = B;
    IL4f(x,y).para = C; 
    [A,B,C] = fit(t10',i10,'smoothingspline')  ;
    IL10f(x,y).func = A;
    IL10f(x,y).err = B;
    IL10f(x,y).para = C;
    
    [A,B,C] = fit(tlps',l,'smoothingspline') ;
    LPSf(x,y).func =A;
    LPSf(x,y).err = B;
    LPSf(x,y).para = C;
    end
end
%% error
%%%% getting errors
MavgErr = zeros(5,1);
I4avgErr = zeros(5,1);
I10avgErr = zeros(5,1);
LPSavgErr = zeros(5,1);
for x = 1:27
    for y = 1:4
    
    MavgErr(1) =  MavgErr(1) + Mof(x,y).err.sse;
    
    MavgErr(3) =  MavgErr(3) + Mof(x,y).err.dfe;
    
    MavgErr(5) =  MavgErr(5) + Mof(x,y).err.rmse;
    
    I4avgErr(1) =  I4avgErr(1) + IL4f(x,y).err.sse;
    
    I4avgErr(3) =  I4avgErr(3) + IL4f(x,y).err.dfe;
    
    I4avgErr(5) =  I4avgErr(5) + IL4f(x,y).err.rmse;
    
    I10avgErr(1) =  I10avgErr(1) + IL10f(x,y).err.sse;
   
    I10avgErr(3) =  I10avgErr(3) + IL10f(x,y).err.dfe;
   
    I10avgErr(5) =  I10avgErr(5) + IL10f(x,y).err.rmse;
    
    LPSavgErr(1) =  LPSavgErr(1) + LPSf(x,y).err.sse;
    LPSavgErr(2) =  LPSavgErr(2) + LPSf(x,y).err.rsquare;
    LPSavgErr(3) =  LPSavgErr(3) + LPSf(x,y).err.dfe;
    LPSavgErr(4) =  LPSavgErr(4) + LPSf(x,y).err.adjrsquare;
    LPSavgErr(5) =  LPSavgErr(5) + LPSf(x,y).err.rmse;
    end
end
MavgErr  = MavgErr./(27*4);
I4avgErr =  I4avgErr./(27*4);
I10avgErr =  I10avgErr./(27*4);
LPSavgErr =  LPSavgErr./(27*4);   
% R had nan so had to be calculated seperate
countm = 0;
count4 = 0;
count10 = 0;
for x = 1:27
    for y = 1:4
    
    if (isnan(Mof(x,y).err.rsquare) == 0);
    MavgErr(2) =  MavgErr(2) + Mof(x,y).err.rsquare;
   
    MavgErr(4) =  MavgErr(4) + Mof(x,y).err.adjrsquare;
    countm = countm + 1;
    end
    
    if(isnan(IL4f(x,y).err.rsquare) == 0);
    I4avgErr(2) =  I4avgErr(2) + IL4f(x,y).err.rsquare;
   
    I4avgErr(4) =  I4avgErr(4) + IL4f(x,y).err.adjrsquare;
    count4 = count4 + 1;
    end
    
    if(isnan(IL10f(x,y).err.rsquare) == 0);
    I10avgErr(2) =  I10avgErr(2) + IL10f(x,y).err.rsquare;
   
    I10avgErr(4) =  I10avgErr(4) + IL10f(x,y).err.adjrsquare;
    count10 = count10 + 1;
    end
    
   
    end
end
I10avgErr(2) = I10avgErr(2)/count10;
I10avgErr(4) = I10avgErr(4)/count10;

I4avgErr(2) = I4avgErr(2)/count10;
I4avgErr(4) = I4avgErr(4)/count10;

MavgErr(2) =  MavgErr(2)/countm;
MavgErr(4) =  MavgErr(4)/countm;
%%
%%%Sampling 
Mo = zeros(50,4,27);
IL4 = zeros(50,4,27);
IL10 = zeros(50,4,27);
LPS = zeros(50,4,27);
a = 0:(9/39):9;
b = (9+1.5):(15/10):24;
t = [a,b];
for x = 1:27
    for y = 1:4
    
       Mo(:,y,x) = Mof(x,y).func(t);
       IL4(:,y,x) = IL4f(x,y).func(t);
       IL10(:,y,x) = IL10f(x,y).func(t);
       LPS(:,y,x) = LPSf(x,y).func(t);
        
    end
end
% Avg sample 
Mo = zeros(50,27);
IL4 = zeros(50,27);
IL10 = zeros(50,27);
LPS = zeros(50,27);
a = 0:(9/39):9;
b = (9+1.5):(15/10):24;
t = [a,b];

for x = 1:27
    
    
       Mo(:,x) = M(x).func(t);
       IL4(:,x) = I4(x).func(t);
       IL10(:,x) = I10(x).func(t);
       LPS(:,x) = LP(x).func(t);
        
   
end


% plot for fitting
subplot(3,1,1);
plot([0 .5 1 3 6 12 24], IL10(:,3,21),'*')
set(gca,'fontsize',16)
hold on
plot(IL10f(21,3).func);
title('2.A')
xlabel('Time','FontSize',16)
ylabel('Expression Level','FontSize',16)

subplot(3,1,2);
plot([0 .5 1 3 6 12 24], LPS(:,1,20),'*')
set(gca,'fontsize',16)
hold on
plot(LPSf(20,1).func);
title('2.B')
xlabel('Time','FontSize',16)
ylabel('Expression Level','FontSize',16)

subplot(3,1,3);
plot([0 .5 1 3 6 12 24], IL4(:,3,5),'*')
set(gca,'fontsize',16)
hold on
plot(IL4f(5,3).func);
title('2.C')
xlabel('Time','FontSize',16)
ylabel('Expression Level','FontSize',16)

%%5 in case it goes below negative 
for x = 1:27
    
    for y = 1:4
        M = Mo(:,y,x);
        I4 = IL4(:,y,x);
        I10 = IL10(:,y,x);
        L = LPS(:,y,x);
        % large number so it wont affect min
        M(M<=0)= [];
        I4(I4<=0) = [];
        I10(I10<=0) = [];
        L(L<=0) = [];
        
        minM = min(M);
        min4 = min(I4);
        min10 = min(I10);
        minL = min(L);
        
        M = Mo(:,y,x);
        I4 = IL4(:,y,x);
        I10 = IL10(:,y,x);
        L = LPS(:,y,x);
        
        M(M<0) = minM;
        I4(I4<0) = min4;
        I10(I10<0) = min10;
        L(L<0) = minL;
        
        Mo(:,y,x) = M;
        IL4(:,y,x) = I4;
        IL10(:,y,x) = I10;
        LPS(:,y,x) = L;
    end
    
end
   % for avg  get rid of neg 
for x = 1:27
    
  
        M = Mo(:,x);
        I4 = IL4(:,x);
        I10 = IL10(:,x);
        L = LPS(:,x);
        % large number so it wont affect min
        M(M<=0)= [];
        I4(I4<=0) = [];
        I10(I10<=0) = [];
        L(L<=0) = [];
        
        minM = min(M);
        min4 = min(I4);
        min10 = min(I10);
        minL = min(L);
        
        M = Mo(:,x);
        I4 = IL4(:,x);
        I10 = IL10(:,x);
        L = LPS(:,x);
        
        M(M<0) = minM;
        I4(I4<0) = min4;
        I10(I10<0) = min10;
        L(L<0) = minL;
        
        Mo(:,x) = M;
        IL4(:,x) = I4;
        IL10(:,x) = I10;
        LPS(:,x) = L;

    
end    
%%
%%% binarization
%%% Avg this is for comparison 
Mob = zeros(50,27);
IL4b = zeros(50,27);
IL10b = zeros(50,27);
LPSb = zeros(50,27);
for x = 1:27

    M = Mo(:,x);
    I4 = IL4(:,x);
    I10 = IL10(:,x);
    L = LPS(:,x);
    
    [A,B] = Bi(M,4);
    Mob(:,x) = A;
    [A,B] = Bi(I4,4);
    IL4b(:,x) = A;
    [A,B] = Bi(I10,4);
    IL10b(:,x) = A;
    [A,B] = Bi(L,4);
    LPSb(:,x) = A;

end

% 4 experiments
Mob = zeros(50,4,27);
IL4b = zeros(50,4,27);
IL10b = zeros(50,4,27);
LPSb = zeros(50,4,27);
MoMean = zeros(4,27);
IL4Mean = zeros(4,27);
IL10Mean = zeros(4,27);
LPSMean = zeros(4,27);
for x = 1:27
 for y = 1:4
    M = Mo(:,y,x);
    I4 = IL4(:,y,x);
    I10 = IL10(:,y,x);
    L = LPS(:,y,x);
    
    [A,B] = Bi(M,4);
    Mob(:,y,x) = A;
    MoMean(y,x) = sum(B)/2;
    [A,B] = Bi(I4,4);
    IL4b(:,y,x) = A;
    IL4Mean(y,x) = sum(B)/2;
    [A,B] = Bi(I10,4);
    IL10b(:,y,x) = A;
    IL10Mean(y,x) = sum(B)/2;
    [A,B] = Bi(L,4);
    LPSb(:,y,x) = A;
    LPSMean(y,x) = sum(B)/2;
 end
end
%add time point to z and t1 to fit nicer
close all
t = [0 .5 1 3 6 12 24];
plot(t,IL4(:,1,10),'*')
set(gca,'fontsize',18)
hold on
a = 0:(9/39):9;
b = (9+1.5):(15/10):24;
t = [a,b];
plot(IL4f(10,1).func)
plot([ 0 24], [2.38 2.38])
yyaxis right
plot(t,IL4b(:,1,10))
%title('Binization vs. continuous Function for IL-10 in M2a')
xlabel('Time')
ylabel('Gene Expression Level')
legend('Orginal Time Points','Binarization','Curve Fitting Function','Threshold')
%title('Binarization vs. continuous Function for IL-10 in M2a')

%%
%LOOCV
% muactive/muinactive/sdactive

muactive = zeros(27,4);
muinactive = zeros(27,4);

sdactive = zeros(27,4);
sdinactive = zeros(27,4);
Data = IL4;
Mean = IL4mean;
for x = 1:27
   
    for y = 1:4
       active = [];
       inactive = [];
       a = Data(x,y,:);
       b = Mean(x,y);
       
        for z = 1:50
           
           if(a(z) >= b)
           
               active = [active;a(z)];
               
           else
               
               inactive = [inactive;a(z)];
               
           end
        muactive(x,y) = mean(active);
        muinactive(x,y) = mean(inactive);
        sdactive(x,y) = std(active);
        sdinactive(x,y) = std(inactive);
              
            
        end
        
    end
    
end
