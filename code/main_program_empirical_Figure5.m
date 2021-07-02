
clear all
rng('default')
rng(101)
data=xlsread('city_data_9_12_v14.xlsx');
%Wuhan is listed as the 162th city
%column 11=newly confirmed case excluding importation in March
%column 12=cumulative case
%column 13=population flow within city
%column 14=smoothed population flow within city
%column 15=abroad index
%column 16=smoothed abroad index
%column 17=current temperature
%column 18=1 week lagged temperature
%column 19=2 week lagged temperature
%column 20=3 week lagged temperature
%column 21=4 week lagged temperature
%column 22=current min temperature
%column 23=1 week lagged min temperature
%column 24=2 week lagged min temperature
%column 25=3 week lagged min temperature
%column 26=4 week lagged min temperature
%column 27=current max temperature
%column 28=1 week lagged maximum temperature
%column 29=2 week lagged maximum temperature
%column 30=3 week lagged maximum temperature
%column 31=4 week lagged maximum temperature
%cloumn 32=current wind speed
%column 33=1 week lagged wind speed
%column 34=2 week lagged wind speed
%column 35=3 week lagged wind speed
%column 36=4 week lagged wind speed
%column 37=current maximum wind speed
%column 38=1 week lagged maximum wind speed
%column 39=2 week lagged maximum wind speed
%column 40=3 week lagged maximum wind speed
%column 41=4 week lagged maximum wind speed
%column 42=precipitation
%column 43=1 week lagged precipitation
%column 44=2 week lagged precipitation
%column 45=3 week lagged precipitation
%column 46=4 week lagged precipitation
%column 47=total number of hospital
%column 48=total number of sanjia hospital
%column 49=total number of zhuanke hospital
%column 50=total number of chuanran hospital
%column 51=total number of yiliaobaojian hospital
%column 52=total number of weishengyuan
%column 53=total number of kouqiang
%column 54=total number of fuke hospital
%column 55=total number of jijiu hospital
%column 56=total number of zheng xing hospital
%column 57=total number of jibingfangyu hospital
%column 58=total number of yanke
%column 59=total number of jingshen
%column 60=total number of zonghe
%column 61=total number of erbihou
%column 62=total number of zhongliu
%column 63=total number of xiongke
%column 64=total number of Naoke
%column 65=total number of zhenshuo
%column 66=total number of guke
%column 67=population density
%column 68=GDP per capita
%column 70-96 provincial effect
%column 97  5 days lag abroad index I (not used)
%column 98 smoothed abroad index I for 2020 (not used)
%column 99 smoothed move-in index 2020
%column 100 smoothed move-in index 2019
%column 101 smoothed abroad index normalized by city 2020
%column 102 smoothed abroad index normalized by national 2020
%column 103 smoothed abroad index without normalization 2020
%column 104 log smoothed abroad index without normalization 2020 
%column 105 smoothed abroad index normalized by city 2019
%column 106 smoothed abroad index normalized by national 2019
%column 107 smoothed abroad index without normalization 2019
%column 108 log smoothed abroad index without normalization 2019 
%column 109 7 day average log smoothed abroad index 2020
%column 110 14 day lagged move-in index for 2020
%column 111 7 day lagged move-in index for 2019
%column 112 14 day lagged move-in index for 2019
%column 113 log smoothed abroad index without normalization 2020 (5 days)
%column 114 7 day log smoothed abroad index without normalization 2020 
%column 115 14 day log smoothed abroad index without normalization 2020 
%column 116 log smoothed abroad index without normalization 2019 (5 days)
%column 117 7 day log smoothed abroad index without normalization 2019 
%column 118 14 day log smoothed abroad index without normalization 2019 
%column 119 14 day average log smoothed abroad index 2020, t-14 to t-1
%column 120 9 day average log smoothed abroad index 2020, t-9 to t-1
%column 121 average log smoothed abroad index from t-8 to t-14
%column 122 importation cases for different cities
%column 123 total newly confirm cases for cities including road
%importation
%column 124 old total newly confirm cases for cities
%column 125 domestic newly confirm cases for cities
%column 126 air importation cases for different cities
%column 127 total newly confirmed cases excluding road importation
%column 128 is heilongjiang indicator
%column 129 is new cumulative case
%column 130 is indicator for cities with quarantine facilities

yo=data(301:30000,127); %daily confirm cases for all cities,begins at Jan 20th, the 6th date
yl=data(1:29700,127); % lagged daily confirm cases
ycum=data(301:30000,129); % cumulative begins at Jan 21st
X13=data(301:30000,13); %5 days lagged population flow within city
X14=data(301:30000,14); % smoothed population flow within city
bo1=data(301:30000,104); %log abroad index for all cities
bo2=data(301:30000,108); % smoothed abroad index for all cities
X99=data(301:30000,99); % 2020 5 days lag inflow index for all cities
X100=data(301:30000,100); % 2020 smoothed inflow index for all cities

X17=data(301:30000,17); % current temperature 
X18=data(301:30000,18); X19=data(301:30000,19);  X20=data(301:30000,20);    X21=data(301:30000,21);    % 1 week to 4 week lagged temperature
X22=data(301:30000,22); % current min temperature
X23=data(301:30000,23); X24=data(301:30000,24); X25=data(301:30000,25); X26=data(301:30000,26); % 1 week to 4 week min temperature
X27=data(301:30000,27); % current max temperature
X28=data(301:30000,28); X29=data(301:30000,29); X30=data(301:30000,30); X31=data(301:30000,31); % 1 week to 4 week max temperature
X32=data(301:30000,32); % current wind speed
X33=data(301:30000,33); X34=data(301:30000,34); X35=data(301:30000,35); X36=data(301:30000,36); % 1 week to 4 week wind speed
X37=data(301:30000,37); % current maximum wind speed
X38=data(301:30000,38); X39=data(301:30000,39); X40=data(301:30000,40); X41=data(301:30000,41); % 1 week to 4 week max wind speed
X42=data(301:30000,42); % current precipitation
X43=data(301:30000,43); X44=data(301:30000,44); X45=data(301:30000,45); X46=data(301:30000,46); % 1 week to 4 week precipitation
% number of hospitals
X47=data(301:30000,47); X48=data(301:30000,48); X49=data(301:30000,49); X50=data(301:30000,50); X51=data(301:30000,51); X52=data(301:30000,52); X53=data(301:30000,53); X54=data(301:30000,54); X55=data(301:30000,55);
X54=data(301:30000,54); X55=data(301:30000,55); X56=data(301:30000,56); X57=data(301:30000,57); X58=data(301:30000,58); X59=data(301:30000,59); X60=data(301:30000,60); X61=data(301:30000,61); X62=data(301:30000,62);
X63=data(301:30000,63); 

X64=data(301:30000,64);  X65=data(301:30000,65); X66=data(301:30000,66); X67=data(301:30000,67); X68=data(301:30000,68);
Xp=data(1:300,70:96); % provincial effect

xx=data(301:600,67);
[r1,c]=find(isnan(xx));

n=300; T=99; 

r2=[r1;30;65];% delete hulunbeier and mudanjiang

 n1=300-length(r2);
r=[r2;162];

Y=zeros(n1-1,T);  Yw=zeros(T,1); YL=zeros(n1-1,T); X1=zeros(n1-1,30,T);   % X1=zeros(n-1,20,T);  
Z1=zeros(n1-1,48); % exclude the total number of hospital
Xw=zeros(T,31);  YwL=zeros(T,1); index=zeros(n1-1,T); indexc=zeros(n1-1,T);
Ycum=zeros(n1-1,T);




for q=1:T
 yy=yo((q-1)*n+1:q*n,1);
 Yw(q,1)=yy(162); yy(r)=[]; %eliminate Wuhan
 Y(:,q)=yy;
 
 yyl=yl((q-1)*n+1:q*n,1); 
 YwL(q,1)=yyl(162); yyl(r)=[]; YL(:,q)=yyl;
 
  
  yyc=ycum((q-1)*n+1:q*n,1);
 yyc(r)=[];  Ycum(:,q)=yyc;
 

  x17=X17((q-1)*n+1:q*n,1); Xw(q,1)=x17(162); x17(r)=[]; X1(:,1,q)=x17;  x18=X18((q-1)*n+1:q*n,1); Xw(q,2)=x18(162); x18(r)=[]; X1(:,2,q)=x18;  
  x19=X19((q-1)*n+1:q*n,1); Xw(q,3)=x19(162); x19(r)=[]; X1(:,3,q)=x19; x20=X20((q-1)*n+1:q*n,1); Xw(q,4)=x20(162); x20(r)=[]; X1(:,4,q)=x20; 
  x21=X21((q-1)*n+1:q*n,1); Xw(q,5)=x21(162); x21(r)=[]; X1(:,5,q)=x21;
  
  x22=X22((q-1)*n+1:q*n,1); Xw(q,6)=x22(162); x22(r)=[]; X1(:,6,q)=x22; 
   x23=X23((q-1)*n+1:q*n,1); Xw(q,7)=x23(162); x23(r)=[]; X1(:,7,q)=x23; 
   x24=X24((q-1)*n+1:q*n,1); Xw(q,8)=x24(162); x24(r)=[]; X1(:,8,q)=x24; 
  x25=X25((q-1)*n+1:q*n,1); Xw(q,9)=x25(162); x25(r)=[]; X1(:,9,q)=x25; 
  x26=X26((q-1)*n+1:q*n,1); Xw(q,10)=x26(162); x26(r)=[]; X1(:,10,q)=x26; 
   x27=X27((q-1)*n+1:q*n,1); Xw(q,11)=x27(162); x27(r)=[]; X1(:,11,q)=x27; 
   x28=X28((q-1)*n+1:q*n,1); Xw(q,12)=x28(162); x28(r)=[]; X1(:,12,q)=x28; 
   x29=X29((q-1)*n+1:q*n,1); Xw(q,13)=x29(162); x29(r)=[]; X1(:,13,q)=x29;
   x30=X30((q-1)*n+1:q*n,1); Xw(q,14)=x30(162); x30(r)=[]; X1(:,14,q)=x30;  
   x31=X31((q-1)*n+1:q*n,1); Xw(q,15)=x31(162); x31(r)=[]; X1(:,15,q)=x31; 
   x32=X32((q-1)*n+1:q*n,1); Xw(q,16)=x32(162); x32(r)=[]; X1(:,16,q)=x32; 
    x33=X33((q-1)*n+1:q*n,1); Xw(q,17)=x33(162); x33(r)=[]; X1(:,17,q)=x33;
    x34=X34((q-1)*n+1:q*n,1); Xw(q,18)=x34(162); x34(r)=[]; X1(:,18,q)=x34;  
    x35=X35((q-1)*n+1:q*n,1); Xw(q,19)=x35(162); x35(r)=[]; X1(:,19,q)=x35; 
    x36=X36((q-1)*n+1:q*n,1); Xw(q,20)=x36(162); x36(r)=[]; X1(:,20,q)=x36;
     x37=X37((q-1)*n+1:q*n,1); Xw(q,21)=x37(162); x37(r)=[]; X1(:,21,q)=x37;
      x38=X38((q-1)*n+1:q*n,1); Xw(q,22)=x38(162); x38(r)=[]; X1(:,22,q)=x38;
      x39=X39((q-1)*n+1:q*n,1); Xw(q,23)=x39(162); x39(r)=[]; X1(:,23,q)=x39;
       x40=X40((q-1)*n+1:q*n,1); Xw(q,24)=x40(162); x40(r)=[]; X1(:,24,q)=x40;
      x41=X41((q-1)*n+1:q*n,1); Xw(q,25)=x41(162); x41(r)=[]; X1(:,25,q)=x41;
      x42=X42((q-1)*n+1:q*n,1); Xw(q,26)=x42(162); x42(r)=[]; X1(:,26,q)=x42;
      x43=X43((q-1)*n+1:q*n,1); Xw(q,27)=x43(162); x43(r)=[]; X1(:,27,q)=x43;
       x44=X44((q-1)*n+1:q*n,1); Xw(q,28)=x44(162); x44(r)=[]; X1(:,28,q)=x44;
        x45=X45((q-1)*n+1:q*n,1); Xw(q,29)=x45(162); x45(r)=[]; X1(:,29,q)=x45;
         x46=X46((q-1)*n+1:q*n,1); Xw(q,30)=x46(162); x46(r)=[]; X1(:,30,q)=x46;
%    
     x48=X48((q-1)*n+1:q*n,1);  x48(r)=[]; Z1(:,1)=x48;  
    x49=X49((q-1)*n+1:q*n,1);  x49(r)=[]; Z1(:,2)=x49;   x50=X50((q-1)*n+1:q*n,1);  x50(r)=[]; Z1(:,3)=x50;  
    x51=X51((q-1)*n+1:q*n,1);  x51(r)=[]; Z1(:,4)=x51;   x52=X52((q-1)*n+1:q*n,1);  x52(r)=[]; Z1(:,5)=x52;
     x53=X53((q-1)*n+1:q*n,1);  x53(r)=[]; Z1(:,6)=x53;   x54=X54((q-1)*n+1:q*n,1);  x54(r)=[]; Z1(:,7)=x54; 
    x55=X55((q-1)*n+1:q*n,1);  x55(r)=[]; Z1(:,8)=x55;    x56=X56((q-1)*n+1:q*n,1);  x56(r)=[]; Z1(:,9)=x56;  
    x57=X57((q-1)*n+1:q*n,1);  x57(r)=[]; Z1(:,10)=x57;   x58=X58((q-1)*n+1:q*n,1);  x58(r)=[]; Z1(:,11)=x58;
    x59=X59((q-1)*n+1:q*n,1);  x59(r)=[]; Z1(:,12)=x59;   x60=X60((q-1)*n+1:q*n,1);  x60(r)=[]; Z1(:,13)=x60;  
    x61=X61((q-1)*n+1:q*n,1);  x61(r)=[]; Z1(:,14)=x61;   x62=X62((q-1)*n+1:q*n,1);  x62(r)=[]; Z1(:,15)=x62;   
    x63=X63((q-1)*n+1:q*n,1);  x63(r)=[]; Z1(:,16)=x63;   x64=X64((q-1)*n+1:q*n,1);  x64(r)=[]; Z1(:,17)=x64;
    x65=X65((q-1)*n+1:q*n,1);  x65(r)=[]; Z1(:,18)=x65;       x66=X66((q-1)*n+1:q*n,1);  x66(r)=[]; Z1(:,19)=x66; 
     x67=X67((q-1)*n+1:q*n,1);  x67(r)=[]; Z1(:,20)=x67;    x68=X68((q-1)*n+1:q*n,1);  x68(r)=[]; Z1(:,21)=x68;  
 %X1(:,:,q)=[X11((q-1)*m+1:q*m,1) X12((q-1)*m+1:q*m,1) X13((q-1)*m+1:q*m,1) X14((q-1)*m+1:q*m,1) X15((q-1)*m+1:q*m,1) X16((q-1)*m+1:q*m,1) X17((q-1)*m+1:q*m,1)];
end

Xp(r,:)=[];
Z1(:,22:48)=Xp;


T=99;

% load spatial weights based upon population flow across cities
load spatial_weights_9_29.mat;
load spatial_weights_lag_9_29.mat;


W5=cell(T,1); WW5=cell(T,1); % population flow from wuhan to other cities, 5 days lag
W5L=cell(T,1);  


B1=zeros(n1-1,T);  Xwfl=zeros(T,1);  B2=zeros(n1-1,T);

   for q=1:T
       
       bb1=bo1((q-1)*n+1:q*n,1); bb1(r)=[]; B1(:,q)=bb1;
       bb2=bo2((q-1)*n+1:q*n,1); bb2(r)=[]; B2(:,q)=bb2;
     x99=X99((q-1)*n+1:q*n,1); Xw(q,31)=x99(162); x99(r)=[]; index(:,q)=x99;
     x100=X100((q-1)*n+1:q*n,1); Xw(q,31)=x100(162); x100(r)=[]; indexc(:,q)=x100;
     
   end
   

   
   
  for t=1:T
    %row-normalized population flow from Wuhan
  ww1=WWm1{t};
  
  wwr=ww1(:,162);
  wwr(r)=[];
  WW5{t}=wwr;
  
    ww1(r,:)=[]; ww1(:,r)=[];
  W5{t}=ww1;
  
    
 end

for t=1:T
    %lagged row-normalization matrices
   w1=WWm1L{t};

    w1(r,:)=[]; w1(:,r)=[];
     
   W5L{t}=w1;
end
   
   
   
    




taul0=38;
taur0=20; 
taud0=[49; 74];
taug0=35;
taui0=16;


lambda10=0.26; lambda20=0.05; 
betai10=0.18; betai20=0;
rho10=0.6;  rho20=0.2; 
delta10=0; delta20=1; delta30=0.4;
mu0=0;
gamma10=0.08; gamma20=0.02; 
 a20=0.4; 
beta10=0.5*ones(30,1);  beta20=0.5*ones(48,1); a30=0.5*ones(31,1);
sigmav0=1; sigmaw0=1; a10=1;
lambda0=[lambda10;lambda20]; rho0=[rho10;rho20];
delta0=[delta10;delta20;delta30];  gamma0=[gamma10;gamma20];
betai0=[betai10;betai20];
Time10=zeros(T,1);

for q=3:T
Time10(q,1)=0.05;    
end





nit=20000;
br=0.2;

nomit=nit*br;
sce=2;
%1.25% domestic transmission after their own change-points
%2.50% domestic transmission after their own change-points
%3.75% domestic transmission after their own change-points
%4.100% domestic transmission after their own change-points




results=mcmcbfnsbr1c5(nit,br, n1,T,Y,YL,Yw,YwL, X1,B1,B2,index,indexc,Z1,Xw,W5,W5L,WW5, lambda0, rho0,mu0, delta0,betai0,beta10,beta20,gamma0, sigmav0, Time10,sigmaw0,taul0,taur0,taud0,taug0,taui0,sce);
 lambda1s= results.lambda1;
  lambda2s= results.lambda2;
 rho1s= results.rho1;
rho2s=results.rho2;
mus=results.mu;
 delta1s= results.delta1;
 delta2s=results.delta2;
  delta3s=results.delta3;
  betai1s=results.betai1;
  betai2s=results.betai2;
  gamma1s=results.gamma1;
gamma2s=results.gamma2;
sigmavs=results.sigmav;
 sigmaws=results.sigmaw;
  tauls=results.taul;
 taurs=results.taur;
 taud1s=results.taud1;
 taud2s=results.taud2;
 taugs=results.taug; 
 tauis=results.taui;  
 Ycs=results.Yc;

 hperc = 0.95;
 



% plot newly confirmed cases
tn = datetime(2020,1,21) + caldays(0:T-1);
dn=datenum(tn);

%**********************************************************
 % generate counterfactual for cumulative cases in figure 3
 %*********************************************************
 
 % compute the actual national cumulative case 
Ycumn=zeros(T,1); 
for t=1:T
  Ycumn(t,1)=sum(Ycum(:,t));  
end

%derive the mcmc draws of the predicted counterfactual national cumulative case
Yct=zeros(T,nit-nomit);
 Yct(1,:)=Ycs(1,:)+ones(1,nit-nomit)*sum(YL(:,1));
 for t=2:T
   Yct(t,:)=Yct(t-1,:)+Ycs(t,:);
    
 end

 %compute the mean counterfactual national cumulative case, and CI
Yccumn=zeros(T,1);

Ycm=mean(Ycs,2);  

Yccumn(1,1)=Ycm(1,1)+sum(YL(:,1));
for t=2:T
   Yccumn(t,1)=Yccumn(t-1,1)+Ycm(t,1);
    
end

 
 
 
 
 YctU=zeros(T,1);  YctL=zeros(T,1);
 
 
 for t=1:T
     bound=hpdi(Yct(t,:)',hperc);
     YctU(t,1)=bound(1);  YctL(t,1)=bound(2);
 end   
 



   plot(dn,Ycumn(1:T),'k');
 datetick('x','yyyy-mm-dd');
 hold on;
 plot(dn,Yccumn(1:T),'--k');
 datetick('x','yyyy-mm-dd');  
 hold on;
 plot(dn,YctU,'o');
 hold on;
 plot(dn,YctL,'*');
 
 
