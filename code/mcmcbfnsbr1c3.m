function [results] =mcmcbfnsbr1c3(nit,br, n,T,Y,YL,Yw,YwL, X1,B1,B2,index,indexc,Z1,Xw,W1,W1L,WWr, lambda0, rho0,mu0, delta0,betai0,beta10,beta20,gamma0, sigmav0, Time10,sigmaw0,taul0,taur0,taud0,taug0,taui0,sce)
% =========================================================================
% This function conducts the MCMC estimation of a higher-order SAR model with spatial errors.
%==========================================================================

nomit=nit*br; %burn-in sample
  
%store all parameters
lambda1s=zeros(nit,1); lambda2s=zeros(nit,1); 
rho1s=zeros(nit,1);  rho2s=zeros(nit,1); 
mus=zeros(nit,1);   
delta1s=zeros(nit,1);   delta2s=zeros(nit,1); delta3s=zeros(nit,1);
betai1s=zeros(nit,1);  betai2s=zeros(nit,1);
beta1s=zeros(30,nit);  beta2s=zeros(48,nit); 
gamma1s=zeros(nit,1); gamma2s=zeros(nit,1); 
sigmavs=zeros(nit,1);  sigmaws=zeros(nit,1); a2s=zeros(nit,1); a3s=zeros(31,nit); 
  Time1s=zeros(T,nit);  %store factor for all periods 
  a1s=zeros(nit,1);
tauls=zeros(nit,1);
taurs=zeros(nit,1); 
taud1s=zeros(nit,1); taud2s=zeros(nit,1);
taugs=zeros(nit,1);  tauis=zeros(nit,1);



%store transform data
Yr=zeros(n-1,T);     Yc=zeros(T,nit);   
 Res=zeros(n-1,T); 

%=======================================
%prepare for the AM algorithm for lambda 
%=======================================
KKt=5;   st1=0.1/sqrt(KKt); ratio=0.05;
sum1=zeros(5,1);  sum2=zeros(5,5);

%==========================
% value of prior parameters
%===========================
TT=200; %initial period of AM algorithm 
P=10;  % for beta1 and beta2   
a=0.001;b=0.001;  %for sigma

theta0=[lambda0;rho0;mu0];



for i=1: nit
    
     if mod(i,100)==0
        fprintf('%d\n',i)
     end
    
     
     %**************************************
     % Gibbs step for a2 and a3 in the AR(1) of Wuhan
     %***************************************
      %sample beta0=[delta0; beta20]
 sumb=zeros(33,1); sumbv=zeros(33,33);
 
for t=1:T
 ZZ=[1,YwL(t), Xw(t,:)];
 sumb=sumb+ZZ'*(Yw(t))/sigmaw0;
 sumbv=sumbv+ZZ'*ZZ/sigmaw0;

end

sumv=(sumbv+eye(33)/(P))\eye(33);
%sumv=mean(cat(3,sumv,sumv'),3); 

Tr=sumv*sumb;

at=mvnrnd(Tr,sumv); 
a0=at;
a0=a0';
a10=a0(1);
a20=a0(2); 
a30=a0(3:33);

a1s(i,1)=a10;
a2s(i,1)=a20;  a3s(:,i)=a30;


% sampling sigmaw in AR(1) for wuhan
 sums=0;
 for t=1:T
   ZZ=[1,YwL(t), Xw(t,:)];
  ss=Yw(t)-ZZ*a0;
  sums=sums+ss'*ss;   
 end
 
 
 ap=a+T/2;
 bp=b+0.5*sums;
sigmaw0= 1./gamrnd(ap,1/bp);  
sigmaws(i,1)=sigmaw0;   
     
     
     

   
   %AM algorithm for lambda, rho and mu
   
   if (i<=2*TT)
      accept=0;
      theta1=mvnrnd(theta0,st1^2*eye(5));
      theta1=theta1';
   lambda1=theta1(1:2);
      rho1=theta1(3:4);
      mu1=theta1(5);
     
      
      lambdam=max(abs(lambda1));
      rhom=max(abs(rho1));
      mum=max(abs(mu1));
     
     
      
      while (accept==0) %reject bounds on lambda1 and rho1
          if (lambdam<1 ) && (lambdam+rhom+mum<1)
              accept=1;
          else
            theta1=mvnrnd(theta0,st1^2*eye(5));
      theta1=theta1';
   lambda1=theta1(1:2);
      rho1=theta1(3:4);
      mu1=theta1(5);
     
     
      
      lambdam=max(abs(lambda1));
      rhom=max(abs(rho1));
      mum=max(abs(mu1));
     
          end
      end
  end
  
  if (i>2*TT)
      accept=0;
      if (i<=nomit)
      vvarr=(1-ratio)^2*2.38^2*varr/KKt+ratio^2*st1^2*eye(5);
      end
      theta1=mvnrnd(theta0,vvarr);
      theta1=theta1';
 
    lambda1=theta1(1:2);
      rho1=theta1(3:4);
      mu1=theta1(5);
     
      
      lambdam=max(abs(lambda1));
      rhom=max(abs(rho1));
      mum=max(abs(mu1));
     
      
      while (accept==0) %reject bounds on lambda1
          if (lambdam<1 ) && (lambdam+rhom+mum<1)
              accept=1;
          else
      theta1=mvnrnd(theta0,vvarr);
    theta1=theta1';
 
  lambda1=theta1(1:2);
      rho1=theta1(3:4);
      mu1=theta1(5);
     
      
      lambdam=max(abs(lambda1));
      rhom=max(abs(rho1));
      mum=max(abs(mu1));
     
          end
      end
  end
  
 frpp=1;
 

  
  for t=1: T
      
    lambdan0=lambda0(1)*(t<=taul0)+lambda0(2)*(1-(t<=taul0));   
   rhon0=rho0(1)*(t<=taur0)+rho0(2)*(1-(t<=taur0));  
    mun0=mu0;
     deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1))-(t>taud0(2)))+delta0(3)*(t>taud0(2));  
     betain0=betai0(1)*(t<=taui0)+betai0(2)*(1-(t<=taui0));
     gamman0=gamma0(1)*(t<=taug0)+gamma0(2)*(1-(t<=taug0)); 
    
 lambdan1=lambda1(1)*(t<=taul0)+lambda1(2)*(1-(t<=taul0));   
   rhon1=rho1(1)*(t<=taur0)+rho1(2)*(1-(t<=taur0));  
  mun1=mu1;
    
    
    %deltan1=delta1(1)*(t<=tau20)+delta1(2)*(1-(t<=tau20));  
     %gamman1=gamma1(1)*(t<=tau10)+gamma1(2)*(1-(t<=tau10)-(t>=tau20))+gamma1(3)*(t>=tau20) ;
     
      
      
      
      
      
  frpp=likn(Y(:,t), YL(:,t), Yw(t), X1(:,:,t), B1(:,t),index(:,t), Z1,  W1{t},W1L{t},  WWr{t}, Time10(t),lambdan0,rhon0,mun0,lambdan1,rhon1,mun1,deltan0,deltan0,betain0,betain0,beta10,beta20,gamman0,gamman0,sigmav0)*frpp;
  end
    

    
   % Determine the transition probability
    Acceptr=min(1,frpp);
   
    % Draw from uniform(0,1)
    u1=rand(1,1); 
    
   % Transition to candidate delta11 with probability acceptr
     if (Acceptr>u1) 
        lambda1_1=lambda1(1);
        lambda2_1=lambda1(2);
        rho1_1=rho1(1);
        rho2_1=rho1(2);
    
        mu_1=mu1;
       
       
     else
       lambda1_1=lambda0(1);
        lambda2_1=lambda0(2);
        rho1_1=rho0(1);
        rho2_1=rho0(2);
     
       mu_1=mu0;
       
       
         
     end
    
   % Store the value and continue
       rho10=rho1_1;
     rho1s(i,1)=rho10;
     
       rho20=rho2_1;
     rho2s(i,1)=rho20;

       
     mu0=mu_1;
     mus(i,1)=mu0;
     
  
     lambda10=lambda1_1;
     lambda1s(i,1)=lambda10;
     
      lambda20=lambda2_1;
      lambda2s(i,1)=lambda20;  
      
      
      lambda0=[lambda10;lambda20];
      rho0=[rho10;rho20];
      
      
      
%calculate the empirical covariance
    theta0=[lambda0;rho0; mu0];
    
    if (i<=nomit)
        sum1=sum1+theta0;
        sum2=sum2+theta0*theta0';
    end
    
    if (i>1)&&(i<=nomit)
        mean1=sum1/i;
        varr=sum2/i-mean1*mean1';
    end
   
for t=1:T
    lambdan0=lambda0(1)*(t<=taul0)+lambda0(2)*(1-(t<=taul0));   
    rhon0=rho0(1)*(t<=taur0)+rho0(2)*(1-(t<=taur0));  
    
    SS=eye(n-1)-lambdan0*W1{t};

   Yr(:,t)=SS*Y(:,t)-rhon0*YL(:,t)-mu0*W1L{t}*YL(:,t);   
end
    

   


  %sample beta0=[delta0; beta10; beta20]
 sumb=zeros(37,1); sumbv=zeros(37,37);
 
for t=1:T
 ZZ=[B1(:,t).*(t<=taud0(1)), B1(:,t).*(1-(t<=taud0(1))-(t>taud0(2))), B1(:,t).*(t>taud0(2)), index(:,t)*(t<=taui0), index(:,t)*(t>taui0), WWr{t}*Yw(t).*(t<=taug0),  WWr{t}*Yw(t).*(1-(t<=taug0)),  X1(:,:,t)];
 sumb=sumb+ZZ'*(Yr(:,t)-Z1*beta20-ones(n-1,1)*Time10(t))/sigmav0;
 sumbv=sumbv+ZZ'*ZZ/sigmav0;

end

sumv=(sumbv+eye(37)/(P))\eye(37);
%sumv=mean(cat(3,sumv,sumv'),3); 

Tr=sumv*sumb;

betat=mvnrnd(Tr,sumv); 


beta0=betat';


delta0=beta0(1:3);
betai0=beta0(4:5);
gamma0=beta0(6:7);
beta10=beta0(8:37);

delta1s(i,1)=delta0(1);
delta2s(i,1)=delta0(2);
delta3s(i,1)=delta0(3);
betai1s(i,1)=betai0(1);
betai2s(i,1)=betai0(2);
gamma1s(i,1)=gamma0(1);
gamma2s(i,1)=gamma0(2);

beta1s(:,i)=beta10;




%Sample sigmav
 sums=0;
 for t=1:T
 ZZ=[B1(:,t).*(t<=taud0(1)), B1(:,t).*(1-(t<=taud0(1))-(t>taud0(2))), B1(:,t).*(t>taud0(2)), index(:,t)*(t<=taui0), index(:,t)*(t>taui0), WWr{t}*Yw(t).*(t<=taug0),  WWr{t}*Yw(t).*(1-(t<=taug0)),  X1(:,:,t)];
  ss=Yr(:,t)-ZZ*beta0-Z1*beta20-ones(n-1,1)*Time10(t);
  sums=sums+ss'*ss;   
 end
 
 
 ap=a+(n-1)*T/2;
 bp=b+0.5*sums;
sigmav0= 1./gamrnd(ap,1/bp);  
sigmavs(i,1)=sigmav0;   


 %Gibbs sampling for time-fixed effects
    
    for t=2:T
   ZZ=[B1(:,t).*(t<=taud0(1)), B1(:,t).*(1-(t<=taud0(1))-(t>taud0(2))), B1(:,t).*(t>taud0(2)), index(:,t)*(t<=taui0), index(:,t)*(t>taui0), WWr{t}*Yw(t).*(t<=taug0),  WWr{t}*Yw(t).*(1-(t<=taug0)),  X1(:,:,t)];
        Yf=Yr(:,t)-ZZ*beta0-Z1*beta20;
        
        Sigmaf=(1+(n-1)/(sigmav0))^(-1);
       
        Tf=Sigmaf*ones(n-1,1)'*Yf/(sigmav0);
        
        time0=mvnrnd(Tf,Sigmaf);
        
        Time10(t,1)=time0;
       
    end
    
      Time1s(:,i)=Time10;
    



    
 %sample beta20
 sumr=zeros(48,1); sumrv=zeros(48,48);
 for t=1:T
   
     sumrv=sumrv+Z1'*Z1/sigmav0;
 ZZ=[B1(:,t).*(t<=taud0(1)), B1(:,t).*(1-(t<=taud0(1))-(t>taud0(2))), B1(:,t).*(t>taud0(2)), index(:,t)*(t<=taui0), index(:,t)*(t>taui0), WWr{t}*Yw(t).*(t<=taug0),  WWr{t}*Yw(t).*(1-(t<=taug0)),  X1(:,:,t)];
Tc=Yr(:,t)-ZZ*beta0-ones(n-1,1)*Time10(t);
     sumr=sumr+Z1'*Tc/sigmav0;
 end
 
 sumv=(sumrv+eye(48)/(P))\eye(48);
%sumv=mean(cat(3,sumv,sumv'),3); 

Tr=sumv*sumr;

beta2t=mvnrnd(Tr,sumv); 
beta20=beta2t';
beta2s(:,i)=beta20;

% M-H step for taur0


   taur1=14+unidrnd(84); 

          
    frpt=1;
  
  for t=1: T
      
   lambdan0=lambda0(1)*(t<=taul0)+lambda0(2)*(1-(t<=taul0));   
   rhon0=rho0(1)*(t<=taur0)+rho0(2)*(1-(t<=taur0));  
    mun0=mu0;
     deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1))-(t>taud0(2)))+delta0(3)*(t>taud0(2));  
     gamman0=gamma0(1)*(t<=taug0)+gamma0(2)*(1-(t<=taug0)); 
     betain0=betai0(1)*(t<=taui0)+betai0(2)*(1-(t<=taui0));
     
     rhon1=rho0(1)*(t<=taur1)+rho0(2)*(1-(t<=taur1));  
   
     
  frpt=likn(Y(:,t), YL(:,t), Yw(t), X1(:,:,t), B1(:,t),index(:,t), Z1,  W1{t},W1L{t},  WWr{t}, Time10(t),lambdan0,rhon0,mun0,lambdan0,rhon1,mun0,deltan0,deltan0,betain0,betain0,beta10,beta20,gamman0,gamman0,sigmav0)*frpt;
  end
    
   
           
    % Determine the transition probability
    Acceptr=min(1,frpt);
   
    % Draw from uniform(0,1)
    u1=rand(1,1); 
    
   % Transition to candidate delta11 with probability acceptr
     if (Acceptr>u1)
       
         taur_1=taur1;
        
   
         
     else
     
         taur_1=taur0;
       
    
     end
     
  
  

      taur0=taur_1;
   
     taurs(i,1)=taur0;
     
   
 % M-H step for taul
  
 taul1=14+unidrnd(84);

          
    frpt=1;
  
  for t=1: T
      
   lambdan0=lambda0(1)*(t<=taul0)+lambda0(2)*(1-(t<=taul0));   
   rhon0=rho0(1)*(t<=taur0)+rho0(2)*(1-(t<=taur0));  
    mun0=mu0;
     deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1))-(t>taud0(2)))+delta0(3)*(t>taud0(2));  
     gamman0=gamma0(1)*(t<=taug0)+gamma0(2)*(1-(t<=taug0)); 
    betain0=betai0(1)*(t<=taui0)+betai0(2)*(1-(t<=taui0));
      lambdan1=lambda0(1)*(t<=taul1)+lambda0(2)*(1-(t<=taul1));   

    frpt=likn(Y(:,t), YL(:,t), Yw(t), X1(:,:,t), B1(:,t),index(:,t), Z1,  W1{t},W1L{t},  WWr{t}, Time10(t),lambdan0,rhon0,mun0,lambdan1,rhon0,mun0,deltan0,deltan0,betain0,betain0,beta10,beta20,gamman0,gamman0,sigmav0)*frpt;
  end
    
   
           
    % Determine the transition probability
    Acceptr=min(1,frpt);
   
    % Draw from uniform(0,1)
    u1=rand(1,1); 
    
   % Transition to candidate delta11 with probability acceptr
     if (Acceptr>u1)
         taul_1=taul1;
        
     
         
     else
       taul_1=taul0;
     
    
     end
     
     taul0=taul_1;
     tauls(i,1)=taul0;
  

    % M-H step for taug
   taug1=14+unidrnd(84);

 frpt=1;
  
  for t=1: T
      
   lambdan0=lambda0(1)*(t<=taul0)+lambda0(2)*(1-(t<=taul0));   
   rhon0=rho0(1)*(t<=taur0)+rho0(2)*(1-(t<=taur0));  
    mun0=mu0;
     deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1))-(t>taud0(2)))+delta0(3)*(t>taud0(2));  
     gamman0=gamma0(1)*(t<=taug0)+gamma0(2)*(1-(t<=taug0)); 
    betain0=betai0(1)*(t<=taui0)+betai0(2)*(1-(t<=taui0));
     
     gamman1=gamma0(1)*(t<=taug1)+gamma0(2)*(1-(t<=taug1)); 
     
     
  frpt=likn(Y(:,t), YL(:,t), Yw(t), X1(:,:,t), B1(:,t), index(:,t),Z1,  W1{t},W1L{t},  WWr{t}, Time10(t),lambdan0,rhon0,mun0,lambdan0,rhon0,mun0,deltan0,deltan0,betain0,betain0,beta10,beta20,gamman0,gamman1,sigmav0)*frpt;
  end
    

           
    % Determine the transition probability
    Acceptr=min(1,frpt);
   
    % Draw from uniform(0,1)
    u1=rand(1,1); 
    
   % Transition to candidate delta11 with probability acceptr
     if (Acceptr>u1)
       
         taug_1=taug1;
    
         
     else
      
         taug_1=taug0;
    
     end
     
   
     
   
      taug0=taug_1;
     taugs(i,1)=taug0;
     
     % M-H step for taui
   
  taui1=14+unidrnd(84);

 frpt=1;
  
  for t=1: T
      
   lambdan0=lambda0(1)*(t<=taul0)+lambda0(2)*(1-(t<=taul0));   
   rhon0=rho0(1)*(t<=taur0)+rho0(2)*(1-(t<=taur0));  
    mun0=mu0;
     deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1))-(t>taud0(2)))+delta0(3)*(t>taud0(2));  
     gamman0=gamma0(1)*(t<=taug0)+gamma0(2)*(1-(t<=taug0)); 
    betain0=betai0(1)*(t<=taui0)+betai0(2)*(1-(t<=taui0));
     
     betain1=betai0(1)*(t<=taui1)+betai0(2)*(1-(t<=taui1)); 
      
  frpt=likn(Y(:,t), YL(:,t), Yw(t), X1(:,:,t), B1(:,t), index(:,t),Z1,  W1{t},W1L{t},  WWr{t}, Time10(t),lambdan0,rhon0,mun0,lambdan0,rhon0,mun0,deltan0,deltan0,betain0,betain1,beta10,beta20,gamman0,gamman0,sigmav0)*frpt;
  end
    

           
    % Determine the transition probability
    Acceptr=min(1,frpt);
   
    % Draw from uniform(0,1)
    u1=rand(1,1); 
    
   % Transition to candidate delta11 with probability acceptr
     if (Acceptr>u1)
       
         taui_1=taui1;
     
         
     else
      
         taui_1=taui0;
    
     end
     
   
     
   
      taui0=taui_1;
     tauis(i,1)=taui0;   
     
     
     
     
     
     % M-H for taud
  

  taud11=40+unidrnd(30);  
  taud12=71+unidrnd(27);
  
  frpt=1;
  
  for t=1: T
      
   lambdan0=lambda0(1)*(t<=taul0)+lambda0(2)*(1-(t<=taul0));   
   rhon0=rho0(1)*(t<=taur0)+rho0(2)*(1-(t<=taur0));  
    mun0=mu0;
     deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1))-(t>taud0(2)))+delta0(3)*(t>taud0(2));  
     gamman0=gamma0(1)*(t<=taug0)+gamma0(2)*(1-(t<=taug0)); 
     betain0=betai0(1)*(t<=taui0)+betai0(2)*(1-(t<=taui0)); 
     
     deltan1=delta0(1)*(t<=taud11)+delta0(2)*(1-(t<=taud11)-(t>taud12))+delta0(3)*(t>taud12);  
     
  frpt=likn(Y(:,t), YL(:,t), Yw(t), X1(:,:,t), B1(:,t),index(:,t), Z1,  W1{t},W1L{t},  WWr{t}, Time10(t),lambdan0,rhon0,mun0,lambdan0,rhon0,mun0,deltan0,deltan1,betain0,betain0,beta10,beta20,gamman0,gamman0,sigmav0)*frpt;
  end
    

           
    % Determine the transition probability
    Acceptr=min(1,frpt);
   
    % Draw from uniform(0,1)
    u1=rand(1,1); 
    
   % Transition to candidate delta11 with probability acceptr
     if (Acceptr>u1)
         
         taud_1=[taud11;taud12];
        
     
     else
    
         taud_1=taud0;
   
     end
     
    
     
     taud0=taud_1;
     taud1s(i,1)=taud0(1);
     taud2s(i,1)=taud0(2);
     
     if (i>nomit)
         %*******************************
         % obtain the residual estimates
         %*******************************
         for t=1:T
             lambdan0=lambda0(1)*(t<=taul0)+lambda0(2)*(t>taul0);
             rhon0=rho0(1)*(t<=taur0)+rho0(2)*(1-(t<=taur0));
             gamman0=gamma0(1)*(t<=taug0)+gamma0(2)*(t>taug0);
             deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1))-(t>taud0(2)))+delta0(3)*(t>taud0(2));
             betain0=betai0(1)*(t<=taui0)+betai0(2)*(1-(t<=taui0));
             SS=eye(n-1)-lambdan0*W1{t};
             Yr=SS*Y(:,t)-rhon0*YL(:,t)-mu0*W1L{t}*YL(:,t);
             ZZ=[B1(:,t), index(:,t), WWr{t}*Yw(t), X1(:,:,t),Z1];
             
             Res(:,t)=Yr-ZZ*[deltan0; betain0; gamman0;beta10; beta20]-Time10(t,1)*ones(n-1,1);
            
         end
         
         YLc=zeros(n-1,T);   YLc(:,1)=YL(:,1);
         
         
         if sce==1
             
             for t=1:T
                 lambdan0=lambda0(1)*(t<=taul0)+lambda0(2)*(1-(t<=taul0));
                 rhon0=rho0(1)*(t<=taur0)+rho0(2)*(t>taur0);
                 gamman0=gamma0(1)*(t<=taug0)+gamma0(2)*(t>taug0);
                 deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1)));
                 betain0=betai0(1)*(t<=taui0)+betai0(2)*(1-(t<=taui0));
                 SS=eye(n-1)-lambdan0*W1{t};
                 SSi=SS\eye(n-1);
                 ZZ=[B1(:,t),index(:,t), WWr{t}*Yw(t), X1(:,:,t),Z1];
                 xs=(rhon0*eye(n-1)+mu0*W1L{t})*YLc(:,t)+ZZ*[deltan0;betain0;gamman0;beta10; beta20]+Time10(t,1)*ones(n-1,1)+Res(:,t);
                 ylc=SSi*xs;
                 Yc(t,i)=sum(ylc);
                 if (t<T)
                     YLc(:,t+1)=ylc;
                 end
                 
             end
             
         elseif sce==2
             
             for t=1:T
                 lambdan0=lambda0(1)*(t<=taul0)+lambda0(2)*(1-(t<=taul0));
                 rhon0=rho0(1)*(t<=taur0)+rho0(2)*(t>taur0);
                 gamman0=gamma0(1)*(t<=taug0)+gamma0(2)*(t>taug0);
                 deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1))-(t>taud0(2)))+delta0(3)*(t>taud0(2));
                 betain0=betai0(1)*(t<=taui0)+betai0(2)*(1-(t<=taui0));
                 SS=eye(n-1)-lambdan0*W1{t};
                 SSi=SS\eye(n-1);
                 ZZ=[B2(:,t), index(:,t), WWr{t}*Yw(t), X1(:,:,t),Z1];
                 
                 xs=(rhon0*eye(n-1)+mu0*W1L{t})*YLc(:,t)+ZZ*[deltan0;betain0;gamman0;beta10; beta20]+Time10(t,1)*ones(n-1,1)+Res(:,t);
                 
                 ylc=SSi*xs;
                 Yc(t,i)=sum(ylc);
                 if (t<T)
                     YLc(:,t+1)=ylc;
                 end
             end
         elseif sce==3
             for t=1:T
                 lambdan0=lambda0(1)*(t<=taul0)+lambda0(2)*(1-(t<=taul0));
                 rhon0=rho0(1)*(t<=taur0)+rho0(2)*(t>taur0);
                 gamman0=gamma0(1)*(t<=taug0)+gamma0(2)*(t>taug0);
                 deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1)));
                 betain0=betai0(1)*(t<=taui0)+betai0(2)*(1-(t<=taui0));
                 SS=eye(n-1)-lambdan0*W1{t};
                 SSi=SS\eye(n-1);
                 ZZ=[B2(:,t),index(:,t), WWr{t}*Yw(t), X1(:,:,t),Z1];
                 xs=(rhon0*eye(n-1)+mu0*W1L{t})*YLc(:,t)+ZZ*[deltan0;betain0;gamman0;beta10; beta20]+Time10(t,1)*ones(n-1,1)+Res(:,t);
                 ylc=SSi*xs;
                 Yc(t,i)=sum(ylc);
                 if (t<T)
                     YLc(:,t+1)=ylc;
                 end
             end
         elseif sce==4
             for t=1:T
                 lambdan0=lambda0(1);
                 rhon0=rho0(1);
                 gamman0=gamma0(1);
                 deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1))-(t>taud0(2)))+delta0(3)*(t>taud0(2));
                 betain0=betai0(1);
                 SS=eye(n-1)-lambdan0*W1{t};
                 SSi=SS\eye(n-1);
                 ZZ=[B1(:,t), index(:,t), WWr{t}*Yw(t), X1(:,:,t),Z1];
                 
                 xs=(rhon0*eye(n-1)+mu0*W1L{t})*YLc(:,t)+ZZ*[deltan0;betain0;gamman0;beta10; beta20]+Time10(t,1)*ones(n-1,1)+Res(:,t);
                 
                 ylc=SSi*xs;
                 Yc(t,i)=sum(ylc);
                 if (t<T)
                     YLc(:,t+1)=ylc;
                 end
             end
             
             
         elseif sce==5
             for t=1:T
                 lambdan0=lambda0(1)*(t<=taul0)+lambda0(2)*(1-(t<=taul0));
                 rhon0=rho0(1)*(t<=taur0)+rho0(2)*(t>taur0);
                 gamman0=gamma0(1)*(t<=taug0)+gamma0(2)*(t>taug0);
                 deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1))-(t>taud0(2)))+delta0(3)*(t>taud0(2));
                 betain0=betai0(1)*(t<=taui0)+betai0(2)*(1-(t<=taui0));
                 SS=eye(n-1)-lambdan0*W1{t};
                 SSi=SS\eye(n-1);
                 ZZ=[B1(:,t),indexc(:,t), WWr{t}*Yw(t), X1(:,:,t),Z1];
                 xs=(rhon0*eye(n-1)+mu0*W1L{t})*YLc(:,t)+ZZ*[deltan0;betain0;gamman0;beta10; beta20]+Time10(t,1)*ones(n-1,1)+Res(:,t);
                 ylc=SSi*xs;
                 Yc(t,i)=sum(ylc);
                 if (t<T)
                     YLc(:,t+1)=ylc;
                 end
                 
             end
             
         else
             for t=1:T
                 lambdan0=lambda0(1);
                 rhon0=rho0(1);
                 gamman0=gamma0(1);
                 deltan0=delta0(1)*(t<=taud0(1))+delta0(2)*(1-(t<=taud0(1))-(t>taud0(2)))+delta0(3)*(t>taud0(2));
                 betain0=betai0(1);
                 SS=eye(n-1)-lambdan0*W1{t};
                 SSi=SS\eye(n-1);
                 ZZ=[B1(:,t), indexc(:,t), WWr{t}*Yw(t), X1(:,:,t),Z1];
                 
                 xs=(rhon0*eye(n-1)+mu0*W1L{t})*YLc(:,t)+ZZ*[deltan0;betain0;gamman0;beta10; beta20]+Time10(t,1)*ones(n-1,1)+Res(:,t);
                 
                 ylc=SSi*xs;
                 Yc(t,i)=sum(ylc);
                 if (t<T)
                     YLc(:,t+1)=ylc;
                 end
             end
             
             
         end
         
         
         
     end

end
    

  
results.Time1=Time1s(:,nomit+1:nit);

 results.lambda1=lambda1s(nomit+1:nit,1);
  results.lambda2=lambda2s(nomit+1:nit,1);
 results.rho1=rho1s(nomit+1:nit,1);
  results.rho2=rho2s(nomit+1:nit,1);
 
  results.mu=mus(nomit+1:nit,1);
   results.delta1=delta1s(nomit+1:nit,1);
   results.delta2=delta2s(nomit+1:nit,1);
    results.delta3=delta3s(nomit+1:nit,1);
    results.betai1=betai1s(nomit+1:nit,1);
    results.betai2=betai2s(nomit+1:nit,1);
  results.beta1=beta1s(:,nomit+1:nit);
 results.gamma1=gamma1s(nomit+1:nit,1);
 results.gamma2=gamma2s(nomit+1:nit,1);
 results.beta2=beta2s(:,nomit+1:nit);
 results.sigmav=sigmavs(nomit+1:nit,1);
  results.sigmaw=sigmaws(nomit+1:nit,1);
  results.a2=a2s(nomit+1:nit,1);
  results.a1=a1s(nomit+1:nit,1);
  results.a3=a3s(:,nomit+1:nit);

  results.taul=tauls(nomit+1:nit,1);
  results.taui=tauis(nomit+1:nit,1);
results.taur=taurs(nomit+1:nit,1);
 results.taug=taugs(nomit+1:nit,1);
results.taud1=taud1s(nomit+1:nit,1);
results.taud2=taud2s(nomit+1:nit,1);
results.Yc=Yc(:,nomit+1:nit,1);
end

% =========================================================================
% support functions below
% =========================================================================
 

function [fvalue ] = likn(Yq, Ylq, Ywq, X1q, B1q,index,Z1,  Wrq,WrLq, WWrq, Time1q,  lambdan0,rhon0,mun0,lambdan1,rhon1,mun1,delta0,delta1,betai0, betai1,beta10,beta20,gamma0,gamma1,sigmav0)
n=length(Yq);

S0=eye(n)-lambdan0*Wrq;
S1=eye(n)-lambdan1*Wrq;

C0=S0*Yq-rhon0*Ylq-mun0*WrLq*Ylq-B1q*delta0-index*betai0-X1q*beta10-Z1*beta20-WWrq*Ywq*gamma0-ones(n,1)*Time1q;
C1=S1*Yq-rhon1*Ylq-mun1*WrLq*Ylq-B1q*delta1-index*betai1-X1q*beta10-Z1*beta20-WWrq*Ywq*gamma1-ones(n,1)*Time1q;

CC0=C0'*C0/(2*sigmav0);
CC1=C1'*C1/(2*sigmav0);

fvalue=(det(S1)/det(S0))*exp(-CC1+CC0);


end
