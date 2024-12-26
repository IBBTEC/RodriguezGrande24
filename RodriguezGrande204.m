% G¡llespie Model for conjugation under Holling´s Type II dynamics
% Dec 2024
%REACTIONS
    %R---->R+1 T=alpha0*R Receptor division
    %D---->D+1 T=alpha*D Donor division
    %T---->T+1 T=alpha2*T Transconjugant division
%
%   D,R---->C+1,D-1,R-1  T=gamma*D*R Formation of conjugative pair
%   T,R---->CC+1,T-1,R-1 T=gamma*T*R Formation of conjugative pair from T
%   C--->T+1,D+1,C-1     T= omega*C Resolution of conjugative pair
%   CC-->D+1,D+1,CC-1    T=omega*CC Resolution of conjugative pair from D

% In this model there is continuous growth fromD and T R is assumed
% Concentrations are expressed in terms of a volume of 10 microliters
% constant
%-----------------------------------------------------------------------

clear all
 %Initial conditions 
 t=0;
 R=4E4;
 D=3E3;
 C=0;
 CC=0;
 T=0;
 
 % Propensities
 alpha0=0; %Growth Rate for R
 alpha=1; % Growth rate for  D
 alpha2=0.95; % Growth rate for T
 kon=0.02; % Pair formation rate
 omega=0.05; % Pair resolution rate
 dd=10e-6;% dilution level 

% Iterations, time and statistic vectors
muestreo=1000;% This parameter is used to record the species number in the x matrix every muestreo iterations. 

N=1000000; % Number of iterations 

x=zeros(floor(N/muestreo),6); % Space allocation for the results matrix x
tau=0;

ii=1;
step=muestreo/(R+D+T)/(alpha); 
x(ii,1:6)=[t R D T C CC];
rates=[alpha0*R alpha*D alpha2*T kon*dd*D*R kon*dd*T*R omega*C omega*CC];
y(ii,1:7)=(rates);
for i=1:N % N iterations
   
    % rates vector
rates= [alpha0*R alpha*D alpha2*T kon*dd*D*R kon*dd*T*R omega*C omega*CC];

kT=sum(rates); 

% Record time statistics


u_time=rand; % Pick random numbers
u_reaction=rand;

tau=-1/kT*log(1-u_time); % Update time statistics
t=t+tau;

if t>=ii*step
    ii=ii+1;
    x(ii,1:6)=[t R D T C CC];
    y(ii,1:7)=(rates);
end

av=0; % Guillespie pointers
j=1;

while av<kT*u_reaction
    av=av+rates(j);
    j=j+1;
end

   switch j-1
       case 1
           % R birth
          R=R+1;
          
 disp('R birth')

 case 2
           % D birth
          D=D+1;       
       disp('D birth')   
 case 3
            % T birth
           T=T+1;
           
            disp('T birth')
 case 4
           % C formation
           C=C+1;
          D=D-1;
          R=R-1;
         disp('C formation')
  case 5 
              % CC formation
            CC=CC+1;
            T=T-1;
            R=R-1;
          disp('CC formation')
  case 6
                % T from D
            C=C-1;
            T=T+1;
            D=D+1;
          disp('D Conjugation')
  case 7
                % T from T
            CC=CC-1;
            T=T+2;
          disp('T Conjugation')
            otherwise
           disp('error')
   
      end


end

