function [Best_Cost,Best_X,Convergence_curve]=LEE(nP,MaxIt,lb,ub,Dim,fobj)
%---------------------------------------------------------------------------------------------------------------------------

% LangEvin Equation (LEE)
% LEE?A Physics-Inspired Optimizer based on LangEvin Equation
% Codes of LEE:http://imanahmadianfar.com/codes/
% Website of LEE:http://www.aliasgharheidari.com/LEE.html

% Huiling Chen, Iman Ahmadianfar, Ali asghar Heidari , Marjan Kordani, Arvin Samadi koucheksarae, and Guoxi Liang 

%  Last update: 12-07-2025

%  e-Mail: im.ahmadian@gmail.com,i.ahmadianfar@bkatu.ac.ir.
%  e-Mail: as_heidari@ut.ac.ir, aliasghar68@gmail.com,
%  e-Mail: kordani.1@buckeyemail.osu.edu,
%  e-Mail: Arvinsamadi.k@gmail.com,

% After use, please refer to the main paper:
% Huiling Chen, Iman Ahmadianfar, Ali asghar Heidari , Marjan Kordani, Arvin Samadi koucheksarae, and Guoxi Liang   
% LEE?A Physics-Inspired Optimizer based on LangEvin Equation
% Neurocomputing, 2025,
%---------------------------------------------------------------------------------------------------------------------------
% You can also follow the paper for related updates in researchgate: https://www.researchgate.net/profile/Iman_Ahmadianfar
% Researchgate: https://www.researchgate.net/profile/Ali_Asghar_Heidari.

%  Website of LEE:%  http://www.aliasgharheidari.com/LEE.html

% You can also use and compare with our other new optimization methods:
                                                                       %(GBO)-2020- http://www.imanahmadianfar.com/codes.
                                                                       %(RUN)-2021- http://www.aliasgharheidari.com/LEE.html
                                                                       %(INFO)-2022-http://www.aliasgharheidari.com/INFO.html
                                                                       %(HGS)-2021- http://www.aliasgharheidari.com/HGS.html
                                                                       %(SMA)-2020- http://www.aliasgharheidari.com/SMA.html
%---------------------------------------------------------------------------------------------------------------------------
%% Define Parameters and Objective Function
upr = 0.3;           
LB=lb;                                       % Lower Bound of Variables
UB=ub;                                       % Upper Bound of Variables

inT = 1;                                    % Initial Temperature
     
%% Initialization
X = LB+(UB-LB).*rand(nP, Dim);              % Random initial positions     Eq.(3)
Z = rand(nP, Dim);                          % Random initial positions
V = zeros(nP, Dim);                         % Initial velocities
Cost = zeros(nP, 1);                        % Fitness
CostNew = zeros(nP, 1);                     % Fitness
T = inT;

% Create Initial Population
for i=1:nP
      Cost(i) = fobj(X(i,:));
end
[Best_Cost,ind] = min(Cost);
Best_X = X(ind,:);

%% Main loop
for iter = 1:MaxIt
    
    Lambda = exp(-(iter/MaxIt)^(2));        % Eq.(12)
    beta = sin(((pi/2)*iter + pi*(iter/MaxIt)))*(1-iter/MaxIt) .* exp(-2*(iter/MaxIt));  % Eq.(14-2)       

    [~,ind]=sort(Cost);
      
      a0=1:nP;
      [a1, a2] = Gna1a2(nP,a0);           % Generate Random Indices                                    
      b = randn(nP,1).*beta;              % Eq.(14-1) 
    
      pr0 = upr+0.1.*rand(nP,1);
      pr1 = pr0+0.1.*sinh(randn(nP,1));   % Eq.(18-6)  
      pr2 = pr0+0.1.*sinh(randn(nP,1));   % Eq.(21) 
      F = 0.5+0.1*sinh(randn(nP,1));      % Eq.(17)
      Xpb = (X(ind(randi(4,1,nP)),:));
      
    for i=1:nP
        
        % Update particles' velocities and positions based on Langevin equation
        Xavg = (X(a1(i),:)+X(a2(i),:))/2; % Eq.(11)
        randomForce = sqrt(2 * T * Lambda) * randn.*(X(i,:)-Xavg);         % Eq.(10)
        V(i,:) = V(i,:) - Lambda * V(i,:) + randomForce;                   % Eq.(14)
        U = Xavg + F(i).*(Xpb(i,:)-X(i,:))+ b(i).*V(i,:);                  % Eq.(16)

        Stp1 = Best_X - Xavg;                                              % Eq.(18-2)                                            
        Stp2 = X(i,:) - Xavg;                                              % Eq.(18-3) 
        
        UL = rand*(UB-LB);
        Stp = Stp1+(Stp2-Stp1)/UL;                                         % Eq.(18-1) 
        
        r=rand;
        L=rand(1,Dim)<pr1(i);Sigma=(randn*(r)+(1-rand)^2);            
        Z(i,:) = L.*U + (1-L).*(X(i,:) +  Stp.*Sigma);                     % Eq.(18)            

        if rand<pr2(i)
            rdn1 = randn; rdn2 = randn;
            if rand<(1-iter/MaxIt)
                Z(i,:) = X(i,:) + F(i).*(rdn1).*(Best_X-X(a1(i),:))*Lambda;% Eq.(19) 
            else
                Z(i,:) = Best_X + (F(i).*rdn1.*(b(i).*Xpb(i,:)-X(i,:))*Lambda+F(i).*(rdn2).*(b(i).*Best_X-X(i,:))*Lambda); % Eq.(20) 
            end
        end


    %  Boundary handling 
    Flag4ub=Z(i,:)>UB;
    Flag4lb=Z(i,:)<LB;    
    Z(i,:)=(Z(i,:).*(~(Flag4ub+Flag4lb)))+UB.*Flag4ub+LB.*Flag4lb; 

        CostNew(i) = fobj(Z(i,:));
        if CostNew(i)<Cost(i)
            X(i,:) = Z(i,:);
            Cost(i) = CostNew(i);
            if Cost(i)<Best_Cost
                Best_X = X(i,:);
                Best_Cost = Cost(i);
            end
        end
    end


    T = inT * exp(-(iter/MaxIt));

    Convergence_curve(iter) = Best_Cost;
    disp(['it: ' num2str(iter) ' , BestCost: ' num2str(Best_Cost)])

end







