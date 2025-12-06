% LangEvin Equation (LEE)
% LEE:A Physics-Inspired Optimizer based on LangEvin Equation
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
% LEE:A Physics-Inspired Optimizer based on LangEvin Equation
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

clear 
close all
clc
nP=20;          % Number of Population

Func_name='F1'; % Name of the test function

MaxIt=750;      % Maximum number of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Func_name);

[Best_fitness,BestPositions,Convergence_curve] = LEE(nP,MaxIt,lb,ub,dim,fobj);

%% Draw objective space

figure,
hold on
semilogy(Convergence_curve,'Color','r','LineWidth',4);
title('Convergence curve')
xlabel('Iteration');
ylabel('Best fitness obtained so far');
axis tight
grid off
box on
legend('LEE')


