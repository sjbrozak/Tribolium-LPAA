%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bifurcation and lyapunov exponent diagrams for parameter mu_A.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc

back_step = 100;
Nmax = 50000;

nparam = 500;
parameter_vec = linspace(0,1,nparam);

%% LPA MODEL
L_plot1 = zeros(nparam,back_step+1);
P_plot1 = zeros(nparam,back_step+1);
A_plot1 = zeros(nparam,back_step+1);

% using medians of experimental fittings
              % b           c1              c2          c3                  mu_l    mu_a            
param_list = [20            0.0179          0.0003      1.076E-13           0.6053  0.0842];

disp("Generating bifurcation diagrams...")
for i=1:length(parameter_vec)

    param_temp = parameter_vec(i);

    param_list(end) = param_temp;

    % generate bifurcation plot
    [L,P,A,~] = LPA(param_list,Nmax);

    L_plot1(i,:) = L(end-back_step:end);
    P_plot1(i,:) = P(end-back_step:end);
    A_plot1(i,:) = A(end-back_step:end);
    
end

% Lyapunov exponent
parameter_vec = linspace(0,1,nparam);

lyap_LPA = zeros(length(parameter_vec),1);

disp("Generating Lyapunov diagrams...")
parfor(j = 1:length(parameter_vec),8)
    disp("j = "+j)
    lyap_params = param_list;
    lyap_params(end) = parameter_vec(j);

    % iterate map
    [L,P,A,J,lyap] = LPA(lyap_params,Nmax);
    lyap_LPA(j) = lyap;

%     lyap_LPA(j) = (1/Nmax)*log(norm(J,'fro'));
    
end

save LPA_mua.mat lyap_LPA L_plot1 P_plot1 A_plot1

%% LPAA MODEL
L_plot2 = zeros(nparam,back_step+1);
P_plot2 = zeros(nparam,back_step+1);
A1_plot2 = zeros(nparam,back_step+1);
A2_plot2 = zeros(nparam,back_step+1);

%             b             c1              c2          mu_l    mu_p      mu_a
param_list = [6.4232        0.0099          0.0028      0.6053  2.635E-12 0.0358];


disp("Generating bifurcation diagrams...")
for i=1:length(parameter_vec)

    param_temp = parameter_vec(i);

    param_list(end) = param_temp;

    % generate bifurcation plot
    [L,P,A1,A2,~,~] = LPAA(param_list,Nmax);

    L_plot2(i,:) = L(end-back_step:end);
    P_plot2(i,:) = P(end-back_step:end);
    A1_plot2(i,:) = A1(end-back_step:end);
    A2_plot2(i,:) = A2(end-back_step:end);
    
end

% Lyapunov exponent
lyap_LPAA = zeros(length(parameter_vec),1);

disp("Generating Lyapunov diagrams...")
parfor(j = 1:length(parameter_vec),8)
    disp("j = "+j)
    lyap_params = param_list;
    lyap_params(end) = parameter_vec(j);

    % iterate map
    [L,P,A1,A2,J,lyap] = LPAA(lyap_params,Nmax);
    lyap_LPAA(j) = lyap;

%     lyap_LPAA(j) = (1/Nmax)*log(norm(J,'fro'));
    
end

save LPAA_mua.mat lyap_LPAA L_plot2 P_plot2 A1_plot2 A2_plot2





%% PLOT
set(groot, 'defaultAxesTickLabelInterpreter','none');
set(groot, 'defaultLegendInterpreter','none');
set(0,'defaulttextInterpreter','none','defaultAxesFontSize',18) 

figure()
    t = tiledlayout(2,2);

    nexttile
    plot(parameter_vec,L_plot1+P_plot1+A_plot1,'k.','MarkerSize',5);
    ylabel("Total final population size","FontSize",24,'FontName','Times New Roman')
    set(gca,'FontName','Times New Roman')
    title("(a) LPA")
    ax = gca;ax.FontSize = 24;


    nexttile
    plot(parameter_vec,L_plot2+P_plot2+A1_plot2+A2_plot2,'k.','MarkerSize',5);
    ylabel("Total final population size","FontSize",24,'FontName','Times New Roman')
    set(gca,'FontName','Times New Roman')
    title("(b) LPAA")
    ax = gca;ax.FontSize = 24;


    nexttile
    plot(parameter_vec,lyap_LPA,'k-','LineWidth',2);
    yline(0,'r:','LineWidth',2);
    ylabel("Lyapunov exponent",'FontName','Times New Roman')
    ylim([-0.04 0.04])
    set(gca,'FontName','Times New Roman')
    ax = gca;ax.FontSize = 24;


    nexttile
    plot(parameter_vec,lyap_LPAA,'k-','LineWidth',2);
    yline(0,'r:','LineWidth',2)
    ylabel("Lyapunov exponent",'FontName','Times New Roman')
    ylim([-0.04 0.04])
    set(gca,'FontName','Times New Roman')
    ax = gca;ax.FontSize = 24;

    xlabel(t,"$\mu_a$ (Natural adult mortality)",'interpreter','latex',"FontSize",24)






%% functions
function [L,P,A1,A2,J,lyap] = LPAA(params,Nmax)
    
    b = params(1);
    c1 = params(2);
    c2 = params(3);
    mu_l = params(4);
    mu_p = params(5);
    mu_a = params(6);


    L = zeros((Nmax),1);
    P = zeros((Nmax),1);
    A1 = zeros((Nmax),1);
    A2 = zeros((Nmax),1);

    L(1) = 0;
    P(1) = 0;
    A1(1) = 20;
    A2(1) = 0;


    J1 = evalJ_LPAA(params,L(1),P(1),A1(1),A2(1));
    s1 = norm(J1,'fro');
    S1 = J1/s1;
    Stemp = S1;
    Stemp_prod = Stemp;
    norm(Stemp_prod);

    lyap = 0;

    for n = 1:(Nmax-1)       % Main computational loop
        L(n+1) = A2(n)*b*exp(-c1*A2(n));
        P(n+1) = (1 - mu_l)*L(n);
        A1(n+1) = (1 - mu_p)*P(n);
        A2(n+1) = A1(n)*exp(-c2*A2(n)) + (1 - mu_a)*A2(n);

        Jtemp = evalJ_LPAA(params,L(n+1),P(n+1),A1(n+1),A2(n+1));
        stemp = norm(Jtemp*Stemp_prod,'fro');
        Stemp = Jtemp/stemp;
        Stemp_prod = Stemp*Stemp_prod;
        
        
        lyap = lyap + log(stemp);
    end

    lyap = 1/Nmax*lyap;
    J = Jtemp;

end

function [L,P,A,J,lyap] = LPA(params,Nmax)
    % standard LPA model

    b = params(1);
    c1 = params(2);
    c2 = params(3);
    c3 = params(4);
    mu_l = params(5);
    mu_a = params(6);

    L = zeros(Nmax,1);
    P = zeros(Nmax,1);
    A = zeros(Nmax,1);

    L(1) = 0;
    P(1) = 0;
    A(1) = 20;

    J1 = evalJ_LPA(params,L(1),P(1),A(1));
    J = J1;

    s1 = norm(J1,'fro');
    S1 = J1/s1;
    Stemp = S1;
    Stemp_prod = Stemp;
    norm(Stemp_prod);

    lyap = log(s1);
    
    for n = 1:Nmax-1       % Main computational loop
        L(n+1) = A(n)*b*exp(-(c1*A(n)) - c2*L(n));
        P(n+1) = (1 - mu_l)*L(n);
        A(n+1) = P(n)*exp(-c3*A(n)) + (1 - mu_a)*A(n);

        Jtemp = evalJ_LPA(params,L(n+1),P(n+1),A(n+1));
    %     J = Jtemp*J;
        stemp = norm(Jtemp*Stemp_prod,'fro');
        Stemp = Jtemp/stemp;
        Stemp_prod = Stemp*Stemp_prod;
        norm(Stemp_prod);
        
        lyap = lyap + log(stemp);
    end

    lyap = 1/Nmax*lyap;
    J = Jtemp;

end




function matrix = evalJ_LPA(params,L,P,A)
    b = params(1);
    c1 = params(2);
    c2 = params(3);
    c3 = params(4);
    mu_l = params(5);
    mu_a = params(6);   

    g1 = exp(-c2*L-c1*A);
    g2 = exp(-c3*A);

    matrix = [-c2*b*A*g1 0 b*(1-c1*A)*g1;
              (1-mu_l) 0 0;
              0 g2 -c3*P*g2+(1-mu_a)];


end

function matrix = evalJ_LPAA(params,L,P,A1,A2)
    b = params(1);
    c1 = params(2);
    c2 = params(3);
    mu_l = params(4);
    mu_p = params(5);
    mu_a = params(6);   

    g1 = exp(-c1*A2);
    g2 = exp(-c2*A2);


    matrix = [ 0 0 0 b*(1-c1*A2)*g1;
              (1-mu_l) 0 0 0;
               0 (1-mu_p) 0 0;
               0 0 g2 -c2*A1*g2+(1-mu_a)];
end