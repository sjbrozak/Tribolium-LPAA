%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fits experimental subsets of the data to the LPAA model through
% minimizing the residuals of one-step forecasts.
%
% To use:
%   (1)   Specify experimental subgroup on line 18.
%   (2)   Hit "Run".
%   (2.5) This program uses parallel processing. To turn this off, set
%   parallelflag = 0 on line 62.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc
rng('default');
set(0, 'DefaultAxesFontName', 'Times New Roman');
%% 

% specify experimental group and parse out info for plotting
group = "P_ub_05";
exp_info = split(group,"_");

% read data
warning('off','all')
Data = readtable("tribolium_data_pilot.xlsx","Sheet",group);
warning('on','all')

%% Data processing

TimeStep = 1:max(Data.Time);
Week = (TimeStep-1)*2;  % Convert to weeks, starting at week 0.



% -------------- MEANS ----------------
Means = table; % Preallocation of table of means for each age class

Means.Larvae = zeros(max(Data.Time),1);
Means.Pupae = zeros(max(Data.Time),1);
Means.Adult = zeros(max(Data.Time),1);

for i = TimeStep % Calculate the means for each age class
    Means.Larvae(i) = mean(Data.Larvae(Data.Time == i));
    Means.Pupae(i) = mean(Data.Pupae(Data.Time == i), "omitnan"); % Omit 1 missing data point at step 5
    Means.Adult(i) = mean(Data.Adult(Data.Time == i));
end

% ------ compute mortality for larvae from data ------
mu_l = mean(1 - (Means.Pupae(3:end)./Means.Larvae(2:end-1)));

%% OPTIMIZATION

%     b    c1   c2   mu_p/c3  mu_a    
lb = [1    0    0    0        0   ];
ub = [20   1    1    1        1   ];

% approx. initial parameter guesses based on Constantino et al. (1998)
p0 = [10 0.01 0.01 0.2 0.2];

% weights for data: [larvae pupae adult]
weights = [1/3 1/3 1/3];

% use parallel computing
parallelflag = 1;

% ---------------- LPA MODEL -----------------
% model flag LPA = 0
problemLPA = createOptimProblem('fmincon','x0',p0,...
    'objective',@(params)obj_fun(params,Means.Larvae,Means.Pupae,Means.Adult,mu_l,0,weights), ...
                            'lb',lb,'ub',ub);

msLPA = MultiStart("Display","iter","UseParallel",parallelflag);
[best_params_LPA,fvalLPA] = run(msLPA,problemLPA,5000);

% ---------------- LPAA MODEL -----------------
% model flag LPAA = 1
problemLPAA = createOptimProblem('fmincon','x0',p0,...
    'objective',@(params)obj_fun(params,Means.Larvae,Means.Pupae,Means.Adult,mu_l,1,weights), ...
                            'lb',lb,'ub',ub);

msLPAA = MultiStart("Display","iter","UseParallel",parallelflag);
[best_params_LPAA,fvalLPAA] = run(msLPAA,problemLPAA,5000);


%% SIMULATE AND PLOT

LPA_onestep_pairs_L = zeros(max(Data.Time)-1,2);
LPA_onestep_pairs_P = zeros(max(Data.Time)-1,2);
LPA_onestep_pairs_A = zeros(max(Data.Time)-1,2);

LPAA_onestep_pairs_L = zeros(max(Data.Time)-1,2);
LPAA_onestep_pairs_P = zeros(max(Data.Time)-1,2);
LPAA_onestep_pairs_A = zeros(max(Data.Time)-1,2);

for k = 2:length(Means.Larvae)-1
    % LPA Model
    IC_LPA_temp = [Means.Larvae(k) Means.Pupae(k) Means.Adult(k)];
    [La,Pa,Aa] = LPA(best_params_LPA,2,mu_l,IC_LPA_temp);
    LPA_onestep_pairs_L(k,:) = La; LPA_onestep_pairs_P(k,:) = Pa; LPA_onestep_pairs_A(k,:) = Aa;

    % LPAA Model
    IC_LPAA_temp = [Means.Larvae(k) Means.Pupae(k) (Means.Adult(k)-Means.Adult(k-1)) Means.Adult(k-1)];
    [Lb,Pb,A1b,A2b] = LPAA(best_params_LPAA,2,mu_l,IC_LPAA_temp);
    LPAA_onestep_pairs_L(k,:) = Lb; LPAA_onestep_pairs_P(k,:) = Pb; LPAA_onestep_pairs_A(k,:) = A1b+A2b;
end

best_params_LPA = best_params_LPA'
disp("------")
best_params_LPAA = best_params_LPAA'

% plot
compareplots(Week,Means.Larvae,Means.Pupae,Means.Adult, ...
             LPA_onestep_pairs_L,LPA_onestep_pairs_P,LPA_onestep_pairs_A, ...
             LPAA_onestep_pairs_L,LPAA_onestep_pairs_P,LPAA_onestep_pairs_A, ...
             exp_info,fvalLPA,fvalLPAA)

%% FUNCTIONS
function err = obj_fun(params,Ldata,Pdata,Adata,mu_l,flag,weights)
    
    if flag == 0 % simulate LPA model
        LPA_onestep_pairs_L = [];
        LPA_onestep_pairs_P = [];
        LPA_onestep_pairs_A = [];
        for k = 2:length(Ldata)-1
            IC_LPA_temp = [Ldata(k) Pdata(k) Adata(k)];
            [La,Pa,Aa] = LPA(params,2,mu_l,IC_LPA_temp);
            LPA_onestep_pairs_L = [LPA_onestep_pairs_L; La']; 
            LPA_onestep_pairs_P = [LPA_onestep_pairs_P; Pa'];  
            LPA_onestep_pairs_A = [LPA_onestep_pairs_A; Aa']; 
        end
        Q(1) = sum(( Ldata(3:end) - LPA_onestep_pairs_L(:,end) ).^2);
        Q(2) = sum(( Pdata(3:end) - LPA_onestep_pairs_P(:,end) ).^2);
        Q(3) = sum(( Adata(3:end) - LPA_onestep_pairs_A(:,end) ).^2);

    else % simulate LPAA model
        LPAA_onestep_pairs_L = [];
        LPAA_onestep_pairs_P = [];
        LPAA_onestep_pairs_A = [];
        
        for k = 2:length(Ldata)-1
            IC_LPAA_temp = [Ldata(k) Pdata(k) (Adata(k)-Adata(k-1)) Adata(k-1)];
            [Lb,Pb,A1b,A2b] = LPAA(params,2,mu_l,IC_LPAA_temp);
            LPAA_onestep_pairs_L = [LPAA_onestep_pairs_L; Lb']; 
            LPAA_onestep_pairs_P = [LPAA_onestep_pairs_P; Pb'];  
            LPAA_onestep_pairs_A = [LPAA_onestep_pairs_A; A1b'+A2b'];  
        end
    
        Q(1) = sum(( Ldata(3:end) - LPAA_onestep_pairs_L(:,end) ).^2);
        Q(2) = sum(( Pdata(3:end) - LPAA_onestep_pairs_P(:,end) ).^2);
        Q(3) = sum(( Adata(3:end) - LPAA_onestep_pairs_A(:,end) ).^2);
    end

    err = weights(1)*Q(1) + weights(2)*Q(2) + weights(3)*Q(3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These iteration functions for the LPA and LPAA model aren't strictly
% necessary because we aren't fitting the entire dataset, but we're keeping
% it for generality.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L,P,A1,A2] = LPAA(params,Nmax,mu_l,IC)

    b = params(1);
    c1 = params(2);
    c2 = params(3);
    mu_p = params(4);
    mu_a = params(5);

    L = zeros(Nmax,1);
    P = zeros(Nmax,1);
    A1 = zeros(Nmax,1);
    A2 = zeros(Nmax,1);

    L(1) = IC(1);
    P(1) = IC(2);
    A1(1) = IC(3);
    A2(1) = IC(4);

    for n = 1:Nmax-1       % Main computational loop
        L(n+1) = A2(n)*b*exp(-c1*A2(n));
        P(n+1) = (1 - mu_l)*L(n);
        A1(n+1) = (1 - mu_p)*P(n);
        A2(n+1) = A1(n)*exp(-c2*A2(n)) + (1 - mu_a)*A2(n);
    end
end


function [L,P,A] = LPA(params,Nmax,mu_l,IC)

    b = params(1);
    c1 = params(2);
    c2 = params(3);
    c3 = params(4);
    mu_a = params(5);

    L = zeros(Nmax,1);
    P = zeros(Nmax,1);
    A = zeros(Nmax,1);

    L(1) = IC(1);
    P(1) = IC(2);
    A(1) = IC(3);

    for n = 1:Nmax-1       % Main computational loop
        L(n+1) = A(n)*b*exp(-(c1*A(n)) - c2*L(n));
        P(n+1) = (1 - mu_l)*L(n);
        A(n+1) = P(n)*exp(-c3*A(n)) + (1 - mu_a)*A(n);
    end
end

function compareplots(Week,Ldata,Pdata,Adata,La,Pa,Aa,Lb,Pb,Ab,titlestr,fval1,fval2)
    stepWeek = Week;
    
    NP = titlestr(1);
    if titlestr(2) == "b"
        flour = "bleached";
    else 
        flour = "unbleached";
    end

    if titlestr(3) == "05"
        perc = "0.5%";
    else
        perc = "1%";
    end

    figure()
        tl = tiledlayout(3,2);
        title(tl,NP+" "+perc+" "+flour,'fontsize',30,'FontName','Times New Roman');
        xlabel(tl,'Week','fontsize',24,'FontName','Times New Roman')

        % ----------- LARVAE ---------------
        nexttile; hold on; 
        plot(Week(1:length(Ldata)), Ldata, 'o-', 'MarkerSize', 8, 'MarkerEdgeColor', 'b',...
        'MarkerFaceColor', [0.7 0.7 1], 'LineWidth', 2.5,'Color', 'b');
        for j = 2:length(La)
            plot([stepWeek(j),stepWeek(j+1)], La(j,:),'b:o', 'LineWidth', 2.5)
        end
        ylabel("Larvae");
        xlim([Week(1) Week(end)]);
        ax = gca;ax.FontSize = 24;
        title(sprintf("(a) LPA (Error: %.2e)",fval1))
        ylim([0 500]);

        nexttile; hold on; 
        plot(Week(1:length(Ldata)), Ldata, 'o-', 'MarkerSize', 8, 'MarkerEdgeColor', 'b',...
        'MarkerFaceColor', [0.7 0.7 1], 'LineWidth', 2.5,'Color', 'b');
        for j = 2:length(Lb)
            plot([stepWeek(j),stepWeek(j+1)], Lb(j,:),'b:o', 'LineWidth', 2.5)
        end
        xlim([Week(1) Week(end)]);
        ax = gca;ax.FontSize = 24;
        title(sprintf("(a) LPAA (Error: %.2e)",fval2))
        ylim([0 500]);

        % ----------- PUPAE ---------
        nexttile; hold on; 
        plot(Week(1:length(Ldata)), Pdata, 'o-', 'MarkerSize', 8, 'MarkerEdgeColor', 'r',...
        'MarkerFaceColor', [1 0.7 0.7], 'LineWidth', 2.5,'Color', 'r');
        for j = 2:length(Pa)
            plot([stepWeek(j),stepWeek(j+1)], Pa(j,:),'r:o', 'LineWidth', 2.5)
        end
        ylabel("Pupae")
        xlim([Week(1) Week(end)]);
        ax = gca;ax.FontSize = 24;
        ylim([0 300]);

        nexttile; hold on; 
        plot(Week(1:length(Ldata)), Pdata, 'o-', 'MarkerSize', 8, 'MarkerEdgeColor', 'r',...
        'MarkerFaceColor', [1 0.7 0.7], 'LineWidth', 2.5,'Color', 'r');
        for j = 2:length(Pb)
            plot([stepWeek(j),stepWeek(j+1)], Pb(j,:),'r:o', 'LineWidth', 2.5)
        end
        xlim([Week(1) Week(end)]);
        ax = gca;ax.FontSize = 24;
        ylim([0 300]);


        % ----------- ADULTS ---------
        nexttile; hold on; 
        plot(Week(1:length(Ldata)), Adata, 'o-', 'MarkerSize', 8, 'MarkerEdgeColor', 'k',...
        'MarkerFaceColor', [0.5 0.5 0.5], 'LineWidth', 2.5,'Color', 'k');
        for j = 2:length(Aa)
            plot([stepWeek(j),stepWeek(j+1)], Aa(j,:),'k:o', 'LineWidth', 2.5)
        end
        ylabel("Adults")
        xlim([Week(1) Week(end)]);
        legend("Data","Model",'location','northwest')
        ax = gca;ax.FontSize = 24;
        ylim([0 600]);

        nexttile; hold on; 
        plot(Week(1:length(Ldata)), Adata, 'o-', 'MarkerSize', 8, 'MarkerEdgeColor', 'k',...
        'MarkerFaceColor', [0.5 0.5 0.5], 'LineWidth', 2.5,'Color', 'k');
        for j = 2:length(Ab)
            plot([stepWeek(j),stepWeek(j+1)], Ab(j,:),'k:o', 'LineWidth', 2.5)
        end
        xlim([Week(1) Week(end)]);
        legend("Data","Model",'location','northwest')
        ax = gca;ax.FontSize = 24;
        ylim([0 600]);


        figure()
            tq = tiledlayout(1,3);
            ax1 = nexttile;
            h1 = qqplot(Ldata(2:end) - Lb(:,2));
            h1(1).LineWidth = 2.5; h1(2).LineWidth = 2.5; h1(3).LineWidth = 2.5;
            h1(1).MarkerSize = 8;
            ax = gca;ax.FontSize = 16;
            ylabel(ax1,'Larvae residuals')
            title(ax1,'(a)')
        
        
            ax2 = nexttile;
            h2 = qqplot(Pdata(2:end) - Pb(:,2));
            h2(1).LineWidth = 2.5; h2(2).LineWidth = 2.5; h2(3).LineWidth = 2.5;
            h2(1).MarkerSize = 8;
            ax = gca;ax.FontSize = 16;
            ylabel(ax2,'Pupae residuals')
            title(ax2,'(b)')
        
            ax3 = nexttile;
            h3 = qqplot(Adata(2:end) - Ab(:,2));
            h3(1).LineWidth = 2.5; h3(2).LineWidth = 2.5; h3(3).LineWidth = 2.5;
            h3(1).MarkerSize = 8;
            ax = gca;ax.FontSize = 16;
            ylabel(ax3,'Adult residuals')
            title(ax3,'(c)')
        
            title(tq,"QQ plot of LPAA residuals ("+NP+" "+perc+" "+flour+")",'FontName','Times New Roman','fontSize',20);
end


