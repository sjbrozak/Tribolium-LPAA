%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stability plot.
%
% The bounds are obtained by separating R0 into its components.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc

% parameters from fitting
c1 = 0.0099;
c2 = 0.0028;
mu_l = 0.6053;
mu_p = 2.635E-12;
k = (1 - mu_l) * (1 - mu_p);

mu_a = linspace(0.0000001,1,200);

% E0 GAS upper bound: R0<1
E0GASub = log(mu_a) - log(k);

% Estar GAS upper bounds: R0 < min{...}:
EstarGASub1 = 1 + log(mu_a) - log(k);
EstarGASub2 = log(exp(1)*c1/c2*(1-mu_a)./mu_a) + log(mu_a) - log(k);
EstarGAS = min([EstarGASub1; EstarGASub2]);

% Estar LAS upper bounds: R0 < min{...}
EstarLASub1 = 1 + c2/c1 + log(mu_a) - log(k);
EstarLASub2 = (1-mu_a)./(mu_a)*(c1/c2) + log(mu_a) - log(k);
EstarLAS = min([EstarLASub1; EstarLASub2]);
%% simulate

%             b         c1              c2          mu_l    mu_p      mu_a
param_list = [exp(-1)	0.0099          0.0028      0.6053  2.635E-12 0.6];
case1 = LPAA(param_list,500);

param_list = [exp(1)	0.0099          0.0028      0.6053  2.635E-12 0.6];
case2 = LPAA(param_list,500);

param_list = [exp(1.6)	0.0099          0.0028      0.6053  2.635E-12 0.6];
case3 = LPAA(param_list,500);

param_list = [exp(2.9)	0.01513003467	0	0.6053	0	0.6];
case4 = LPAA(param_list,500);

xpts = [0.6, 0.6, 0.6, 0.6];
ypts = [-1, 1, 1.6, 2.9];



%% plot

figure()
    t = tiledlayout(4,5);

    nexttile(t,[4,3])
    hold on; box on
    grid on; ylim([-4 3])

    % plot Estar LAS region
    mask2 = EstarLAS > E0GASub;
    fx2 = [mu_a(mask2), fliplr(mu_a(mask2))];
    fy2 = [E0GASub(mask2), fliplr(EstarLAS(mask2))];
    fill_color2 = [0.7 0.7 1];
    fh2 = fill(fx2,fy2,fill_color2);

    plot(mu_a, max([EstarLASub1; EstarLASub2]),'--','Color',[0, 0, 1, 0.3],'LineWidth',3)
    plot(mu_a, EstarLAS,'--','Color',[0, 0, 1, 1],'LineWidth',3)

    % plot Estar GAS region
    mask2 = EstarGAS > E0GASub;
    fx2 = [mu_a(mask2), fliplr(mu_a(mask2))];
    fy2 = [E0GASub(mask2), fliplr(EstarGAS(mask2))];
    fill_color2 = [1 0.7 0.7];
    fh2 = fill(fx2,fy2,fill_color2);
    plot(mu_a, max([EstarGASub1; EstarGASub2]),'--','Color',[1, 0, 0, 0.3],'LineWidth',3)
    plot(mu_a, EstarGAS,'--','Color',[1, 0, 0, 1],'LineWidth',3)

    % plot E0 GAS region
    fakelb = -100*ones(1,200);
    mask1 = E0GASub > fakelb;
    fx1 = [mu_a(mask1), fliplr(mu_a(mask1))];
    fy1 = [fakelb(mask1), fliplr(E0GASub(mask1))];
    fill_color1 = [.9 0.9 0.9];
    fh1 = fill(fx1,fy1,fill_color1);
    plot(mu_a,E0GASub,'--','Color',"k",'LineWidth',3)

    str = {'  (e)','  (d)','  (c)','  (b)'};
    plot(xpts, ypts,'.','Color',"black",'MarkerSize',20)
    text(xpts,ypts,str,'FontSize',20,'FontName','Times New Roman')

    ylim([-2 4]);
    str1 = {'$E^*$ \textbf{G.A.S.}'};
    text(.7,1.2,str1,'FontSize',22,'Interpreter','latex');

    str3 = {'\textbf{Hypothesized upper bound for L.A.S.}'};
    text(0.1,2.8,str3,'FontSize',22,'Interpreter','latex');

    str2 = {'$E^0$ \textbf{G.A.S.}'};
    text(.7,-01,str2,'FontSize',22,'Interpreter','latex');


    ax = gca;
    ax.FontSize = 22;
    set(gca,'FontName','Times New Roman')
    title('(a)','FontSize',22)

    xlabel('$\mu_a$','Interpreter','latex','FontSize',22)
    ylabel('$\ln{b}$','Interpreter','latex','FontSize',22)



    nexttile(t,[1,2])
    plot(case4,'k','LineWidth',2); axis tight
    ax = gca;
    ax.FontSize = 18;
    set(gca,'FontName','Times New Roman')
    title('(b)  ln(b) = 2.9, \mu_a = 0.6','FontSize',22,'Interpreter','tex')
    axis padded
    xticks(0:100:500);xticklabels([0 200 400 600 800 1000])
    xlim([0,500])
    xlabel('Week')
    

    nexttile(t,[1,2])
    plot(case3,'k','LineWidth',2); axis tight
    ax = gca;
    ax.FontSize = 18;
    set(gca,'FontName','Times New Roman')
    title('(c)  ln(b) = 1.6, \mu_a = 0.6','FontSize',22,'Interpreter','tex')
    xticks(0:100:500);xticklabels([0 200 400 600 800 1000])
    axis padded
    xlim([0,500])

    nexttile(t,[1,2])
    plot(case2,'k','LineWidth',2); 
    ax = gca;
    ax.FontSize = 18;
    set(gca,'FontName','Times New Roman')
    title('(d) ln(b) = 1, \mu_a = 0.6','FontSize',22,'Interpreter','tex')
    xticks(0:100:500);xticklabels([0 200 400 600 800 1000])
    axis padded
    xlim([0,500])


    nexttile(t,[1,2])
    plot(case1,'k','LineWidth',2); 
    ax = gca;
    ax.FontSize = 18;
    set(gca,'FontName','Times New Roman')
    title('(e) ln(b) = -1, \mu_a = 0.6','FontSize',22,'Interpreter','tex')
    xticks(0:100:500);xticklabels([0 200 400 600 800 1000])
    axis padded
    xlim([0,500])





%%
function [total] = LPAA(params,Nmax)
    
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

    for n = 1:(Nmax-1)       % Main computational loop
        L(n+1) = A2(n)*b*exp(-c1*A2(n));
        P(n+1) = (1 - mu_l)*L(n);
        A1(n+1) = (1 - mu_p)*P(n);
        A2(n+1) = A1(n)*exp(-c2*A2(n)) + (1 - mu_a)*A2(n);
    end

    total = L+P+A1+A2;
end
    