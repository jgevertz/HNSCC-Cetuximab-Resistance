%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Santiago Cardenas 04/26/2021                                          %
% Numerically fits Model 1.1 with no pre-existing or acquired           %
% resistance to each treatment mouse (using censored data). Algorithm   %
% is run 15 times, and saves the best parameter set per mouse from the  %
% 15 runs. Parallel toolbox required, with access to 29 parallel pools. %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
path = 'Best_Fit_Model1_1';
if exist(path, 'dir') ~= 7
    mkdir(path)
end

%% Variables
filename = '../Censored_Data/CTXData_censored.xlsx';
A = xlsread(filename);
A(A == 0) = NaN; %first row is not numbers
A(1,1) = 0; %fills in first cell as Day 0
num_mice = size(A,2)-1; 
numreps = 15;
Volume = {};
Days = {};

for n = 2:num_mice+1
    Volume_orig = A(:,n);
    Days_orig = A(:,1);
    for vol = 1:size(Volume_orig,1)
        if (isnan(Volume_orig(vol,1))==1)
            Days_orig(vol,1) = NaN;
        end
    end
    Volume{n-1}(:) = rmmissing(Volume_orig); %removes any entry that contains missing data
    Days{n-1}(:) = rmmissing(Days_orig);
end

% Best parameters for each mouse
Best_r_s = zeros(num_mice,numreps + 1); Best_r_r = zeros(num_mice,numreps + 1);
Best_g = zeros(num_mice,numreps + 1); Best_lambda_s = zeros(num_mice,numreps + 1); gamma = log(0.5)/-4.75;
Best_lambda_r = zeros(num_mice,numreps + 1); Best_SSE = zeros(num_mice,numreps + 1); Best_y0 = zeros(num_mice,numreps + 1);
Best_y0r = zeros(num_mice,numreps + 1);
Figures = gobjects(num_mice,numreps);

%% Sobol parameters
min_r_s = 0; max_r_s = .2;
min_r_r = 0; max_r_r = .1;
min_g = 0; max_g = 0;
min_lambda_s = 0; max_lambda_s = .1;
min_lambda_r = 0; max_lambda_r = .1; 
n_Sbpts = 1500000; % 1.5e6
n_SApts = 500000; % 5e5
SA_SSEs = {};

%% Initialize parallel pool
numcores = num_mice;
parpool('local',numcores); % Call to open the distributed processing % 29 on cluster
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0
else
    poolsize = poolobj.NumWorkers
end

parfor m = 1:num_mice
    Volume_loop = Volume{m};
    Days_loop = Days{m};
    n = size(Volume{m},2);
    
    for rep = 1:numreps
        %% Set parameters for mouse m
        dosage =  ceil((Days_loop(end)/7)); %calculates dosage for each mouse based on days
        min_y0 = 0;
        max_y0 = 2*Volume_loop(1);
        min_y0r = 0; % min_y0;
        max_y0r = 0; % max_y0;
        
        a = clock;
        rng(1000*a(5)*a(6));
        n_skip = randi([1000000, 10000000]); % 1e6 to 1e7
        n_leap = 0;
        
        %% Create sobol points and scale them
        % clearvars uniform_sobol_scaled
        uniform_sobol = sobolset(7,'Skip',n_skip,'Leap',n_leap)
        uniform_sobol = net(uniform_sobol,n_Sbpts);
        uniform_sobol_scaled = zeros(n_Sbpts,7);
        
        uniform_sobol_scaled(:,1) = (max_r_s-min_r_s)*uniform_sobol(:,1) + min_r_s; %r_s
        uniform_sobol_scaled(:,2) = (max_r_r-min_r_r)*uniform_sobol(:,2) + min_r_r; %r_r
        uniform_sobol_scaled(:,3) = (max_g-min_g)*uniform_sobol(:,3) + min_g; %g
        uniform_sobol_scaled(:,4) = (max_lambda_s-min_lambda_s)*uniform_sobol(:,4) + min_lambda_s; %lambda_s
        uniform_sobol_scaled(:,5) = (max_lambda_r-min_lambda_r)*uniform_sobol(:,5) + min_lambda_r; %lambda_r
        uniform_sobol_scaled(:,6) = (max_y0-min_y0)*uniform_sobol(:,6) + min_y0; % S0
        uniform_sobol_scaled(:,7) = (max_y0r-min_y0r)*uniform_sobol(:,7) + min_y0r; % R0
        SSE_QMC = [];
        
        for u = 1:n_Sbpts % n_pts  % Solve DE at all biologically-relevant Sobol points
            if ((uniform_sobol_scaled(u,4) > uniform_sobol_scaled(u,5))&&(uniform_sobol_scaled(u,1)>=uniform_sobol_scaled(u,2)))
                sb_soln = {}; % sobol solution vector
                d1 = 1; %initialize d
                
                %create first soln vector
                sb_soln{d1} = ode23s(@(t,x) ExpDrugModel(t,x,uniform_sobol_scaled(u,1),...
                    uniform_sobol_scaled(u,2),uniform_sobol_scaled(u,3), uniform_sobol_scaled(u,4),...
                    uniform_sobol_scaled(u,5),gamma),0:7,...
                    [uniform_sobol_scaled(u,6),uniform_sobol_scaled(u,7),5]);
                t1 = 7;
                t2 = 14;
                
                %loop through d times to create d sol vecs
                for d = 2:dosage
                    sb_soln{d} = ode23s(@(t,x) ExpDrugModel(t,x,uniform_sobol_scaled(u,1),...
                        uniform_sobol_scaled(u,2),uniform_sobol_scaled(u,3), uniform_sobol_scaled(u,4),...
                        uniform_sobol_scaled(u,5),gamma),t1:t2,...
                        [sb_soln{d-1}.y(1,end),sb_soln{d-1}.y(2,end),sb_soln{d-1}.y(3,end)+5]);
                    t1 = t1 + 7;
                    t2 = t2 + 7;
                end
                
                %create initial time and solution vecs to concatenate
                d1 = 1;
                sb_time = sb_soln{d1}.x; % sobol time vector
                sb_solution = sb_soln{d1}.y(1,:)+sb_soln{d1}.y(2,:); % sobol total solution vector
                for d = 2:dosage
                    temp = sb_soln{d}.y(1,2:end)+sb_soln{d}.y(2,2:end); %concatenate the 2 y soln vecs together
                    sb_time = [sb_time sb_soln{d}.x(2:end);]; %concatenate time
                    sb_solution = [sb_solution temp]; % solution vec being concatenated
                end
                
                newz = 0;
                for day = 1:length(Days_loop)
                    [c , index] = min(abs(sb_time-Days_loop(day)));
                    error = (sb_solution(index) - Volume_loop(day))^2;
                    newz = newz + error;
                end
                SSE_QMC(u) = newz;
                
            else % Not biologically allowable parameter set
                SSE_QMC(u) = 10^10;
            end
            
        end
        
        % Identify the Sobol point with the lowest SSE
        [c, index] = min(SSE_QMC);
        
        %% Simulated annealing starts at best-fit Sobol point
        % these variables only exist in parfor loop and are recreated each loop
        % each is a vector that logs SA changes, values are put into SAsteps
        % variables at the end
        newParam_r_s = []; newParam_r_r = []; newParam_g = [];
        newParam_lambda_s = []; newParam_lambda_r = [];
        y0 = []; y0r = []; SSE = [];
        
        i = 1;
        newParam_r_s(i) = uniform_sobol_scaled(index,1);
        newParam_r_r(i) = uniform_sobol_scaled(index,2);
        newParam_g(i) = uniform_sobol_scaled(index,3);
        newParam_lambda_s(i) = uniform_sobol_scaled(index,4);
        newParam_lambda_r(i) = uniform_sobol_scaled(index,5);
        y0(i)=  uniform_sobol_scaled(index,6);
        y0r(i) = uniform_sobol_scaled(index,7);
        SSE(i) = c;
        

        %% Simulated annealing
        order_r_s = floor(log10(newParam_r_s(1)))-1;
        order_r_r = floor(log10(newParam_r_r(1)))-1;
        order_lambda_s = floor(log10(newParam_lambda_s(1)))-1;
        order_lambda_r = floor(log10(newParam_lambda_r(1)))-1;
        order_y0 = floor(log10(y0(1)))-1;
        order_y0r = floor(log10(y0r(1)))-1;
        
        for j = 1:n_SApts
            %% Step 1: perturb parameters
            a = -1;
            b = 1;
            r = ((b-a)*(rand)+a)*10^(order_r_s);
            s = ((b-a)*(rand)+a)*10^(order_r_r);
            s2 = ((b-a)*(rand)+a)*10^(order_lambda_s);
            q2 = ((b-a)*(rand)+a)*10^(order_lambda_r);
            v = ((b-a)*(rand)+a)*10^(order_y0);
            
            tempnewparam_r_s = newParam_r_s(i)+r;
            tempnewparam_r_r = newParam_r_r(i) + s;
            tempnewparam_g = 0; %newParam_g{m}(i)+q;
            tempnewparam_lambda_s = newParam_lambda_s(i)+s2;
            tempnewparam_lambda_r  = newParam_lambda_r(i) + q2;
            tempnewy0 = y0(i) + v;
            tempnewy0r = 0; % y0r(i) + v2;            
            
            while tempnewy0r < 0 || tempnewy0 < 0 || tempnewparam_lambda_r <0 || tempnewparam_lambda_s < 0 ...
                    || tempnewparam_g < 0 || tempnewparam_r_r < 0 || tempnewparam_r_s < 0 ...
                    || tempnewparam_lambda_r > tempnewparam_lambda_s || tempnewparam_r_s < tempnewparam_r_r
                
                a = -1;
                b = 1;
                r = ((b-a)*(rand)+a)*10^(order_r_s);
                s = ((b-a)*(rand)+a)*10^(order_r_r);
                %q = ((b-a)*(rand)+a)*0.01;
                s2 = ((b-a)*(rand)+a)*10^(order_lambda_s);
                q2 = ((b-a)*(rand)+a)*10^(order_lambda_r);
                v = ((b-a)*(rand)+a)*10^(order_y0);
                % v2 = ((b-a)*(rand)+a)*10^(order_y0r);
                
                tempnewparam_r_s = newParam_r_s(i)+r;
                tempnewparam_r_r = newParam_r_r(i) + s;
                tempnewparam_g = 0; %newParam_g{m}(i)+q;
                tempnewparam_lambda_s = newParam_lambda_s(i)+s2;
                tempnewparam_lambda_r  = newParam_lambda_r(i) + q2;
                tempnewy0 = y0(i) + v;
                tempnewy0r = 0; % y0r(i) + v2;
                
            end
            
            %% Step 2: Solve ODEs at perturbed parameters
            d1 = 1;
            sa_soln = {};
            sa_soln{d1} = ode23s(@(t,x) ExpDrugModel(t,x,tempnewparam_r_s,tempnewparam_r_r,...
                tempnewparam_g, tempnewparam_lambda_s,tempnewparam_lambda_r,gamma),0:7,...
                [tempnewy0,tempnewy0r,5]);
            t1 = 7;  t2 = 14;
            
            for d = 2:dosage
                sa_soln{d} = ode23s(@(t,x) ExpDrugModel(t,x,tempnewparam_r_s,tempnewparam_r_r,...
                    tempnewparam_g, tempnewparam_lambda_s,tempnewparam_lambda_r,gamma),t1:t2,...
                    [sa_soln{d-1}.y(1,end),sa_soln{d-1}.y(2,end),sa_soln{d-1}.y(3,end)+5]);
                t1 = t1 + 7;
                t2 = t2 + 7;
            end
            
            d1 = 1;
            sa_time = sa_soln{d1}.x;
            sa_solution = sa_soln{d1}.y(1,:)+sa_soln{d1}.y(2,:);
            for d = 2:dosage
                temp = sa_soln{d}.y(1,2:end)+sa_soln{d}.y(2,2:end);
                sa_time = [sa_time sa_soln{d}.x(2:end);];
                sa_solution = [sa_solution temp];
            end
            
            initialSSE = 0;
            newz = 0;
            for day = 1:length(Days_loop)
                [c , index] = min(abs(sa_time-Days_loop(day)));
                error = initialSSE + (sa_solution(index) - Volume_loop(day))^2;
                newz = newz + error;
            end
            
            %% Step 3: prob of acceptance
            b = 1000; %beta value
            delta = newz - SSE(i); %delta value
            if delta < 0
                probability = 1;
            else
                probability = 0; %exp(-b*d);
            end
            
            newr = rand;
            if (newr <= probability)
                i = i + 1;
                SSE(i) = newz;
                newParam_r_s(i) = tempnewparam_r_s;
                newParam_r_r(i) = tempnewparam_r_r;
                newParam_g(i) = tempnewparam_g;
                newParam_lambda_s(i) = tempnewparam_lambda_s;
                newParam_lambda_r(i) = tempnewparam_lambda_r;
                y0(i) = tempnewy0;
                y0r(i) = tempnewy0r;
            end
            
        end
        
        %% Optimal parameter set is the last accepted parameter set from simulated annealing
        Best_r_s(m,rep) = newParam_r_s(end);
        Best_r_r(m,rep) = newParam_r_r(end);
        Best_g(m,rep) = newParam_g(end);
        Best_lambda_s(m,rep) = newParam_lambda_s(end);
        Best_lambda_r(m,rep) = newParam_lambda_r(end);
        Best_SSE(m,rep) = SSE(end);
        Best_y0(m,rep) = y0(end);
        Best_y0r(m,rep) = y0r(end);
        
        d1 = 1;
        final_soln = {}
        final_soln{d1} = ode23s(@(t,x) ExpDrugModel(t,x,newParam_r_s(end),newParam_r_r(end),...
            newParam_g(end), newParam_lambda_s(end),newParam_lambda_r(end),gamma),...
            0:7,[y0(end),y0r(end),5]);
        t1 = 7; t2 = 14;
        for d = 2:dosage
            final_soln{d} = ode23s(@(t,x) ExpDrugModel(t,x,newParam_r_s(end),newParam_r_r(end),...
                newParam_g(end), newParam_lambda_s(end),newParam_lambda_r(end),gamma),...
                t1:t2,[final_soln{d-1}.y(1,end),final_soln{d-1}.y(2,end),final_soln{d-1}.y(3,end)+5]);
            t1 = t1 + 7;
            t2 = t2 + 7;
        end
        
        d1 = 1;
        final_time = final_soln{d1}.x;
        final_solution = final_soln{d1}.y(1,:)+final_soln{d1}.y(2,:);
        for d = 2:dosage
            temp = final_soln{d}.y(1,2:end)+final_soln{d}.y(2,2:end);
            final_time = [final_time final_soln{d}.x(2:end);];
            final_solution = [final_solution temp];
        end
        
        figure('visible','off');
        plot(final_time,final_solution,'b','LineWidth',2); hold on;
        plot(Days_loop(:),Volume_loop(:),'ob','MarkerSize', 10,'MarkerFaceColor','b');
        hold off;
        xlabel('Days','FontSize',16);
        ylabel('Volume','FontSize',16);
        title(['Model 1.1: Mouse ' num2str(m) ' rep ' num2str(rep)],'FontSize',18)
        
        fprintf('For Mouse %d rep %d, Best_r_s = %f, Best_r_r = %f, Best_g = %f\nBest_lambda_s = %f , Best_gamma = %f, Best_lambda_r = %f,\nBest_SSE = %f, Best_S0 = %f, Best_R0 = %f\n',...
            m, rep, Best_r_s(m,rep), Best_r_r(m,rep), Best_g(m,rep), Best_lambda_s(m,rep),...
            gamma, Best_lambda_r(m,rep), Best_SSE(m,rep),Best_y0(m,rep),Best_y0r(m,rep));
        
        Figures(m,rep) = gcf;
        SA_SSEs{m,rep} = SSE(:);
        fprintf("Mouse %d rep %d done\n\n",m,rep)
        
    end
    fprintf("Mouse %d all done\n",m) % check
end
delete(gcp('nocreate'))

BestRep = zeros(num_mice,1);
for mouse2 = 1:num_mice
    [minSSE, minrep] = min(Best_SSE(mouse2,1:numreps));
    Best_r_s(mouse2,numreps + 1) = Best_r_s(mouse2,minrep);
    Best_r_r(mouse2,numreps + 1) = Best_r_r(mouse2,minrep);
    Best_g(mouse2,numreps + 1) = Best_g(mouse2,minrep);
    Best_lambda_s(mouse2,numreps + 1) = Best_lambda_s(mouse2,minrep);
    Best_lambda_r(mouse2,numreps + 1) = Best_lambda_r(mouse2,minrep);
    Best_SSE(mouse2,numreps + 1) = Best_SSE(mouse2,minrep);
    Best_y0(mouse2,numreps + 1) = Best_y0(mouse2,minrep);
    Best_y0r(mouse2,numreps + 1) = Best_y0r(mouse2,minrep);
    SA_SSEs{mouse2,numreps + 1} = SA_SSEs{mouse2,minrep};
    BestRep(mouse2) = minrep;
    figure('visible','off');
    plot(1:size(SA_SSEs{mouse2,numreps + 1}),SA_SSEs{mouse2,numreps + 1}(:));
    xlabel('SA step','FontSize',16);
    ylabel('SSE','FontSize',16);
    title(['Model 1.2: Mouse ' num2str(mouse2) ' rep ' num2str(minrep)],'FontSize',18)
    Figures(mouse2,numreps + 1) = gcf;
    fname_fig = [path '/Mouse' num2str(mouse2) 'Rep' num2str(minrep)]; % how to put a number in string
    saveas(Figures(mouse2,minrep),[fname_fig,'.fig']);
    saveas(Figures(mouse2,minrep),[fname_fig,'.png']);
    fname_fig2 = [path '/Mouse' num2str(mouse2) 'Rep' num2str(minrep) 'SSEplot']; % how to put a number in string
    saveas(Figures(mouse2,numreps + 1),[fname_fig2,'.fig']);
    saveas(Figures(mouse2,numreps + 1),[fname_fig2,'.png']);
end

Final_r_s = Best_r_s(:,numreps + 1)'; Final_r_r = Best_r_r(:,numreps + 1)'; Final_g = Best_g(:,numreps + 1)';
Final_lambda_s = Best_lambda_s(:,numreps + 1)'; Final_lambda_r = Best_lambda_r(:,numreps + 1)';
Final_SSE = Best_SSE(:,numreps + 1)'; Final_y0 = Best_y0(:,numreps + 1)'; Final_y0r = Best_y0r(:,numreps + 1)';

filename = [path '/BestParams.mat'];
save(filename,'Best_r_s','Best_r_r','Best_g','Best_lambda_s','Best_lambda_r',...
    'Best_SSE','Best_y0','Best_y0r','BestRep','Final_r_s','Final_r_r',...
    'Final_g','Final_lambda_s','Final_lambda_r','Final_SSE','Final_y0',...
    'Final_y0r','SA_SSEs')


