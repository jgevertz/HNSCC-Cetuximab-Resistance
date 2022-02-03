%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Santiago Cardenas 09/27/2021                                          %
% Numerically fits allee curve to each control mouse (using censored    %
% data). Algorithm is run 15 times, and saves the best parameter set    %
% per mouse from the 15 runs. Parallel toolbox required, with access to %
% 25 parallel pools.                                                    %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% Clock function and saves best fit each time
cl = clock; clN = 0;
for ii = 2:5
    clN = floor(100*clN + cl(ii));
end
path = 'Best_Fit_Allee';
if exist(path, 'dir') ~= 7
    mkdir(path)
end

%% Variables
filename = '../Censored_Data/Control_data_censored.xlsx';
A = xlsread(filename);
A(A == 0) = NaN; %first row is not numbers
A(1,1) = 0; %fills in first cell as Day 0
num_mice = size(A,2)-3; % ignore columns with average and STD
numreps = 15;
Volume = {};
Days = {};
nump = zeros(30);

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
Best_r = zeros(num_mice,numreps + 1); Best_K = zeros(num_mice,numreps + 1);
Best_c = zeros(num_mice,numreps + 1); Best_y0 = zeros(num_mice,numreps + 1);
Best_SSE = zeros(num_mice,numreps + 1);
Figures = gobjects(num_mice,numreps + 1);

%% Sobol parameters
min_r = -0.2; max_r = .2; % change to 0.2
n_Sbpts = 1500000; % 1.5e6
n_SApts = 500000; % 5e5
SA_SSEs1 = {};
SA_SSEs2 = {};

%% Initialize parallel pool
numcores = num_mice;
parpool('local',numcores); % Call to open the distributed processing
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
        min_K = Volume_loop(1); max_K = 100000;
        min_c = 0; max_c = 10*Volume_loop(1);
        
        a = clock;
        rng(1000*a(5)*a(6));
        n_skip = randi([1000000, 10000000]); % 1e6 to 1e7
        n_leap = 0;
        
        %% Create sobol points and scale them
        % clearvars uniform_sobol_scaled
        uniform_sobol = sobolset(4,'Skip',n_skip,'Leap',n_leap)
        uniform_sobol = net(uniform_sobol,n_Sbpts);
        uniform_sobol_scaled = zeros(n_Sbpts,7);
        
        uniform_sobol_scaled(:,1) = (max_r-min_r)*uniform_sobol(:,1) + min_r;
        uniform_sobol_scaled(:,2) = (max_K-min_K)*uniform_sobol(:,2) + min_K;
        uniform_sobol_scaled(:,3) = (max_c-min_c)*uniform_sobol(:,3) + min_c;
        uniform_sobol_scaled(:,4) = (max_y0-min_y0)*uniform_sobol(:,4) + min_y0;
        SSE_QMC = [];
        
        for u = 1:n_Sbpts % n_pts  % Solve DE at all biologically-relevant Sobol points
            % need K>c always. Need y0<c when r<0. No restrictions when r>0
            pass_check = -1;
            if(uniform_sobol_scaled(u,2)>uniform_sobol_scaled(u,3)) % K > c
                if (uniform_sobol_scaled(u,1)>=0) %r > 0 is good-to-go
                    pass_check = 1;
                else % when r<0 also need y0<K
                    if(uniform_sobol_scaled(u,4)<uniform_sobol_scaled(u,2))
                        pass_check = 1;
                    end
                end
            end
            if pass_check == 1                
                %create first soln vector
                sb_soln = ode23s(@(t,x) allee(t,x,uniform_sobol_scaled(u,1),...
                    uniform_sobol_scaled(u,2),uniform_sobol_scaled(u,3)),...
                    Days_loop(1):Days_loop(end),uniform_sobol_scaled(u,4));
                
                newz = 0;
                for day = 1:length(Days_loop)
                    [c , index] = min(abs(sb_soln.x-Days_loop(day)));
                    error = (sb_soln.y(index) - Volume_loop(day))^2;
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
        newParam_r = []; newParam_K = []; newParam_c = [];
        y0 = []; SSE = []; SSEstep = [];
        
        i = 1;
        newParam_r(i) = uniform_sobol_scaled(index,1);
        newParam_K(i) = uniform_sobol_scaled(index,2);
        newParam_c(i) = uniform_sobol_scaled(index,3);
        y0(i)=  uniform_sobol_scaled(index,4);
        SSE(i) = c;

        %% Simulated annealing
        order_r = floor(log10(abs(newParam_r(1))))-1;
        order_K = floor(log10(newParam_K(1)))-1;
        order_c = floor(log10(newParam_c(1)))-1;
        order_y0 = floor(log10(y0(1)))-1;
        
        for j = 1:n_SApts
            %% Step 1: perturb parameters
            a = -1;
            b = 1;
            r = ((b-a)*(rand)+a)*10^(order_r);
            s = ((b-a)*(rand)+a)*10^(order_K);
            q = ((b-a)*(rand)+a)*order_c;
            v = ((b-a)*(rand)+a)*10^(order_y0);
            
            tempnewparam_r = newParam_r(i)+r;
            tempnewparam_K = newParam_K(i) + s;
            tempnewparam_c = newParam_c(i) + q;
            tempnewy0 = y0(i) + v;
            
            % no negative K, c or y0, and cannot have y0>K whenn r<0
            while ( (tempnewparam_K<0) || (tempnewparam_c<0) || (tempnewy0<0) || (tempnewparam_c > tempnewparam_K) || ((tempnewparam_r<0) && (tempnewy0>tempnewparam_K)) )
                r = ((b-a)*(rand)+a)*10^(order_r);
                s = ((b-a)*(rand)+a)*10^(order_K);
                q = ((b-a)*(rand)+a)*10^(order_c);
                v = ((b-a)*(rand)+a)*10^(order_y0);
                tempnewparam_r = newParam_r(i) + r;
                tempnewparam_K = newParam_K(i) + s;
                tempnewparam_c = newParam_c(i) + q;
                tempnewy0 = y0(i) + v;
            end
            
            %% Step 2: Solve ODEs at perturbed parameters
            sa_soln = ode23s(@(t,x) allee(t,x,tempnewparam_r,tempnewparam_K,...
                tempnewparam_c), Days_loop(1):Days_loop(end),tempnewy0);
            
            initialSSE = 0;
            newz = 0;
            for day = 1:length(Days_loop)
                [c , index] = min(abs(sa_soln.x-Days_loop(day)));
                error = initialSSE + (sa_soln.y(index) - Volume_loop(day))^2;
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
                SSEstep(i-1) = j;
                newParam_r(i) = tempnewparam_r;
                newParam_K(i) = tempnewparam_K;
                newParam_c(i) = tempnewparam_c;
                y0(i) = tempnewy0;
            end
        end
        
        fprintf("Mouse %d rep %d SA done\n",m,rep) % check SA done
        
        %% Optimal parameter set is the last accepted parameter set from simulated annealing
        Best_r(m,rep) = newParam_r(end);
        Best_K(m,rep) = newParam_K(end);
        Best_c(m,rep) = newParam_c(end);
        Best_SSE(m,rep) = SSE(end);
        Best_y0(m,rep) = y0(end);
        
        final_soln = ode23s(@(t,x) allee(t,x,newParam_r(end),newParam_K(end),...
            newParam_c(end)),Days_loop(1):Days_loop(end),y0(end));
        
        figure('visible','off');
        plot(final_soln.x,final_soln.y,'b','LineWidth',2); hold on;
        plot(Days_loop(:),Volume_loop(:),'ob','MarkerSize', 10,'MarkerFaceColor','b');
        hold off;
        xlabel('Days','FontSize',16);
        ylabel('Volume','FontSize',16);
        title(['Allee: Mouse ' num2str(m) 'rep ' num2str(rep)],'FontSize',18)
        
        fprintf('For Mouse %d rep %d, Best_r = %f, Best_K = %f, Best_c = %f\nBest_SSE = %f, Best_S0 = %f\n',...
            m, rep, Best_r(m,rep), Best_K(m,rep), Best_c(m,rep),...
            Best_SSE(m,rep),Best_y0(m,rep));
        
        Figures(m,rep) = gcf;
        SA_SSEs1{m,rep} = SSE(:);
        SA_SSEs2{m,rep} = SSEstep(:);
        fprintf("Mouse %d rep %d done\n\n",m,rep)
        
    end
    
    fprintf("Mouse %d all done\n",m) % check
    
    
end
delete(gcp('nocreate'))

BestRep = zeros(num_mice,1);
for mouse2 = 1:num_mice
    [minSSE, minrep] = min(Best_SSE(mouse2,1:numreps));
    Best_r(mouse2,numreps + 1) = Best_r(mouse2,minrep);
    Best_K(mouse2,numreps + 1) = Best_K(mouse2,minrep);
    Best_c(mouse2,numreps + 1) = Best_c(mouse2,minrep);
    Best_SSE(mouse2,numreps + 1) = Best_SSE(mouse2,minrep);
    Best_y0(mouse2,numreps + 1) = Best_y0(mouse2,minrep);
    SA_SSEs1{mouse2,numreps + 1} = SA_SSEs1{mouse2,minrep};
    SA_SSEs2{mouse2,numreps + 1} = SA_SSEs2{mouse2,minrep};
    BestRep(mouse2) = minrep;

    fname_fig = [path '/Mouse' num2str(mouse2) 'Rep' num2str(minrep)]; % how to put a number in string
    saveas(Figures(mouse2,minrep),[fname_fig,'.fig']);
    saveas(Figures(mouse2,minrep),[fname_fig,'.png']);
end

Final_r = Best_r(:,numreps + 1)'; Final_K = Best_K(:,numreps + 1)'; Final_c = Best_c(:,numreps + 1)';
Final_SSE = Best_SSE(:,numreps + 1)'; Final_y0 = Best_y0(:,numreps + 1)';

filename = [path '/BestParams.mat'];
save(filename,'Best_r','Best_K','Best_c',...
    'Best_SSE','Best_y0','BestRep','Final_r','Final_K',...
    'Final_c','Final_SSE','Final_y0',...
    'SA_SSEs1','SA_SSEs2')

