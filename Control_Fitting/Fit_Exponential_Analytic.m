%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% 7/7/20: Analytically fits an exponential growth/decay curve to each   %
% control mouse. Not there are no numerical aproximations here, so the  %
% results are exact. Data is transformed using a logarithm to be        %
% approximately linear, the best-fit line is found, and the slope and   %
% intercept of the line are transformed to get the base (b) and scaling %   
% factor (k) of the line precisely. This gives a best fit curve of the  %
% form: k*b^t.                                                          %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
path = 'Best_Fit_Exponential';
if exist(path, 'dir') ~= 7
    mkdir(path)
end

filename = '../Censored_Data/Control_data_censored.xlsx';
A = xlsread(filename);
A(A == 0) = NaN; %first row is not numbers
A(1,1) = 0; %fills in first cell as Day 0
num_mice = size(A,2)-3; % ignore columns with average and STD
scaling_factor = zeros(num_mice,1); 
base = zeros(num_mice,1); 
SSE = zeros(num_mice,1); 

for m = 2:num_mice+1  
	Volume_orig = A(:,m);
	Volume{m-1}(:) = rmmissing(Volume_orig); %removes any entry that contains missing data
    Days_orig = A(:,1);
	Days{m-1}(:) = Days_orig(1:length(Volume{m-1}(:))); 
    
    % Fit linear function to log of data
    log_volume = log(Volume{m-1}(:));  
    p = polyfit(Days{m-1}(:),log_volume,1); %m = p(1), b = p(2)
    scaling_factor(m-1,1) = exp(p(2));
    base(m-1,1) = exp(p(1)); 
    time = Days{m-1}(1):0.1:Days{m-1}(end); 
    exp_fit = scaling_factor(m-1)*base(m-1).^time; 
    figure; 
    plot(time,exp_fit,'LineWidth',2,'Color','b'); hold on;
    plot(Days{m-1}(:),Volume{m-1}(:),'ob','MarkerSize', 10,'MarkerFaceColor', 'blue'); 
    xlabel('Time (days)','FontSize',16);
    ylabel('Volume (mm^3)','FontSize',16); 
    hold off; 
	fname_fig = [path '/mouse' num2str(m-1)];        
	saveas(gcf,[fname_fig,'.fig'])
 	saveas(gcf,[fname_fig,'.png'])
    
    % Compute SSE
    mouse_SSE = 0; 
	for k = 1:size(Days{m-1}(:), 1) 
        [val, idx] = min( abs(Days{m-1}(k) - time) );
        error = (Volume{m-1}(k) - exp_fit(idx))^2;
        mouse_SSE = mouse_SSE + error;
    end
    SSE(m-1,1) = mouse_SSE; 
end

filename = [path '/Exponential_fit.mat']; 
save(filename, 'base', 'scaling_factor','SSE') 