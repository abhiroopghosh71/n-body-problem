clear
close all

%% Strong scaling
n = [1,2,4,8,16,32];
t = dlmread('strong_scaling_omp_1_N5000_T1000');
ideal = [];
for ii=1:6
    ideal = [ideal t(1,2)/(1000 * 2^(ii-1))];
end

figure
loglog(t(:,1),t(:,2)/1000, '-o'), hold on
loglog(t(:,1),ideal)
legend('Actual', 'Ideal')
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
xlabel('Number of processors')
ylabel('Average time per iteration')
grid
title('Strong Scaling')

%% Weak scaling
t = dlmread('weak_scaling_omp_1_1000_T1000');
ideal = [];
for ii=1:6
    ideal = [ideal t(1,2)/1000];
end

figure
loglog(t(:,1),t(:,2)/1000, '-o'), hold on
loglog(t(:,1),ideal)
ylim([1,100])
legend('Actual', 'Ideal')
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
xlabel('Number of processors')
ylabel('Average time per iteration')
grid
title('Weak Scaling')

%% Thread to thread
t = dlmread('th2th_N5000_T1000');
ideal = [];
for ii=1:6
    ideal = [ideal t(1,2)/(1000 * 2^(ii-1))];
end

figure
loglog(t(:,1),t(:,2)/1000, '-o'), hold on
loglog(t(:,1),ideal)
ylim([1,100])
legend('Actual', 'Ideal')
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
xlabel('Number of OpenMP threads')
ylabel('Average time per iteration')
grid
title('Thread-to-Thread Scaling')