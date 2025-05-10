
%%%The program calls EOLSO to estimate the values of parameters R1, R2, C1, and C2.
%%% This program was developed by Wuke Li  , Ying Xiong 
%%%
clc; clear; close all;
rng('default');

%% 读取 HPPC 的电流、电压，时间单位为 1s
load('data_plus.mat'); % 放电数据
load('OCV_SOC.mat');   % OCV-SOC 关系
U = data_plus(3, :);   % 电压
I = data_plus(2, :);   % 电流
SOC = data_plus(4, :); % SOC 真实值 - 安时法计算得到

%% 拟合 OCV-SOC 曲线，利用 8 次多项式拟合
x = OCV_SOC(2, :);
y = OCV_SOC(1, :);
p = polyfit(x, y, 8);
N1 = length(I);
ocv = zeros(1, N1);
ocv(1) = U(1); % 初始 ocv-soc 数据
for i = 2:N1
    ocv(i) = polyval(p, SOC(i));
end

%% 计算 R0
A = []; B = []; C = []; D = [];
for j = 1:length(I) - 1
    if I(j) == 0 && I(j + 1) < 0
        A = [A, j];     % 即将有电流的点
        B = [B, j + 1]; % 刚刚有电流的点
    end
end
for j = 1:length(I) - 1
    if I(j) < 0 && I(j + 1) == 0
        C = [C, j];     % 电流即将变为 0 的点
        D = [D, j + 1]; % 电流刚变为 0 的点
    end
end
M1 = U(A); % 即将有电流的点
M2 = U(B); % 刚刚有电流的点
N1_U = U(C); % 电流即将变为 0 的点
N2 = U(D);   % 电流刚变为 0 的点
Rohm = (M1(1:length(C)) - M2(1:length(C)) + N2 - N1_U) / 2 / 1;

%% 参数设置
timer1 = cputime;
Ts = 1;          % 系统采样时间
model_RC = 2;    % 模型选择
if model_RC == 1 % R int 模型
    % 参数上下界 [R1, C1]
    lb = [1e-03, 1e-02];
    ub = [5e-02, 1e+06];
    dim = 2;
elseif model_RC == 2
    lb = [0.01, 0.01, 1, 1000];  % 下界
    ub = [0.1, 0.1, 10000, 75000]; % 上界
    dim = 4;
elseif model_RC == 3
    lb = [1e-03, 1e-03, 1e-03, 1e-02, 1e-02, 1e-02];
    ub = [5e-02, 3e-02, 3e-02, 1e+06, 1e+06, 1e+06];
    dim = 6;
end

%% 区分放电和静置阶段的数据
M = {}; % 静置部分的时间、电流、电压、SOC 值
for i = 1:length(A) - 1
    idx = D(i):B(i + 1) - 1;
    I_ = I(idx);       % 电流
    U_ = U(idx);       % 电压
    SOC_ = SOC(idx);   % SOC
    t = 0:length(I_) - 1; % 时间
    ocv_ = ocv(idx);
    M{i} = [t; I_; U_; SOC_; ocv_]; % 静置部分的数据
end
N = {}; % 有电流激励部分的时间、电流、电压、SOC 值
for i = 1:length(D)
    idx = B(i):D(i) - 1;
    I_ = I(idx);       % 电流
    U_ = U(idx);       % 电压
    SOC_ = SOC(idx);   % SOC
    t = 0:length(I_) - 1; % 时间
    ocv_ = ocv(idx);
    N{i} = [t; I_; U_; SOC_; ocv_]; % 有电流激励部分的数据
end

%% 优化算法参数设置  
SearchAgents_no = 100; % 种群数量，建议选择 100  
Max_iter = 500;       % 迭代次数  
dim = 4;              % 参数维度  

% 初始化参数和收敛曲线存储 (参数数 x 算法数 x 数据段数)  
parameter_dis = zeros(dim, 6, 9); % 参数维度，算法数量（6），数据段数量  
Convergence_curves = zeros(Max_iter, 6, 9); % 迭代次数，算法数量，数据段数量  

% 初始化计算时间存储
computation_times = zeros(6, 9); % 算法数量, 数据段数量

for ZZ = 1:9  
    U_data = N{ZZ}(3, :)';  
    I_data = N{ZZ}(2, :)';  
    ocv_data = N{ZZ}(5, :);  
    R0_data = Rohm(:, ZZ);  

    % EOLSO
    tic; % 
    [Alpha_score_EOLSO, Alpha_pos_EOLSO, Convergence_curve_EOLSO, ~, ~] = EOLSO(...  
        SearchAgents_no, Max_iter, lb, ub, dim, U_data, I_data, ocv_data, R0_data, model_RC);  
    computation_times(1, ZZ) = toc; % 
    
%     % ESO  
%     tic; % 
%     [Alpha_score_PSO, Alpha_pos_PSO, Convergence_curve_PSO, ~, ~] = ESO_R0(...  
%         SearchAgents_no, Max_iter, lb, ub, dim, U_data, I_data, ocv_data, R0_data, model_RC);  
%     computation_times(2, ZZ) = toc; % 
%     
%     % SO 
%     tic; % 
%     [Alpha_score_SO, Alpha_pos_SO, Convergence_curve_SO, ~, ~] = SO_R0_claude(...  
%         SearchAgents_no, Max_iter, lb, ub, dim, U_data, I_data, ocv_data, R0_data, model_RC);  
%     computation_times(3, ZZ) = toc; % 
%     
%     % GJO  
%     tic; % 
%     [Alpha_score_GJO, Alpha_pos_GJO, Convergence_curve_GJO, ~, ~] = GJO_R0(...  
%         SearchAgents_no, Max_iter, lb, ub, dim, U_data, I_data, ocv_data, R0_data, model_RC);   
%     computation_times(4, ZZ) = toc; % 
%     
%     % HBA 
%     tic; % 
%     [Food_Score_HBA, Alpha_pos_HBA, Convergence_curve_HBA, ~, ~] = HBA_R0(...  
%         SearchAgents_no, Max_iter, lb, ub, dim, U_data, I_data, ocv_data, R0_data, model_RC);   
%     computation_times(5, ZZ) = toc; % 
%     
%     % GWO   
%     tic; % 
%     [Alpha_score_GWO, Alpha_pos_GWO, Convergence_curve_GWO] = GWO_R0(...  
%         SearchAgents_no, Max_iter, lb, ub, dim, U_data, I_data, ocv_data, R0_data, model_RC);  
%     computation_times(6, ZZ) = toc; % 
    
    % 存储所有算法的参数结果到同一个变量中  
    parameter_dis(:, 1, ZZ) = Alpha_pos_EOLSO'; % 第一列对应 DOBL_SO  
%     parameter_dis(:, 2, ZZ) = Alpha_pos_PSO';  % 第二列对应 PSO  
%     parameter_dis(:, 3, ZZ) = Alpha_pos_SO';   % 第三列对应 SO  
%     parameter_dis(:, 4, ZZ) = Alpha_pos_GJO';  % 第四列对应 GJO   
%     parameter_dis(:, 5, ZZ) = Alpha_pos_HBA';  % 第五列对应 HBA  
%     parameter_dis(:, 6, ZZ) = Alpha_pos_GWO';  % 第六列对应 GWO  

    % 存储收敛曲线  
    Convergence_curves(:, 1, ZZ) = Convergence_curve_EOLSO;  
%     Convergence_curves(:, 2, ZZ) = Convergence_curve_PSO;  
%     Convergence_curves(:, 3, ZZ) = Convergence_curve_SO;  
%     Convergence_curves(:, 4, ZZ) = Convergence_curve_GJO;  
%     Convergence_curves(:, 5, ZZ) = Convergence_curve_HBA;  
%     Convergence_curves(:, 6, ZZ) = Convergence_curve_GWO;  
end  

% 计算平均参数（取前8个数据段的平均）  
R0_mean = mean(Rohm(1:8));  

% 初始化平均参数存储  
R1 = zeros(6, 1); % 六个算法的 R1 平均值  
R2 = zeros(6, 1); % 六个算法的 R2 平均值  
C1 = zeros(6, 1); % 六个算法的 C1 平均值  
C2 = zeros(6, 1); % 六个算法的 C2 平均值  

% 计算平均参数  
for alg = 1:1  
    R1(alg) = mean(parameter_dis(1, alg, 1:8));  
    R2(alg) = mean(parameter_dis(2, alg, 1:8));  
    C1(alg) = mean(parameter_dis(3, alg, 1:8));  
    C2(alg) = mean(parameter_dis(4, alg, 1:8));  
end  

% 计算每个算法的平均计算时间（前8个数据段）
avg_computation_times = mean(computation_times(:, 1:8), 2);

% 将结果保存为结构体，方便比较  
results.R1 = R1;  
results.R2 = R2;  
results.C1 = C1;  
results.C2 = C2;  
results.R0 = R0_mean;  
results.computation_times = avg_computation_times;
save('results.mat', 'results');  

% 显示结果  
%algorithm_names = {'EOLSO', 'ESO', 'SO', 'GJO', 'HBA', 'GWO'}; 
algorithm_names = {'EOLSO'}; 
for alg = 1:1  
    fprintf('算法：%s\n', algorithm_names{alg});  
    fprintf('R0 = %.6f Ω\n', R0_mean);  
    fprintf('R1 = %.6f Ω\n', R1(alg));  
    fprintf('R2 = %.6f Ω\n', R2(alg));  
    fprintf('C1 = %.6f F\n', C1(alg));  
    fprintf('C2 = %.6f F\n', C2(alg));  
    fprintf('平均计算时间 = %.4f 秒\n', avg_computation_times(alg));
    fprintf('-----------------------------\n');  
end  

%% 计算并绘制平均收敛曲线  
% 计算每种算法在所有数据段上的平均收敛曲线  
average_convergence_EOLSO = mean(Convergence_curves(:, 1, 1:8), 3);  
% average_convergence_PSO = mean(Convergence_curves(:, 2, 1:8), 3);  
% average_convergence_SO = mean(Convergence_curves(:, 3, 1:8), 3);   
% average_convergence_GJO = mean(Convergence_curves(:, 4, 1:8), 3);  
% average_convergence_HBA = mean(Convergence_curves(:, 5, 1:8), 3);  
% average_convergence_GWO = mean(Convergence_curves(:, 6, 1:8), 3);  

% 创建迭代次数向量  
iterations = 1:Max_iter;  

% 绘制平均收敛曲线  
figure;  
plot(iterations, average_convergence_EOLSO, 'r-', 'LineWidth', 2); hold on;  
% plot(iterations, average_convergence_PSO, 'b--', 'LineWidth', 2);  
% plot(iterations, average_convergence_SO, 'g-.', 'LineWidth', 2);  
% plot(iterations, average_convergence_GJO, 'k:', 'LineWidth', 2);  
% plot(iterations, average_convergence_HBA, 'm-', 'LineWidth', 2);  
% plot(iterations, average_convergence_GWO, 'c--', 'LineWidth', 2);  

xlabel('Iterations', 'FontSize', 14);  
ylabel('Fitness', 'FontSize', 14);  
 
%legend('EOLSO','ESO','SO','GJO','HBA','GWO');  
grid on;  
