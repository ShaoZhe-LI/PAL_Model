%% 清理环境
clear; clc; close all;

%% 获取文件路径 fp_tmp 及 文件名称 fn0
[fp_tmp, fn0, ~] = fileparts(mfilename('fullpath')); % 当前 m 文件路径及名称

%% 数据路径选择及存储设置
data_set = 'Ideal source';
proc_save_dir = fullfile(fp_tmp, 'data', data_set);
if ~exist(proc_save_dir, 'dir')
    mkdir(proc_save_dir);
end

save_data = false; % 是否保存数据
save_fig = false; % 是否保存图片

%% 全局控制量
type_getG = 'Analytic'; % 计算Green函数谱的方法， 'Analytic' or 'FHT'
type_amp = 'Amp'; % 计算声压幅值的方法， 'Amp' or 'SPL'
type_field = 'Ultrasound'; % 声场类型， 'Ultrasound' or 'Audio'
type_profile = 'Vortex-m'; % 声源振动模式， 'Uniform' or 'Focus' or 'Vortex-m'






z_goal_phi = 1; % 绘制 z=1m 处的图像

x_range = [0.2, 1.5]; % x 方向显示范围
y_range = [-0.6, 0.6]; % y 方向显示范围
show_xylabel = false; % 是否显示 xy 轴标签


%% 参数设置
m1_v = [0 1 5 10];

a = 0.15;
v0 = 0.172;     % 默认0.172
c = 343;
rho = 1.21;
beta = 1.2;
fu = 42e3;  % 载波频率
fa = 2e3;   % 音频频率
f1=fu;
f2=fu+fa;
N_FHT = 16384 * 2;        % 径向点数，16384 * 2
delta = c/f2/3;         % z方向取点间隔，c/f2/3
delta = c/fa/6;         % z方向取点间隔，c/fa/3
rho_max = 1.2;
zu_max = 8;        % 18
za_max = 4+delta;   % 4
r_boundary = 0.24;   % xOy平面图像范围

results = cell(1, length(m1_v));  % N 是循环次数
max_results = zeros(1, length(m1_v));
min_results = zeros(1, length(m1_v));