clear; clc; close all;

%% =========================
% 参数
%% =========================
P0 = [0, 0];        % 起点
L0 = 12;            % 初始边长
dL = 0.45;          % 每一段长度递减量
Nseg = 24;          % 总段数
lw = 2.5;           % 线宽

bg_color = [0.92 0.92 0.92];
line_color = 'k';

% 三个方向（形成60度夹角的等边三角形边方向）
angles_deg = [0, 120, -120];

%% =========================
% 生成三角螺旋折线
%% =========================
pts = zeros(Nseg+1, 2);
pts(1,:) = P0;

for k = 1:Nseg
    Lk = L0 - (k-1)*dL;
    if Lk <= 0
        pts = pts(1:k,:);
        break;
    end
    
    ang = angles_deg(mod(k-1,3)+1);
    dir = [cosd(ang), sind(ang)];
    
    pts(k+1,:) = pts(k,:) + Lk * dir;
end

%% =========================
% 绘图
%% =========================
figure('Color','w');
ax = axes;
ax.Color = bg_color;
hold on; axis equal;

plot(pts(:,1), pts(:,2), '-', ...
    'Color', line_color, ...
    'LineWidth', lw);

% 可选标注
text(max(pts(:,1))*0.62, max(pts(:,2))*0.9, '(3,-1)', ...
    'FontSize', 28, ...
    'FontName', 'Times New Roman', ...
    'Color', 'k');

% 外观
margin = 1.0;
xlim([min(pts(:,1))-margin, max(pts(:,1))+margin]);
ylim([min(pts(:,2))-margin, max(pts(:,2))+margin]);

set(gca, 'XTick', [], 'YTick', []);
box on;
set(gca, 'LineWidth', 3);