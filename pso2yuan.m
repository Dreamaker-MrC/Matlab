%% Mr.C
%% I. 清空环境
clc
clear
 
%% II. 绘制目标函数曲线
figure
[x, y] = meshgrid(-5: 0.1: 5, -5: 0.1: 5);
z = x.^2 + y.^2 - 10*cos(2*pi*x) - 10*cos(2*pi*y) + 20;
mesh(x,y,z)
hold on
 
 
%% III. 参数初始化
c1 = 1.49445;          %权重参数c1
c2 = 1.49445;          %权重参数c2
 
maxgen = 100;                                                  % 进化次数 
sizepop = 10;                                                  % 种群规模
dimension = 2;                                                  % 这里因为是二元函数的求解，即二维，故列数为2
 
% 速度的边界
Vmax = 1;
Vmin = -1;
% 种群的边界
popmax = 5;
popmin = -5;
% 用于计算惯性权重，经验值
ws = 0.9;
we = 0.4;
%用于计算直觉模糊熵
Tx=ones(sizepop,1)           %全局计数器
Ti=ones(sizepop,1)          %粒子计数器
E=zeros(maxgen,1);                %每代的直觉模糊熵
L=ones(sizepop,sizepop)         %粒子距离中心点的距离
R=ones(sizepop,1)            %半径 
u=ones(sizepop,1)             %隶属度
pai=ones(sizepop,1)           %非隶属度
r=ones(sizepop,1)           %犹豫度
w=zeros(maxgen,1); 
% 给矩阵预分配内存
pop = zeros(sizepop, dimension);
V = zeros(sizepop, dimension);
fitness = zeros(sizepop, 1);
yy = zeros(maxgen);
%w = zeros(maxgen);
%% IV. 产生初始粒子和速度
for i = 1: sizepop
    % 随机产生一个种群
    pop(i, :) = 5 * rands(1, 2);                                % 初始种群
    V(i, :) = rands(1, 2);                                      % 初始化速度
    % 计算适应度
    fitness(i) = fun(pop(i, :));
end
 h=scatter3(pop(:,1),pop(:,2),fitness,'r');                     %绘制初始位置,红色的点是初始点
%% V. 初始化Personal best和Global best
[bestfitness, bestindex] = max(fitness);
gbest = pop(bestindex, :);                                      % Global best
pbest = pop;                                                    % 个体最佳
fitnesspbest = fitness;                                         % 个体最佳适应度值
fitnessgbest = bestfitness;                                     % 全局最佳适应度值
 
%% VI. 迭代寻优
for i = 1: maxgen
    
    for j=1:sizepop
        for k=1:sizepop
            L(k,j)=norm(pop(j)-pbest(k));
            M=max(L,[],2)
        end
     
        K=0.5*rand
        R(j)=K*M(j)
        
        for k=1:sizepop
            if R(j)>L(k,j)
            Ti(j)=Ti(j)+1    %粒子计数器加一
            else
            Tx(j)=Tx(j)+1    %全局计数器加一
            end
        end
       
        u(j)=Ti(j)/sizepop
        pai(j)=Tx(j)/sizepop
        r(j)=1-u(j)-pai(j)
        
        E(i)=E(i)+(1/sizepop)*(1-(u(j)-r(j)).^2+2*pai(j))/(2-(u(j)-r(j)).^2+pai(j)) %直觉模糊熵表达式
        
    end
    
    w(i)=ws-(ws-we)*(i/maxgen)*(1-E(i))    % 根据直觉模糊熵更新惯性因子w
%     w(i) = ws - (ws - we) * (i / maxgen); 
       
    for j = 1: sizepop
        % 速度更新
        V(j, :) = w(i)*V(j, :) + c1*rand*(pbest(j, :) - pop(j, :)) + c2*rand*(gbest - pop(j, :));
        for k = 1: dimension
            if V(j, k) > Vmax
                V(j, k) = Vmax;
            end
            if V(j, k) < Vmin
                    V(j, k) = Vmin;
            end
        end
         
        % 种群更新(位置更新)
        pop(j, :) = pop(j, :) + V(j, :);
        for k = 1: dimension
            if pop(j, k) > popmax
                pop(j, k) = popmax;
            end
            if pop(j, k) < popmin
                pop(j, k) = popmin;
            end
        end
        % 适应度值更新
        fitness(j) = fun(pop(j, :));
    end
     
    for j = 1: sizepop   
        % 个体最优更新
        if fitness(j) > fitnesspbest(j)
            pbest(j, :) = pop(j, :);
            fitnesspbest(j) = fitness(j);
        end
         
        %群体最优更新
        if fitness(j) > fitnessgbest
            gbest = pop(j, :);
            fitnessgbest = fitness(j);
        end
    end
    yy(i) = fitnessgbest;                                           % 记录每次迭代完毕的群体最优解
end
%% VII.输出结果
[fitnessgbest, gbest]
plot3(gbest(1), gbest(2), fitnessgbest, 'bo','linewidth', 1.5)
 
figure
plot(yy)
title('最优个体适应度', 'fontsize', 12);
xlabel('进化代数', 'fontsize', 12);
ylabel('适应度', 'fontsize', 12);

%% VIII待优化函数
function y = fun(x)
y = x(1).^2 + x(2).^2 - 10*cos(2*pi*x(1)) - 10*cos(2*pi*x(2)) + 20;
end