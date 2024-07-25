% 文件选择和加载
[FileName1, PathName1] = uigetfile('*.txt', '选择监测设备数据的文本文件');
dataFile1 = fullfile(PathName1, FileName1);
C = load(dataFile1);

% 常量定义
c = 340;  % 声速(m/s)
lat_dist = 111.263;  % 纬度间距离(km/deg)
lon_dist = 97.304;  % 经度间距离(km/deg)

% 计算距离函数
function d = dist(x1, y1, z1, x2, y2, z2)
    d = sqrt((x1-x2)^2 * lon_dist^2 + (y1-y2)^2 * lat_dist^2 + (z1-z2)^2);
end

function obj = objectiveFunction(params, C, c, distFunc)
    N = size(C, 1); % 监测设备数量
    M = size(C, 2) - 3; % 残骸数量
    debrisParams = reshape(params(1:4*M), 4, M); % 重塑为4行M列，每列一个残骸参数
    assignment = params(4*M+1:end); % 每个观测时间对应的残骸编号
    
    % 检查 assignment 的长度是否正确
    if length(assignment) ~= M
        error('分配参数的数量必须与残骸数量相同');
    end
    
    obj = 0;
    for i = 1:N
        for j = 1:M
            k = assignment(j); % 对应的残骸编号
            % 检查残骸编号是否在有效范围内
            if k < 1 || k > size(debrisParams, 2)
                error('残骸编号 %d 不在有效范围内 (1-%d) for observation %d', k, size(debrisParams, 2), i);
            end
            % 假设 debrisParams 的第四列是残骸的音爆时间
            tk = debrisParams(4, k); % 残骸的音爆时间
            % 计算观测时间与残骸音爆时间的差
            observedTime = C(i, 4 + j); % 假设第四列之后是观测时间
            obj = obj + (observedTime - tk)^2; % 目标函数计算平方时间差
        end
    end
    
    % 输出 assignment 和 debrisParams 以供调试
    disp('assignment:');
    disp(assignment);
    disp('debrisParams:');
    disp(debrisParams);
end

% 假设 C 已经被正确加载，并且其尺寸是 [N, 3+M]
N = size(C, 1); % 监测设备数量
M = size(C, 2) - 3; % 残骸数量

% 确定参数数量
nvars = 4 * M + M; % 残骸参数（每个4个）和分配参数的总数
% 创建一个初始向量，然后通过重复来匹配 nvars 的长度
initial_vector_debris = [-Inf, -Inf, -Inf, -Inf]; % 残骸参数的下界
initial_vector_assignment = [1]; % 分配参数的下界
repeats_debris = max(1, ceil((4 * M) / length(initial_vector_debris))); % 确保至少重复一次
repeats_assignment = max(1, ceil(M / length(initial_vector_assignment))); % 确保至少重复一次

lb_vector_debris = repmat(initial_vector_debris, 1, repeats_debris);
ub_vector_debris = repmat([Inf, Inf, Inf, Inf], 1, repeats_debris); % 残骸参数的上界
lb_vector_assignment = repmat(initial_vector_assignment, 1, repeats_assignment);
ub_vector_assignment = repmat(M, 1, repeats_assignment); % 分配参数的上界

% 合并界限向量
lb_vector = [lb_vector_debris, lb_vector_assignment];
ub_vector = [ub_vector_debris, ub_vector_assignment];

% 确保界限向量是行向量
lb = lb_vector(:)';
ub = ub_vector(:)';

% 设置粒子群算法选项
options = optimoptions('particleswarm', 'MaxIterations', 1000, 'Display', 'iter');
objFuncWrapper = @(params) objectiveFunction(params, C, c, @dist);

% 执行粒子群算法
[x, fval] = particleswarm(objFuncWrapper, nvars, lb, ub, options);

% 提取最优解
debrisParamsOptimal = reshape(x(1:4*M), 4, M);
assignmentOptimal = x(4*M+1:end);

% 重新排序C中的音爆抵达
C_sorted = C;
for i = 1:size(C, 1)
    for j = 1:M
        C_sorted(i, 4+j) = C(i, 4 + assignmentOptimal(j)); % 使用正确的索引
    end
end

% 输出最优的重新排序后的列表
disp('最优的重新排序后的列表：');
disp(C_sorted);
