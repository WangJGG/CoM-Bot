%%  清空环境变量
warning off             % 关闭报警信息
close all               % 关闭开启的图窗
clear                   % 清空变量
clc                     % 清空命令行

%%  导入数据
% 读取CSV文件
data = readtable('..\static\uploads\nor.csv');

% 获取数据的行数  
numRows = height(data);  
  
% 获取数据的列数  
numColumns = width(data);  
  
% 初始化一个新的矩阵来存储处理后的数据  
res = zeros(0, 100 * numColumns); % 根据列数调整大小  
  
% 重新组织数据  
for i = 1:100:numRows  
    % 提取300行数据  
    rowData = data(i:i+99, :);  
      
    % 初始化新行向量  
    newRow = [];  
      
    % 遍历每一列，并将列数据转换为行向量后添加到newRow中  
    for j = 1:numColumns  
        columnData = rowData{:, j}; % 提取第j列的数据  
        newRow = [newRow; columnData]; % 将列数据转换为行向量并添加到newRow中  
    end  
      
    res(end+1, :) = newRow;  
end

% 保存处理后的数据到新的CSV文件
%writetable(array2table(reshape(newData, [], numel(data.Variables))), '2.csv');;
%%  划分训练集和测试集
P_test = res';
N = size(P_test, 2);

%%  数据归一化
[P_test, ps_input] = mapminmax(P_test, 0, 1);
%%  数据平铺
% 将数据平铺成1维数据只是一种处理方式
% 也可以平铺成2维数据，以及3维数据，需要修改对应模型结构
% 但是应该始终和输入层数据结构保持一致
P_test  =  double(reshape(P_test , 3600, 1, 1, N));

%%  数据格式转换

for i = 1 : N
    p_test{i, 1} = P_test( :, :, 1, i);
end

%%  仿真预测
load('net.mat');
t_sim2 = predict(net, p_test ); 
T_sim2 = vec2ind(t_sim2');
  
% 指定文件夹路径和文件名  
folderPath = '..\Result'; 
fileName = 'output_num.txt';   
  
% 拼接完整的文件路径  
fullFilePath = fullfile(folderPath, fileName);  
  
fid = fopen(fullFilePath, 'w');  
if fid ~= -1  
    fprintf(fid, '%.6f\n', T_sim2);   
    fclose(fid);   
else  
    error('无法打开文件以写入。');  
end
