%%%######################%%%
%%%#Copyleft Liyang 2017#%%%
%%%######################%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%!!!ENCODING CLAIM!!!
%%%This code text use the GB2312 encoding form to
%%%support the Chinese.
%%%If there are some kind of display errors in your IDE, 
%%%please open the code-file without GB2312 in its name.
%%%Or you could use the 'iconv' command in *nix system
%%%to convert the text form by yourself.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%!!!程序功用声明!!!
%%%程序标题: 均匀介质中电磁波测井信号模拟
%%%编程时间: 2017.10.21
%%%预计时长: 5h
%%%目的: 使用电磁学基本的方程, 利用Matlab模拟测井信号的电压值
%%%功能: 数值计算两线圈之间的磁场
%%%注意:由于matlab版本问题, 有的matlab版本可能并不支持将脚本和
%%%     函数写在一起, 只需要将脚本和函数分开或存放, 或是将脚本编成
%%%     即可.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function func_loggin_siglelayer_GB2312()
  %#########################################
  %%%%%%%%%%%%计算前的准备%%%%%%%%%%%%%%%%%%%%
  %#########################################
  clear;   
  format long; %设置计算精度

  %#########################################
  %%%%%%%%%%%%基本物理学参量的计算%%%%%%%%%%%%%%
  %#########################################
  global I_T rho S M L mu omega k t W 

  I_T = 1;                               %发射线圈电流大小,A (T_T)
  rho = 10^(-3);   %rho_T = rho_R = rho  %发射/探测线圈半径,m
  S = pi * rho^2;                        %发射/探测线圈面积
  M = S * I_T;                           %发射线圈磁矩
 
  L = 0.8;                               %发射-探测线圈距离,m

  sigma = 1;                             %介质电导率,S/m
  mu_0 = 1.25663706144 * 10^(-6);        %真空磁导率
  mu_r = 1;                              %介质相对磁导率
  mu = mu_r * mu_0;                      %介质磁导率
  f = 20000;                             %信号频率,Hz
  omega = 2*pi*f;                        %信号角频率
  k = sqrt(1i * omega * mu * sigma);     %电磁波波矢

  t = [-0.9061798459, -0.5384693101,...
        0, ...
        0.5384693101, 0.9061798459];     %高斯积分中的零点
  W = [0.2369268851, 0.4786286705,...
       0.568888889, ...
       0.4786286705, 0.2369268851];      %高斯积分中的比重

  %#########################################
  %%%%%%%%%%%%%%%%%主函数%%%%%%%%%%%%%%%%%%%%
  %#########################################

  %%用两种方法求解磁场强度
  approx_H = cal_approx_H();
  precise_H = cal_precise_H();

  %%用两种方法求解电势大小
  %使用所得磁场求解电势
  approx_V_with_approx_H = cal_approx_V(approx_H);
  approx_V_with_precise_H = cal_approx_V(precise_H);
  %直接精确求解电势
  precise_V = cal_precise_V();

  %%处理所得结果为幅角和幅值的形式
  %幅值
  R_precise_V = abs(precise_V);
  R_approx_V_with_approx_H = abs(approx_V_with_approx_H);
  R_approx_V_with_precise_H = abs(approx_V_with_precise_H);
  %幅角
  theta_precise_V = angle(precise_V)/pi;
  theta_approx_V_with_approx_H = angle(approx_V_with_approx_H)/pi;
  theta_approx_V_with_precise_H = angle(approx_V_with_precise_H)/pi;
  %形式处理
  char_exp_precise_V = [num2str(R_precise_V),...
                       '*e^(iπ*',num2str(theta_precise_V),')'];
  char_exp_approx_V_with_approx_H = ...
                       [num2str(R_approx_V_with_approx_H),...
                       '*e^(iπ*',num2str(theta_approx_V_with_approx_H),')'];
  char_exp_approx_V_with_precise_H = ...
                       [num2str(R_approx_V_with_precise_H),...
                       '*e^(iπ*',num2str(theta_approx_V_with_precise_H),')'];

  %%输出所得结果
  final_print_char = ...
       ['-----------------------------------------------\n',           ...
        '=============均匀介质测井信号(V)==================\n',           ...
        '-----------------------------------------------\n',           ... 
        '\n',                                                          ...
        'V_approx(H_approx)   =',num2str(approx_V_with_approx_H),'\n', ...
        '                     =',char_exp_approx_V_with_approx_H,'\n', ...
        '\n',                                                          ...
        'V_approx(H_precise)  =',num2str(approx_V_with_precise_H),'\n',...
        '                     =',char_exp_approx_V_with_precise_H,'\n',...
        '\n',                                                          ...
        'V_precise            =',num2str(precise_V),'\n',              ...
        '                     =',char_exp_precise_V,'\n',              ...
        '(V)\n',                                                       ...
        '-----------------------------------------------','\n',        ...
        '-------------ProgramComplete!------------------','\n',        ...
        '-----------------------------------------------\n'];
  %将结果输出到屏幕
  fprintf(final_print_char); 
  
  %将结果输出到数据存储文件
  loggin_data_file_handle = fopen('loggin_sigle_layers_data.txt','w');
  fprintf(loggin_data_file_handle, final_print_char);
  fclose(loggin_data_file_handle);
end 
%#########################################
%%%%%%%%%%%%辅助功能函数%%%%%%%%%%%%%%%%%%%%
%#########################################
function approx_V = cal_approx_V(H)
%利用已知的磁场求近似的电势的函数

  global omega mu S;
  approx_V = 1i * omega * mu * S * H; 
end

%%%%%%%%

function precise_V = cal_precise_V()
%由精确公式, 不借助磁场大小, 直接得到电势的函数
%积分使用的是简单的矩形积分法

  global omega mu I_T rho L k
  
  kPhiIntervalSplitNum = 300;                         %定义在phi的离散取值点数
  
  phi = linspace(0, pi, kPhiIntervalSplitNum+1);      %获取取点坐标数据
  
  first_phi = phi;
  last_phi = phi;
  first_phi(end) = [];
  last_phi(1) = [];
  delta_phi = last_phi - first_phi;                   %求解坐标两两之间的距离
  
  R = sqrt(rho^2 + rho^2 - 2*rho*rho .* cos(phi) + L^2);  
  f_phi = exp(1i*k*R) ./ R .* cos(phi);               %被积函数的形式
  
  first_f_phi = f_phi;
  last_f_phi = f_phi;
  first_f_phi(end) = [];
  last_f_phi(1) = [];
  average_f_phi = (first_f_phi + last_f_phi) / 2;     %矩形区域的中心值
  
  inte_f_phi = sum(delta_phi .* average_f_phi);       %两个向量做内积
  
  precise_V = 1i * omega * mu * I_T * rho * rho * inte_f_phi;  
end

%%%%%%%%

function approx_H = cal_approx_H()
%计算近似的磁场大小
  global M k L;
  approx_H = (M/(2*pi*L^3)) * (1-1i*k*L) * exp(1i*k*L);
end

%%%%%%%%

function precise_H = cal_precise_H()
%计算精确的磁场大小
  global M;
  
  kPickIntervalPoints = 300;                           %积分中, 分区间的个数
  k_rho = get_interval_split(kPickIntervalPoints);     %分区间的分界坐标数据
  
  C = get_C(k_rho, kPickIntervalPoints);               %获取C的值
  
  inte_H = 0;                                          %初始化积分值
  n = 3;                                               %积分式子中的n值
  for C_index = 1:4                                    %l = C_index - 1
    for interval_num = 1:kPickIntervalPoints
      inte_H = inte_H + ...
               C(C_index, interval_num) * ...
               cal_inte_J_0_kn(n, C_index-1, ...
                               k_rho(interval_num), ...
                               k_rho(interval_num+1) );   
                                                       %累加计算积分的值
    end
  end
  
  precise_H = 1i*M/(4*pi) * inte_H;
end

%%%%%%%%

function inte_J_0_kn = cal_inte_J_0_kn(n, l, k_rho_j, k_rho_jplus1)
%使用高斯勒让德计算bessel函数乘以幂次n的积分
  global rho t W;
  
  x_j = k_rho_j * rho;
  x_jplus1 = k_rho_jplus1 * rho;                       %将k_rho_j转化为x_j
  
  x = (x_jplus1 + x_j)/2 + ((x_jplus1 - x_j)/2) * t;   %计算零点t对应的x的值
  g_x = besselj(0, x) .* (x.^(n+l));                     
  
  inte_J_0_xn = (x_jplus1 - x_j)/2 * sum(g_x .* W);    %计算简化后的bessel积分
  
  inte_J_0_kn = (1/rho^(n+l+1)) * inte_J_0_xn;         %输出最终积分值
end
  
%%%%%%%%

function C = get_C(k_rho, all_interval_num)
%用于计算系数C的函数
  global k L;
  %k_z = sqrt(k^2 - k_rho.^2);
  f_k_rho = exp(1i*sqrt(k^2 - k_rho.^2)*L) ./ sqrt(k^2 - k_rho.^2);        
                                                       %计算分界坐标上的函数值
  
  %使用对数组掐头去尾的方法计算
  %诸如M_j+1 - M_j的数值
  first_f_k_rho = f_k_rho;                              
  last_f_k_rho = f_k_rho;                              %将数组复制为两份
  first_f_k_rho(end) = [];                             %去掉尾部                      
  last_f_k_rho(1) = [];                                %去掉头部
  delta_f_k_rho = last_f_k_rho - first_f_k_rho;        %两数组相减
  
  first_k_rho = k_rho;
  last_k_rho = k_rho;
  first_k_rho(end) = [];
  last_k_rho(1) = [];
  h = last_k_rho - first_k_rho;                         %取点间隔h
  
  d_numberator_split = delta_f_k_rho ./ h;              %d表达式中分子中的一项
  
  first_d_numberator_split = d_numberator_split;
  last_d_numberator_split = d_numberator_split;
  first_d_numberator_split(end) = [];
  last_d_numberator_split(1) = [];
  d_numberator = 6 * (last_d_numberator_split - first_d_numberator_split); 
                                                        %d的分子
  
  first_h = h;
  last_h = h;
  first_h(end) = [];
  last_h(1) = [];
  h_add = first_h + last_h;                             %h错位相加所得的值
  
  d = d_numberator ./ h_add;                              
  d = [d d(end)];
  d = [d(1) d];                                         %d的值
  
  mu_Nplus1 = 1;
  mu = first_h ./ h_add;
  mu = [mu mu_Nplus1];                               %mu的值
  
  lambda_1 = 1;
  lambda = last_h ./ h_add;
  lambda = [lambda_1 lambda];                           %lambda的值
  
  M_matrix = diag(repmat(2,[1,all_interval_num+1])) + ...
                  diag(mu, -1) + diag(lambda, +1);     %M系数矩阵的值
  
  M_res = M_matrix \ d';
  M_res = M_res';                                       %M的解
  
  first_M_res = M_res;
  last_M_res = M_res;
  first_M_res(end) = [];
  last_M_res(1) = [];                                   %准备M的掐头去尾数组
  
  %计算各个C系数的值
  C_3 = (last_M_res - first_M_res) ./ (6 * h); 
  C_2 = (last_k_rho .* first_M_res - first_k_rho .* last_M_res) ./ (2 * h);
  C_1 = ((first_M_res - last_M_res)/6) .* h + ...
        (-last_k_rho.^2 .* first_M_res + first_k_rho.^2 .* last_M_res +...
        2*(last_f_k_rho - first_f_k_rho))./ (2 * h); 
  C_0 = (first_k_rho.*last_M_res - last_k_rho.*first_M_res)/6 .* h +...
        (last_k_rho.^3.*first_M_res - first_k_rho.^3.*last_M_res - ...
        6*(first_k_rho.*last_f_k_rho - last_k_rho.*first_f_k_rho)) ./ (6*h);
  %将C系数打包为一个2D数组 
  C = [];
  C = [C;C_0];
  C = [C;C_1];
  C = [C;C_2];
  C = [C;C_3];
end

%%%%%%%%

function interval_split = get_interval_split(split_num)
%获取分区间坐标点的函数
  k_1 = 10^(-6);                           %区间的下界
  k_max = 10^8;                            %区间的上界
  delta_k = log(k_max/k_1)/split_num;      %区间的增长因子
  
  interval_split = zeros(1,split_num+1);   %初始化坐标点数据存储数组
  for interval_split_index=1:split_num+1 
    interval_split(interval_split_index) = ...
          k_1*exp((interval_split_index-1)*delta_k);
  end 
end
