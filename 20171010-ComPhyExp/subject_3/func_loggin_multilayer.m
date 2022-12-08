%%%######################%%%
%%%#Copyleft 2017 Liyang#%%%
%%%######################%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%!!!ENCODING CLAIM!!!
%%%This code text use the UTF-8 encoding form to 
%%%support the Chinese.
%%%If there are some kind of display errors in your IDE, 
%%%please open the code-file with GB2312 in its name.
%%%Or you could use the 'iconv' command in *nix system
%%%to convert the text form by yourself.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%!!!程序功用声明!!!
%%%程序名称: 测井信号模拟系统
%%%编程时间: 2017.11.4
%%%预计用时: 6h
%%%目的:测井计算
%%%功能:利用matlab计算测井信号强度
%%%注意:由于matlab版本问题, 有的matlab版本可能并不支持将脚本和
%%%    函数写在一起, 只需要将脚本和函数分开或存放, 或是将脚本编成
%%%    函数即可.
%%%    本代码选取了第二种解决方式.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function func_loggin_multilayer()
  %====================================================
  %%%%%%%%%%%%%%%%%%计算前的准备%%%%%%%%%%%%%%%%%%%%%%%%%
  %====================================================
  clear; clc;  
  format long;                              %设置计算精度
  
  disp('########################');
  disp('******测井理论计算程式*****');
  disp('########################');
  disp('Version 1.3.1');
  disp('last modified 2017-10-11')
  disp('Copyleft 2017 Liyang');
  disp(' ');
  %====================================================
  %%%%%%%%%%%%%%%%%%计算模式的选择%%%%%%%%%%%%%%%%%%%%%%%%
  %====================================================
  %选择测试模式开关, 1:开启, 0:关闭
  global kOpenTestMode
  choice_cal_method = input('您是否要进入测试模式(Y/N):', 's');
  if strcmpi('Y', choice_cal_method) || strcmpi('y', choice_cal_method)
    kOpenTestMode = 1;
    disp('您选择进入了>>测试模式<<!');
  else
    kOpenTestMode = 0; 
    disp('您选择进入了>>计算模式<<!');
  end
  disp(' ');
  
  %计算模式(采用截断法或直接求解法求V)
  global kUseFMinusI
  choice_cal_method = input('您是否要使用[F-I]方法计算(Y/N):', 's');
  if strcmpi('N', choice_cal_method) || strcmpi('n', choice_cal_method)
    kUseFMinusI = 0;
    disp('您选择进入了>>[F]计算模式<<!');
  else
    kUseFMinusI = 1; 
    disp('您选择进入了>>[F-I]计算模式<<!');
  end
  disp(' ');
  
  %选择测井模型的层数数据
  kTotalLayerNum = 3;                       
  data_list_char = ['1. 一层模型';'3. 三层模型'; '5. 五层模型'];
  data_list_number = [1, 3, 5];
  
  fprintf('您可以调用的数据标签有:\n');
  sizeof_data_list_char = size(data_list_char);
  for print_step = 1:sizeof_data_list_char(1)
    fprintf('-->');
    fprintf('%c', data_list_char(print_step,:));
    fprintf('\n');
  end
  fprintf('\n');
  while 1
    kTotalLayerNum = input('请输入要计算的模型的数据标签对应的数字:');
    
    if ismember(kTotalLayerNum, data_list_number)
      break;
    end
    
    disp('!!!输入的数据不在列表当中, 请重新选择!!!');
  end
  
  tic;                                        %程序计时器, 计算计时开始
  disp('Timer Begin...');
  
  %====================================================
  %%%%%%%%%%%%%%%%%%基本物理参数的定义%%%%%%%%%%%%%%%%%%%%
  %====================================================
  %%%-----------基本的公用参数------------%%%
  global omega I_T z_R rho_T rho_R t W 
  global kPickIntervalNumVref kPickIntervalNumVdir kPickIntervalNumEphi
  f = 2000000;                                %信号频率,Hz
  omega = 2*pi*f;                             %信号角频率
  I_T = 1;                                    %发射线圈电流大小,A (T_T)
  z_R = 0.762;                                %发射-探测线圈距离,m
  
  %%%使rho_T与rho_R有轻微不同, 可使E_phi快速收敛%%%
  %%%进而可以用直接法(不利用[F-I]替换法)求V%%%
  rho_T = 0.1143;                             %接受线圈的半径,m
  rho_R = rho_T+0.0001;                       %探测线圈的半径,m

  kPickIntervalNumVref = 1000;                %V_ref积分取点个数
  kPickIntervalNumVdir = 300;                 %V_dir积分取点个数
  kPickIntervalNumEphi = 1000;                %测试Ephi绘图的取点数

  t = [-0.9061798459, -0.5384693101,...
        0, ...
        0.5384693101, 0.9061798459];          %L5高斯积分中的零点
  W = [0.2369268851, 0.4786286705,...
       0.568888889, ...
       0.4786286705, 0.2369268851];           %L5高斯积分中的比重
  %%%+++++++++基本的公用参数结束++++++++++%%%
  
  %%%--------测井层参数数据库-------------%%%
  global mu rho_layer sigma epsilon_star k
  
  mu_0 = 1.25663706144 * 10^(-6);             %真空介质磁导率
  epsilon_0 =  8.854187817*10^(-12);          %真空介电常数
  
  %------5层数据的电磁学性质参数---------
  if kTotalLayerNum == 5     
    mu_r = [1, 1, 1, 1, 1];                   %由里到外各层的相对磁导率
    rho_layer = [0.1016, 0.1118, 0.1270, 0.3048];
                                              %各层的分界半径值 
    epsilon_r = [1, 1, 80, 1, 1];             %由里到外各层的相对介电常数
    sigma = [10^8, 10^(-5), 10, 1, 1];        %由里到外各层的电导率, S/m
  
  %------3层数据的电磁学性质参数---------
  elseif kTotalLayerNum == 3
    mu_r = [1, 1, 1];                         %由里到外各层的相对磁导率
    rho_layer = [0.1016, 0.1270];             %各层的分界半径值
    epsilon_r = [1, 80, 1];                   %由里到外各层的相对介电常数
    sigma = [10^8, 10, 1];                    %由里到外各层的电导率, S/m
  %------均匀介质的电学性质参数----------
  elseif kTotalLayerNum == 1
    mu_r = [1, 1, 1];                         %由里到外各层的相对磁导率
    rho_layer = [0.1016, 0.1270];             %各层的分界半径值
    epsilon_r = [1, 1, 1];                    %由里到外各层的相对介电常数
    sigma = [1, 1, 1];                        %由里到外各层的电导率, S/m
  %----------------------------------
  %----------------------------------
  end
  %+++++++++++数据选项卡结束++++++++++++
  
  mu = mu_0 * mu_r;                           %由里到外各层的绝对磁导率
  epsilon = epsilon_0 * epsilon_r;            %由里到外各层的绝对介电常数
  epsilon_star = epsilon+1i*(sigma/omega);    %由里到外各层的复介电常数
  k = sqrt(omega^2*mu.*epsilon_star);         %电磁波波矢数值
  %%%++++++++测井层参数数据库结束+++++++++++++%%%
  
  %%%------------另外的特殊参数--------------%%%
  global m total_layer_num
  total_layer_num = length(mu_r);             %获取总共的层数
  m = get_device_layer();                     %获取探测器在哪一层介质当中
  %%%++++++++++另外的特殊参数结束+++++++++++++%%%
  
  %====================================================
  %%%%%%%%%%%%%%%%%%%主函数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %====================================================
 
  if kOpenTestMode
    %%%进入测试模块
    clc;
    disp('================================');
    disp('!!!您目前处于数据计算正确性测试模式!!!');
    disp('================================');
    disp('正在计算初始指定的sigma下的电压值...');
    
    if kUseFMinusI                              %如果使用了[F-I]求解E_phi, 则需要加V_dir
      test_V_all = cal_V_ref() + cal_V_dir(); 
    else                                        %如果使用[F]求解E_phi则不需要加V_dir
      test_V_all = cal_V_ref();
    end
    disp('V_special = ');
    disp(test_V_all);
    disp('(V)');
    disp(' ');
    
    %------------------------测试E_phi的模块--------------------
    disp('正在生成E_phi的图像...');
    
    %初始化测试k_z数组和E_phi数组
    k_z_test = get_interval_split(10^(-10), 10^15, kPickIntervalNumEphi-1);
    E_ref_array = zeros(1,kPickIntervalNumEphi);
    E_all_array = zeros(1,kPickIntervalNumEphi);
    
    %计算E_phi数据
    for interval_num = 1:kPickIntervalNumEphi
      %计算(F-I)对应的电场值
      E_ref_array(interval_num) = cal_E_phi_kz(k_z_test(interval_num), 1);
      %计算(F)对应的电场值
      E_all_array(interval_num) = cal_E_phi_kz(k_z_test(interval_num), 0);
    end
    
    %绘制图像
    ref_real_pic = loglog(k_z_test, abs(real(E_ref_array)), '--b','LineWidth',2);
    hold on;
    ref_imag_pic = loglog(k_z_test, abs(imag(E_ref_array)), '--r','LineWidth',2);
    hold on;
    all_real_pic = loglog(k_z_test, abs(real(E_all_array)), '-b','LineWidth',2);
    hold on;
    all_imag_pic = loglog(k_z_test, abs(imag(E_all_array)), '-r','LineWidth',2);
    grid on;
    axis([10^(-10),10^15,10^(-100),1]);
    scrsz = get(0,'ScreenSize');  
    set(gcf,'position', ...
        0.6*scrsz+0.1*[scrsz(3), scrsz(3),scrsz(4),scrsz(4)]); 
                                              %设置窗口大小
    legend([ref_imag_pic,ref_real_pic,all_imag_pic,all_real_pic],...
            '|Ref Image|',' |Ref Real|', '|Total Image|','|Total Real|');
                                              %设置图例
    xlabel('k_z');                     
    ylabel('E_{phi}(k)');
    title('Test for E_{phi}');
    saveas(gcf,['loggin_', num2str(kTotalLayerNum), '_layers_Ephi.jpg'])
                                              %存储图像文件
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
     disp('E_phi的图像生成完毕!');
     disp(' ');
     disp('数据测试结束!'); 
     
  else 
    %%%进入正式计算模块
    %--------------进度条显示模块.部分(上).-----------------------------
    all_pixel_num = 60;
    char_progress_bar = ['[',...
                         repmat('#', [1,0]),...
                         repmat(' ', [1,all_pixel_num]),...
                         ']',...
                         '0.0%'];
    clc;
    disp(['Calculate the ', num2str(kTotalLayerNum), ' Layers Model!']);
    disp(' ');
    disp('Program in Progressing, Plseae Wait...');  
    disp(char_progress_bar);
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    %---------------sigma末尾元素数值的改变------------------------------
    kPickSigmaFNum = 79;                          %最外层sigma变化取点个数(+1)
    kSigmaFLowerLimit = 10^(-3);                  %最外层sigma取值下限
    kSigmaFUpperLimit = 10;                       %最外层sigma取值上限

    sigma_end = get_interval_split...
               (kSigmaFLowerLimit, kSigmaFUpperLimit , kPickSigmaFNum);         
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    test_V_all = zeros(1,kPickSigmaFNum+1);                %初始化V存储数组
    
    %%循环求解若干个绘图点
    for step_sigma = 1:kPickSigmaFNum+1  
      %更新各个参数
      sigma(end) = sigma_end(step_sigma);         %由里到外各层的电导率, S/m
      epsilon_star = epsilon+1i*(sigma/omega);    %由里到外各层的复介电常数
      k = sqrt(omega^2*mu.*epsilon_star);         %电磁波波矢数值
      
      %###############核心计算组件#######################################
      if kUseFMinusI                              %如果使用了[F-I]求解E_phi, 则需要加V_dir
        V_all = cal_V_ref() + cal_V_dir();
      else                                        %如果使用[F]求解E_phi则不需要加V_dir
        V_all = cal_V_ref();
      end
      test_V_all(step_sigma) = V_all;             %计算接受线圈的电势
      %################################################################
      
      %--------------进度条显示模块.部分(下).------------------------------
      progress_percent = round(step_sigma / (kPickSigmaFNum+1) * 1000);
      progress_pixel_num = round(step_sigma / (kPickSigmaFNum+1)*all_pixel_num);
      char_progress_bar = ['[',...
                           repmat('#', [1,progress_pixel_num]),...
                           repmat(' ', [1,all_pixel_num-progress_pixel_num]),...
                           ']',...
                           num2str(floor(progress_percent/10)),'.',...
                           num2str(mod(progress_percent,10)),'%'];
      clc;
      disp(['Calculate the ', num2str(kTotalLayerNum), ' Layers Model!']);
      disp(' ');
      disp('Program in Progressing, Plseae Wait...');
      disp(char_progress_bar);
      %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end

    disp('Calculation Complete!');                %计算完成提示               

    %-----------数据曲线绘制模块--------------------------------------------
    imag_pic = semilogx(sigma_end, imag(test_V_all), '-r','LineWidth',2);
                                              %绘制虚部
    hold on;   
    real_pic = semilogx(sigma_end, real(test_V_all), '-b','LineWidth',2);
                                              %绘制实部
    grid on;
    xlabel('sigma(S/m)');                     
    ylabel('Voltage(V)');
    legend([imag_pic,real_pic],'Imag','Real');%设置图例
    title(['Loggin: ', num2str(kTotalLayerNum), ' Layers Model']);
                                              %设置曲线标题
    scrsz = get(0,'ScreenSize');              %获取屏幕分辨率
    set(gcf,'position', ...
        0.6*scrsz+0.1*[scrsz(3), scrsz(3),scrsz(4),scrsz(4)]);           
                                              %设置绘图窗体全屏
    saveas(gcf,['loggin_', num2str(kTotalLayerNum), '_layers_V.jpg'])
                                              %存储图像文件
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    %----------数据存储文件设置---------------------------------------------
    char_file_name = ['loggin_', num2str(kTotalLayerNum), '_layers_data.txt'];
    data_file_handle = fopen(char_file_name,'w');

    fprintf(data_file_handle, '############说明开始##############\n');
    fprintf(data_file_handle, '#本文档是func_loggin_multilayer.m生成的数据文档\n');
    fprintf(data_file_handle, '#您可以根据自己的需求提取数据\n');
    fprintf(data_file_handle, '#不同数据的标识符为 <==:datakind:==> \n');
    fprintf(data_file_handle, '#不同的datakind代表不同物理量\n');
    fprintf(data_file_handle, '#您需要自己编程提取本文档中的数据\n');
    fprintf(data_file_handle, '####copyleft---liyang---2017####\n');
    fprintf(data_file_handle, '############说明结束##############\n\n\n\n');
    fprintf(data_file_handle, '>>>DATA_READ_POINT<<<\n');
    fprintf(data_file_handle, '<==:sigma:==>\n');  
    fprintf(data_file_handle, '%e  ', sigma_end);
    fprintf(data_file_handle, '\n\n<==:V_image:==>\n');  
    fprintf(data_file_handle, '%e  ', imag(test_V_all));
    fprintf(data_file_handle, '\n\n<==:V_real:==>\n');
    fprintf(data_file_handle, '%e  ', real(test_V_all));
    
    fclose(data_file_handle);
    
    disp('Data Saver Complete!');
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end
  
  use_time = toc;                            %程序计时器, 计时结束
  disp('Timer End.');
  disp(['Total Time: ', num2str(use_time),'s']);
  disp(' ');
  disp('Program Complete!');                 %程序完成提示
end%主函数结束

%=======================================================
%%%%%%%%%%%%%%%%%%辅助函数定义%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=======================================================
%%%----------------直接计算V_dir的数值------------------%%%
function V_dir = cal_V_dir()
%由精确公式, 直接得到V_dir值的函数
%积分使用的是简单的矩形积分法

  global omega mu m I_T rho_R rho_T z_R sigma kPickIntervalNumVdir
  
  phi = linspace(0, pi, kPickIntervalNumVdir+1);      %获取取点坐标数据
  
  first_phi = phi;
  last_phi = phi;
  first_phi(end) = [];
  last_phi(1) = [];
  delta_phi = last_phi - first_phi;                   %求解坐标两两之间的距离
  
  R = sqrt(rho_T^2 + rho_R^2 - 2*rho_R*rho_T .* cos(phi) + z_R^2);  
  k_dir = sqrt(1i * omega * mu(m) * sigma(m));
  f_phi = exp(1i*k_dir*R) ./ R .* cos(phi);           %被积函数的形式
  
  first_f_phi = f_phi;
  last_f_phi = f_phi;
  first_f_phi(end) = [];
  last_f_phi(1) = [];
  average_f_phi = (first_f_phi + last_f_phi) / 2;     %矩形区域的中心值
  
  inte_f_phi = sum(delta_phi .* average_f_phi);       %两个向量做内积
  
  V_dir = 1i * omega * mu(m) * I_T * rho_R * rho_T * inte_f_phi;  
end



%%%------利用电场关于kz的函数,利用高斯积分法计算V_ref------%%%
function V_ref = cal_V_ref()
%积分求解电势
  global rho_R kPickIntervalNumVref
  k_z_min = 10^(-10);
  k_z_max = 10^7;
  k_z = get_interval_split(k_z_min, k_z_max , kPickIntervalNumVref);
                                           %积分取点坐标
  
  inte_V = 0;                              %初始化积分因子的值
  for interval_num = 1:kPickIntervalNumVref
    inte_V = inte_V + ... 
             cal_inte_cosE_kz(k_z(interval_num), k_z(interval_num+1));
                                           %累加计算积分                                       
  end
  
  V_ref = 4 * pi * rho_R * inte_V;
end

%%%%%%%%

function inte_cosE_kz = cal_inte_cosE_kz(k_z_j, k_z_jplus1)
%使用高斯—勒让德计算每一个小区间的积分值
  global z_R t W;
  
  k_z = (k_z_jplus1 + k_z_j)/2 + ((k_z_jplus1 - k_z_j)/2) * t;  
                                           %计算零点t对应的k_z的值
  g_kz = array_cal_E_phi_kz(k_z) .* cos(k_z*z_R);
                                           %计算k_z值对应的被积函数值
  inte_cosE_kz = (k_z_jplus1 - k_z_j)/2 * sum(g_kz .* W);    
                                           %计算简化后的积分值
end
  
%%%%%%%%

function interval_split = get_interval_split(lower, upper, split_num)
%获取分区间坐标点的函数
  k_1 = lower;                             %区间的下界
  k_max = upper;                           %区间的上界
  delta_k = log(k_max/k_1)/split_num;      %区间的增长因子
  
  interval_split = zeros(1,split_num+1);   %初始化坐标点数据存储数组
  for interval_split_index=1:split_num+1 
    interval_split(interval_split_index) = ...
          k_1*exp((interval_split_index-1)*delta_k);
  end 
end

%%%-----------------计算phi方向电场值-----------------------%%%
function array_E_phi_kz = array_cal_E_phi_kz(array_kz)
%支持输入k_z数组计算E并将输出数据打包为数组的模块
  global kUseFMinusI
  
  array_size = length(array_kz);          %获取输入数组的长度, 5
  
  array_E_phi_kz = zeros(1, array_size);  %初始化数据存储数组
  
  for kz_index = 1:array_size
    array_E_phi_kz(kz_index) = cal_E_phi_kz(array_kz(kz_index), kUseFMinusI);
  end                                     %逐个计算E的数值
end

%%%%%%%%

function E_phi_kz = cal_E_phi_kz(k_z, use_F_minus_I)
%计算电场数值的函数
  global rho_T rho_R I_T omega mu m
  
  gamma_m = cal_gamma_n(m, k_z);
         
  c_m = 1i * I_T * rho_T / 4.0 * gamma_m * [0;1];%得到c_m的数据  
  
  if rho_R >= rho_T 
    Fm_phi_plus = get_F_phi(use_F_minus_I, +1, gamma_m, k_z);
    E_phi_kz = [0, 1i*omega*mu(m)/gamma_m] * (Fm_phi_plus * c_m *...
                cal_multiply_besseljh1(1, 1, gamma_m*rho_T, gamma_m*rho_R));
  else
    Fm_phi_minus = get_F_phi(use_F_minus_I, -1, gamma_m, k_z);
    E_phi_kz = [0, 1i*omega*mu(m)/gamma_m] * (Fm_phi_minus * c_m *...
                cal_multiply_besseljh1(1, 1, gamma_m*rho_R, gamma_m*rho_T));
  end                                             %判断使用哪个公式计算电场E 
end

%%%%%%%%

function F_phi_pm = get_F_phi(use_F_minus_I, sign, gamma_m, k_z)
%计算F矩阵的函数
  global m rho_layer rho_R rho_T
  
  %通过调用函数, 获得狭义透射反射系数数据
  %GAMMA_out中按顺序存储着GAMMA_12, GAMMA_23,... ,GAMMA_(N-1)N
  %GAMMA_in中按顺序存储着GAMMA21, GAMMA_32 ,... , GAMMA_N(N-1) 
  %T类似
  [GAMMA_out, T_out, GAMMA_in, T_in] = cal_local_reflect_factor(k_z);
  
  %得到所需的两个广义反射系数的值
  tilde_GAMMA_m_mplus1 = cal_general_reflect_factor...
                         (m, m+1, k_z, GAMMA_out, T_out, GAMMA_in, T_in);
  tilde_GAMMA_m_mminus1 = cal_general_reflect_factor...
                         (m, m-1, k_z, GAMMA_out, T_out, GAMMA_in, T_in); 
                                                  
  if sign == 1
    %%计算Fphi+的值  
    F_phi_pm = (eye(2) + tilde_GAMMA_m_mplus1 *...
                         cal_frac_besselh1(0, 1, gamma_m*rho_layer(m), ...
                                                 gamma_m*rho_R)...
                         *...
                         cal_frac_besselj(1, 0, gamma_m*rho_R, ...
                                                gamma_m*rho_layer(m))...
               )... 
               *...
               (...
                 (eye(2) - tilde_GAMMA_m_mminus1 * tilde_GAMMA_m_mplus1 *...
                           cal_frac_besselj(0, 0, gamma_m*rho_layer(m-1), ...
                                                  gamma_m*rho_layer(m))...   
                           *...
                           cal_frac_besselh1(0, 0, gamma_m*rho_layer(m), ...
                                                   gamma_m*rho_layer(m-1))...
                 )...
                 \...
                 (eye(2) + tilde_GAMMA_m_mminus1 * ...
                           cal_frac_besselj(0, 1, gamma_m*rho_layer(m-1), ...
                                                  gamma_m*rho_T)...  
                           *...
                           cal_frac_besselh1(1, 0, gamma_m*rho_T, ...
                                                   gamma_m*rho_layer(m-1))...
                 )...
               );
             
  elseif sign == -1
    %%计算Fphi-的值         
    F_phi_pm = (eye(2) + tilde_GAMMA_m_mminus1 *...
                         cal_frac_besselj(0, 1, gamma_m*rho_layer(m-1), ...
                                                gamma_m*rho_R)... 
                         *...
                         cal_frac_besselh1(1, 0, gamma_m*rho_R, ...
                                                 gamma_m*rho_layer(m-1))... 
                )... 
                *...
                (...
                  (eye(2) - tilde_GAMMA_m_mplus1 * tilde_GAMMA_m_mminus1 * ...
                            cal_frac_besselj(0, 0, gamma_m*rho_layer(m-1), ...
                                                   gamma_m*rho_layer(m))... 
                            *...
                            cal_frac_besselh1(0, 0, gamma_m*rho_layer(m), ...
                                                    gamma_m*rho_layer(m-1))...
                  )... 
                  \...
                  (eye(2) + tilde_GAMMA_m_mplus1 * ...
                            cal_frac_besselj(1, 0, gamma_m*rho_T, ...
                                                   gamma_m*rho_layer(m))... 
                            *...
                            cal_frac_besselh1(0, 1, gamma_m*rho_layer(m), ...
                                                    gamma_m*rho_T)...
                  )...
                );
  else
    disp('F_phi_pluse 和 F_phi_minus计算出现问题,符号选取变量只能是+1或-1.');
    return;
  end
  
  if use_F_minus_I
    F_phi_pm = F_phi_pm - eye(2);           %减去单位矩阵以求解相对电压大小
  end
end

%%%%%%%%

function general_reflect_factor = cal_general_reflect_factor...
         (layer_from, layer_to, k_z, GAMMA_out, T_out, GAMMA_in, T_in)
%计算广义反射系数的函数
  
  global rho_layer total_layer_num
  
  if layer_from>layer_to 
    general_reflect_factor = GAMMA_in{1};         %初始化递归数值
    
    for n = 2:layer_to
      gamma_n = cal_gamma_n(n, k_z); 
      
      mmminus1_frac_besselh1_multiply_frac_besselj = ...   
              cal_frac_besselh1(0, 0, gamma_n*rho_layer(n), ...
                                      gamma_n*rho_layer(n-1))...
              *...
              cal_frac_besselj(0, 0, gamma_n*rho_layer(n-1), ...
                                     gamma_n*rho_layer(n));
                                                  %避免重复运算, 优化计算速度
      
      %新一代透射系数
      general_transmit_factor = ...
          (eye(2) - GAMMA_out{n} * general_reflect_factor *...
                    mmminus1_frac_besselh1_multiply_frac_besselj)\T_in{n};                       
      
      %新一代反射系数
      general_reflect_factor = ...
          T_out{n} * general_reflect_factor * general_transmit_factor *...
           mmminus1_frac_besselh1_multiply_frac_besselj + GAMMA_in{n};                
    end
    
  else
    N = total_layer_num - 1;
    general_reflect_factor = GAMMA_out{N};
    
    for n = N-1:-1:layer_from
      gamma_nplus1 = cal_gamma_n(n+1, k_z);
      
       mmplus1_frac_besselh1_multiply_frac_besselj = ...
               cal_frac_besselj(0, 0, gamma_nplus1*rho_layer(n), ...
                                      gamma_nplus1*rho_layer(n+1))... 
               *...
               cal_frac_besselh1(0, 0, gamma_nplus1*rho_layer(n+1), ...
                                       gamma_nplus1*rho_layer(n));
      
      %新一代的透射系数
      general_transmit_factor = ...
          (eye(2) - GAMMA_in{n} * general_reflect_factor * ...
                    mmplus1_frac_besselh1_multiply_frac_besselj)\T_out{n};                    
      
      %新一代反射系数
      general_reflect_factor = ...
          T_in{n} * general_reflect_factor * general_transmit_factor *...
          mmplus1_frac_besselh1_multiply_frac_besselj + GAMMA_out{n};
    end
  end
end

%%%%%%%%
function [GAMMA_out, T_out, GAMMA_in, T_in] = cal_local_reflect_factor(k_z)
%求解局域反射透射系数的函数  
  global total_layer_num
  
  T_out = cell(1,total_layer_num-1);
  GAMMA_in = cell(1,total_layer_num-1);
  GAMMA_out = cell(1,total_layer_num-1);
  T_in = cell(1,total_layer_num-1);         %初始化狭义反射透射系数矩阵的cell数组
  
  for n = 1:(total_layer_num-1)
    [j_nn, j_nplus1n, h_nn, h_nplus1n] = cal_jh_matrix(n, k_z);
    
    T_out{n} = (j_nn-h_nplus1n)\(j_nn-h_nn);
    GAMMA_in{n} = (j_nn-h_nplus1n)\(j_nplus1n-j_nn);
    GAMMA_out{n} = (h_nplus1n-j_nn)\(h_nn-h_nplus1n);
    T_in{n} = (h_nplus1n-j_nn)\(h_nplus1n-j_nplus1n);    
  end 
end

%%%%%%%%
function [j_nn, j_nplus1n, h_nn, h_nplus1n] = cal_jh_matrix(n, k_z)
%计算j和h矩阵的函数
  global epsilon_star omega rho_layer mu 
  
  gamma_n = cal_gamma_n(n, k_z);
  gamma_nplus1 = cal_gamma_n(n+1, k_z);

  j_nn = [0,mu(n);-epsilon_star(n),0]...
          / gamma_n * 1i * omega * ...
          cal_frac_besselj(1, 0, gamma_n*rho_layer(n),...
                                 gamma_n*rho_layer(n));
  
  j_nplus1n = [0,mu(n+1);-epsilon_star(n+1),0]...
               / gamma_nplus1 * 1i * omega * ...
               cal_frac_besselj(1, 0, gamma_nplus1*rho_layer(n),...
                                      gamma_nplus1*rho_layer(n));
  
  h_nn = [0,mu(n);-epsilon_star(n),0] ...
          / gamma_n * 1i * omega * ...
          cal_frac_besselh1(1, 0, gamma_n*rho_layer(n), ...
                                  gamma_n*rho_layer(n));
      
  h_nplus1n = [0,mu(n+1);-epsilon_star(n+1),0]... 
               / gamma_nplus1 * 1i * omega * ...
               cal_frac_besselh1(1, 0, gamma_nplus1*rho_layer(n), ...
                                       gamma_nplus1*rho_layer(n));
end

%%%%%%%%

function gamma_n = cal_gamma_n(n, k_z)
%计算gamma_n的函数
  global k
  gamma_n = sqrt(k(n)^2 - k_z^2);
    if imag(gamma_n) < 0
      gamma_n = -gamma_n;
    end                                         %gamma_n只取虚部是正数的解
end

%%%%%%%%

function frac_besselj = cal_frac_besselj(up_n, down_n, up_x, down_x)
%计算第一类bessel函数的比值形式的值
%up_x, up_n是 分子 上第一类bessel函数的阶数和函数自变量值
%down_x, down_n是 分母 上第一类bessel函数的阶数和函数自变量值

  if up_n~=1 && down_n~=1 && up_n~=0 && down_n~=0
      disp('cal_frac_besselj不支持1,0之外的n的运算!');
      return;
  end
  
  %数据溢出关键是虚部, 此600界限由实际测试得.
  %此处取abs的原因可由, bessel函数在自变量的模趋于无穷时具体的近似形式
  %和e指数上的指数大于600时会导致的数据溢出,以及bessel函数同时在分子分母上
  %综合分析得到, 以下类似
  if abs(imag(up_x))>600 || abs(imag(down_x))>600                                                
    frac_besselj = sqrt(down_x/up_x)*exp(1i*(down_x-up_x));  
    if down_n>up_n
      frac_besselj = -1i * frac_besselj;    
    elseif down_n<up_n
      frac_besselj = 1i * frac_besselj;   
    end                                     %如果bessel函数自变量虚部abs大于600, 
                                            %则使用近似法处理
  else
    frac_besselj = besselj(up_n, up_x)/besselj(down_n, down_x);    
  end                                       %如果bessel函数自变量虚部abs小于600, 
                                            %则直接调用内部函数计算
end

function frac_besselh1 = cal_frac_besselh1(up_n, down_n, up_x, down_x)
%计算汉克尔函数的比值形式的值
  if up_n~=1 && down_n~=1 && up_n~=0 && down_n~=0
      disp('cal_frac_besselh不支持1,0之外的n,m的运算!');
      return;
  end
  
  %与cal_frac_besselj函数的套路相同
  if abs(imag(up_x))>600 || abs(imag(down_x))>600     
    frac_besselh1 = sqrt(down_x/up_x)*exp(1i*(up_x-down_x));  
    if down_n>up_n
      frac_besselh1 = 1i * frac_besselh1;    
    elseif down_n<up_n
      frac_besselh1 = -1i * frac_besselh1;
    end                                      
  else                                      
    frac_besselh1 = besselh(up_n,1,up_x)/besselh(down_n,1,down_x);    
  end                                      
end

function multiply_besseljh1 = cal_multiply_besseljh1(j_n, h1_n, j_x, h1_x)
%计算汉克尔函数乘以第一类贝塞尔函数的值
  if j_n~=1 && h1_n~=1 && j_n~=0 && h1_n~=0
      disp('cal_multiply_besseljh1不支持1,0之外的n,m的运算!');
      return;
  end
  
  %与cal_frac_besselj函数的套路相同
  if abs(imag(j_x))>600 || abs(imag(h1_x))>600   
    multiply_besseljh1 = (1/(pi*sqrt(j_x*h1_x)))*exp(1i*(h1_x-j_x));  
    if h1_n>j_n
      multiply_besseljh1 = -1i * multiply_besseljh1;    
    elseif h1_n<j_n
      multiply_besseljh1 = 1i * multiply_besseljh1;
    end                                    
  else                                      
    multiply_besseljh1 = besselh(h1_n,1,h1_x)*besselj(j_n,j_x);    
  end                                
end

%%%-----------------其他辅助函数-----------------------%%%
function device_layer = get_device_layer()
%计算装置位于哪一层介质中的函数
  global rho_layer rho_R rho_T total_layer_num
  
  for search_index = 1:(total_layer_num-1)
    if rho_layer(search_index) > rho_R
      device_layer_R = search_index;        %探测位于哪一个区间内
      break;
    else
      device_layer_R = total_layer_num;     %如果不在分界线数据范围内, 那么一定位于最外层
    end
  end
  
  for search_index = 1:(total_layer_num-1)
    if rho_layer(search_index) > rho_T
      device_layer_T = search_index;        %探测位于哪一个区间内
      break;
    else
      device_layer_T = total_layer_num;     %如果不在分界线数据范围内, 那么一定位于最外层
    end
  end
  
  if device_layer_R ~= device_layer_T
    disp('探测线圈半径与测量线圈区域不符, 无法计算!');
    disp('如要升级软件, 请修改get_device_layer()函数');
    return;                                   %终止整个程序
  end
  
  device_layer = device_layer_T;
end

%============================================================
%--------------------------程序日志---------------------------
%============================================================
%--2017.11.5
%  -1- 程序编写完成并进行第一次运行. 出现大量错误.
%  -2- 修正了程序出现的若干语法错误, 但仍有数据溢出, 矩阵奇异的警告
%
%--2017.11.6
%  -1- 使用bessel函数在无穷远处的近似, 修正了矩阵数据溢出的问题
%  -2- 数据量级出现严重问题
%
%--2017.11.6
%  -1- 发现对文献阅读理解有误, 没有使用所谓的[F-I]算法
%  -2- 发现每次都完全跑一遍程序十分费时, 因此创建了'数据测试模式'
%  -3- 实现[F-I]算法后, 成功算出理想结果
%
%--2017.11.8
%  -1- 对算法进行优化, 删除了冗余的运算步骤, 速度提高了2倍左右
%  -2- 添加了欢迎界面和自由复制声明(copyleft)
%  -3- 创建processing bar
%  -4- 创建计时系统
%  -5- 3层计算时间: 21.9822s; 5层计算时间: 43.3241s
%
%--2017.11.9
%  -1- 为了保证画图的平滑和美观, 将40个采样点增加为80个
%  -2- 为了保证计算数据稳定性, 将V_ref积分采样点由550个增加为1000个 
%  -3- 3层计算时间: 90.979s; 5层计算时间: 153.0233s 
%
%--2017.11.11
%  -1- 试图通过将rho_T, rho_R数值错开一个小量的手段, 使E_phi快速收敛, 
%      进而直接用F计算V_all
%  -2- 直接法不可用, 虽然E可通过过调节rho_R使其迅速收敛, 但E_phi测试图明确
%      告诉我们, 不使用[F-I]方法处理的E_phi的虚部远远大于实部(大约两个量级),
%      且一直大于实部. 因此,积分值不可避免的会出现较大差异. 然而这与实际是不符
%      的!(理论问题?编程数值处理问题?)
%