%%%######################%%%
%%%#Copyleft Liyang 2017#%%%
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
%%程序标题: 数值计算镜像法所得电场
%%编程时间: 2017.10.21
%%预计时长: 1h
%%目的: 熟悉matlab的基本语法和运算套路
%%功能: 用电像法计算导电球壳外放置一电荷的电场势场分布
%%%注意:由于matlab版本问题, 有的matlab版本可能并不支持将脚本和
%%%     函数写在一起, 只需要将脚本和函数分开或存放, 或是将脚本编成
%%%     即可.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function func_electron_image_method()
  %#########################################
  %%%%%%%%%%%%计算前的准备%%%%%%%%%%%%%%%%%%%%
  %#########################################
  clear;   
  format long; %设置计算精度
 
  %====================================================
  %%%%%%%%%%%%%%%%%%基本物理参数的定义%%%%%%%%%%%%%%%%%%%%%
  %====================================================
  global q epsi_0 d a kPickThetaNum kTotalThetaSet; 
  q = 1;                              %电荷量, C
  epsi_0 = 8.854187818*10^(-12);      %真空介电常数
  d = 1;                              %点电荷到球心的距离,m
  a = 0.5;                            %球壳的半径,m
  kPickThetaNum = 200;                %取theta值的个数
  kTotalThetaSet = ...
  linspace(0, 2*pi, kPickThetaNum);   %将theta从0到2*pi均匀取点

  %====================================================
  %%%%%%%%%%%%%%%%%%主程序%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %====================================================
  %---------定义默认参数------------
  r_min = 0.5;                        %默认r取值的最小值,m
  r_step = 0.01;                      %默认r取值的步长,m
  r_max = 1.5;                        %默认r取值的最大值,m
  plot_sleep_time = 0.05;             %默认动图 1/帧数 (画图休眠时间),s

  %---------获取用户输入------------
  use_default_parameter = input('是否修改默认参数(y/n)\n>','s');

  if strcmpi('Y', use_default_parameter) || strcmpi('y', use_default_parameter)
    r_min = input(['请输入计算中r的下限(m, 请输入大于',num2str(a),'的值)\n>']);
    r_max = input('请输入计算中r的上限(m)\n>');
    r_step = input('请输入计算中r的步长(m)\n>');
    plot_sleep_time = input('请输入一张图的停留时间(s)\n>');
  end

  %---------计算并获取数据-----------
  %获取球壳上的电荷密度数据
  sigma_e = get_shell_electric_density();
  %获取球外电场的性质数据
  [set_U, set_E_r, set_E_theta, set_E_phi] = ...
      pack_up_electric_field_data(r_min,r_step,r_max);

  %--------数据的绘制与输出-----------
  %输出所得数据到指定文件
  save_data_to_file(sigma_e, ...
                    set_U, set_E_r, set_E_theta, set_E_phi,...
                    r_min, r_max, r_step);
               
  %使用数据绘制图像
  plot_data(sigma_e, ...
            set_U, set_E_r, set_E_theta, set_E_phi,...
            r_min, r_max, r_step, plot_sleep_time);

end
%====================================================
%%%%%%%%%%%%%%%%%%功能函数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%====================================================
function save_data_to_file(sigma_e, ...
                  set_U, set_E_r, set_E_theta, set_E_phi,...
                  r_min, r_max, r_step)
%用于存储数据到文件的函数
  global kTotalThetaSet;  

  data_file_handle = fopen('func_electron_image_method_data.txt', 'w');
  %数据文档的说明
  fprintf(data_file_handle, '############说明开始##############\n');
  fprintf(data_file_handle, '#本文档是func_electron_image_method.m生成的数据文档\n');
  fprintf(data_file_handle, '#您可以根据自己的需求提取数据\n');
  fprintf(data_file_handle, '#不同数据的标识符为 <==:datakind:==> \n');
  fprintf(data_file_handle, '#不同的datakind代表不同物理量\n');
  fprintf(data_file_handle, '#二维数组的横向是不同的theta\n');
  fprintf(data_file_handle, '#二维数组的纵向是不同的r\n');
  fprintf(data_file_handle, '#theta和r的取值均已在文件中给出\n');
  fprintf(data_file_handle, '#提取本文件的数据您需要自行编写相关程序\n');
  fprintf(data_file_handle, '####copyleft---liyang---####\n');
  fprintf(data_file_handle, '############说明结束##############\n\n\n\n');
  fprintf(data_file_handle, '>>>DATA_READ_POINT<<<\n');
  
  %输出theta的数据
  fprintf(data_file_handle, '<==:theta:==>\n');
  fprintf(data_file_handle, '%d  ', kTotalThetaSet);
  
  %输出r的数据
  fprintf(data_file_handle, '\n\n<==:r:==>\n');
  fprintf(data_file_handle, '%d\n', r_min:r_step:r_max);
  
  %输出sigma_e的数据
  fprintf(data_file_handle, '\n\n<==:sigma_e:==>\n');
  fprintf(data_file_handle, '%d  ', sigma_e);
  
  %输出U的数据
  fprintf(data_file_handle, '\n\n<==:U:==>\n');
  fprintf_2_d_array(data_file_handle, set_U);
  
  %输出E_r的数据
  fprintf(data_file_handle, '\n\n<==:E_r:==>\n');
  fprintf_2_d_array(data_file_handle, set_E_r);
  
  %输出E_theta的数据
  fprintf(data_file_handle, '\n\n<==:E_theta:==>\n');
  fprintf_2_d_array(data_file_handle, set_E_theta);
   
  %输出set_E_phi的数据
  fprintf(data_file_handle, '\n\n<==:E_phi:==>\n');
  fprintf_2_d_array(data_file_handle, set_E_phi);
  
  fclose(data_file_handle);
end

%%%%%%%%%

function fprintf_2_d_array(data_file_handle, two_d_array)
%用于向文件中存储2维数组数据的函数
  size_two_d_array = size(two_d_array);
  for row = 1:size_two_d_array(1)
    fprintf(data_file_handle, '%d  ', two_d_array(row,:)); 
    fprintf(data_file_handle, '\n');
  end
end

%%%%%%%%%

function plot_data(sigma_e, ...
                   set_U, set_E_r, set_E_theta, set_E_phi,...
                   r_min, r_max, r_step, plot_sleep_time)
%用于绘制数据图像的函数
  global kTotalThetaSet;  

  %绘制电荷密度的图像
  subplot(3,2,1);
  plot(kTotalThetaSet, sigma_e,...
       'LineWidth',2,...
       'color','red');
  axis([0,2*pi,-0.4,0]);
  set(gca,'XTick',0:pi/2:2*pi,'xtickLabel',{'0','π/2','π','3π/2','2π'});
  title('σ_e(θ) picture');

  scrsz = get(0,'ScreenSize'); %获取屏幕分辨率
  set(gcf,'position', scrsz);  %设置绘图窗体全屏
  
  %绘制电场相关的图像
  r = r_min:r_step:r_max;
  for r_index = 1:length(r)
    %绘制U关于theta的图像
    subplot(3,2,3);
    plot(kTotalThetaSet, set_U(r_index,:),...
         'LineWidth',2,...
         'color','red');
    axis([0,2*pi,0,6*10^11]);
    set(gca,'XTick',0:pi/2:2*pi,'xtickLabel',{'0','π/2','π','3π/2','2π'});
    title(['U(θ) picture, r = ',num2str(r(r_index)),'m']);

    %绘制E_r关于theta的图像
    subplot(3,2,4);
    plot(kTotalThetaSet, set_E_r(r_index,:),...
         'LineWidth',2,...
         'color','red');
    axis([0,2*pi,-5*10^15,0]);
    set(gca,'XTick',0:pi/2:2*pi,'xtickLabel',{'0','π/2','π','3π/2','2π'});
    title(['E_r(θ) picture, r = ',num2str(r(r_index)),'m']);

    %绘制E_theta关于theta的图像
    subplot(3,2,5);
    plot(kTotalThetaSet, set_E_theta(r_index,:),...
         'LineWidth',2,...
         'color','red');
    axis([0,2*pi,-5*10^12,5*10^12]);
    set(gca,'XTick',0:pi/2:2*pi,'xtickLabel',{'0','π/2','π','3π/2','2π'});
    title(['E_θ(θ) picture, r = ',num2str(r(r_index)),'m']);

    %绘制E_phi关于theta的图像
    subplot(3,2,6);
    plot(kTotalThetaSet, set_E_phi(r_index,:),...
         'LineWidth',2,...
         'color','red');
    axis([0,2*pi,-0.1,0.1]);
    set(gca,'XTick',0:pi/2:2*pi,'xtickLabel',{'0','π/2','π','3π/2','2π'});
    title(['E_φ(θ) picture, r = ',num2str(r(r_index)),'m']);

    %一帧图片暂停的时间设置
    pause(plot_sleep_time);
  end
end

%%%%%%%%

function [set_U, set_E_r, set_E_theta, set_E_phi] = ...
         pack_up_electric_field_data(r_min,r_step,r_max)
%电场分量数据打包输出函数(数据封装函数)
  
  %给定各个数据包的初值
  set_U = [];
  set_E_r = [];
  set_E_theta = [];
  set_E_phi = [];
  
  %计算数据内容
  for r = r_min:r_step:r_max
    [U, E_r, E_theta, E_phi] = get_electric_field_characteristic(r);
    
    %用动态增长数组的方法, 存储数据
    %>>此方法会引起计算上速度的降低<<
    set_U = [set_U; U];
    set_E_r = [set_E_r; E_r];
    set_E_theta= [set_E_theta; E_theta];
    set_E_phi = [set_E_phi; E_phi];
  end
end

%%%%%%%%

function [U, E_r, E_theta, E_phi] = get_electric_field_characteristic(r)
%计算给定半径上的势能随theta的分布
%计算给定半径上的电场各分量能随theta的分布
  
  global q epsi_0 d a kPickThetaNum kTotalThetaSet
  
  %定义数据存储数组大小(长度)
  %以防止其采用动态数组
  %动态数组将每一轮都更新自身的长度
  %因而采用动态数组会严重拖慢运算速度.
  U = zeros(1,kPickThetaNum);       
  E_r = zeros(1,kPickThetaNum);
  E_theta = zeros(1,kPickThetaNum);
  E_phi = zeros(1,kPickThetaNum);   
  
  for theta_step = 1:kPickThetaNum
    %每次取一个theta点来计算
    theta = kTotalThetaSet(theta_step); 
    
    %计算确定半径上的不同theta的势场值
    U(theta_step) = q/(4*pi*epsi_0) * ...
                      (1/sqrt(r^2 + d^2 - 2*r*d*cos(theta)) - ...
                       1/sqrt((d^2/a^2)*r^2 + a^2 - 2*r*d*cos(theta)));
    %计算确定半径上的不同theta的电场r分量值
    E_r(theta_step) = q/(4*pi*epsi_0) * ...
                     (r-d*cos(theta)/(r^2 + d^2 - 2*r*d*cos(theta))^(1.5) - ...
                     (d^2/a^2)*r-d*cos(theta)/((d^2/a^2)*r^2 + a^2 - 2*r*d*cos(theta))^(1.5));
    %计算确定半径上的不同theta的电场theta分量值
    E_theta(theta_step) = q*d*sin(theta)/(4*pi*epsi_0) * ...
                          (1/(r^2 + d^2 - 2*r*d*cos(theta))^(1.5) - ...
                          1/((d^2/a^2)*r^2 + a^2 - 2*r*d*cos(theta))^(1.5));
    %电场phi分量值恒为0
    E_phi(theta_step) = 0;
  end 
end

%%%%%%%%

function electric_density = get_shell_electric_density()
%计算球壳上的电荷密度的分布  
  global q  d a kPickThetaNum kTotalThetaSet
  
  %定义数据存储数组大小(长度)
  electric_density = zeros(1,kPickThetaNum);
  
  for theta_step = 1:kPickThetaNum
    theta = kTotalThetaSet(theta_step);
   
    %计算r=a时的电荷密度
    electric_density(theta_step) = -q/(4*pi) * (d^2 - a^2)/...
                        a * (a^2 + d^2 - -2*a*d*cos(theta))^(1.5);
  end
end
