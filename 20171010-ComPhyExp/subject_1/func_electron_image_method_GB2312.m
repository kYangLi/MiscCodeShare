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
%%%!!!����������!!!
%%�������: ��ֵ���㾵�����õ糡
%%���ʱ��: 2017.10.21
%%Ԥ��ʱ��: 1h
%%Ŀ��: ��Ϥmatlab�Ļ����﷨��������·
%%����: �õ��񷨼��㵼����������һ��ɵĵ糡�Ƴ��ֲ�
%%%ע��:����matlab�汾����, �е�matlab�汾���ܲ���֧�ֽ��ű���
%%%     ����д��һ��, ֻ��Ҫ���ű��ͺ����ֿ�����, ���ǽ��ű����
%%%     ����.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function func_electron_image_method_GB2312()
  %#########################################
  %%%%%%%%%%%%����ǰ��׼��%%%%%%%%%%%%%%%%%%%%
  %#########################################
  clear;   
  format long; %���ü��㾫��
 
  %====================================================
  %%%%%%%%%%%%%%%%%%������������Ķ���%%%%%%%%%%%%%%%%%%%%%
  %====================================================
  global q epsi_0 d a kPickThetaNum kTotalThetaSet; 
  q = 1;                              %�����, C
  epsi_0 = 8.854187818*10^(-12);      %��ս�糣��
  d = 1;                              %���ɵ����ĵľ���,m
  a = 0.5;                            %��ǵİ뾶,m
  kPickThetaNum = 200;                %ȡthetaֵ�ĸ���
  kTotalThetaSet = ...
  linspace(0, 2*pi, kPickThetaNum);   %��theta��0��2*pi����ȡ��

  %====================================================
  %%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %====================================================
  %---------����Ĭ�ϲ���------------
  r_min = 0.5;                        %Ĭ��rȡֵ����Сֵ,m
  r_step = 0.01;                      %Ĭ��rȡֵ�Ĳ���,m
  r_max = 1.5;                        %Ĭ��rȡֵ�����ֵ,m
  plot_sleep_time = 0.05;             %Ĭ�϶�ͼ 1/֡�� (��ͼ����ʱ��),s

  %---------��ȡ�û�����------------
  use_default_parameter = input('�Ƿ��޸�Ĭ�ϲ���(y/n)\n>','s');

  if strcmpi('Y', use_default_parameter) || strcmpi('y', use_default_parameter)
    r_min = input(['�����������r������(m, ���������',num2str(a),'��ֵ)\n>']);
    r_max = input('�����������r������(m)\n>');
    r_step = input('�����������r�Ĳ���(m)\n>');
    plot_sleep_time = input('������һ��ͼ��ͣ��ʱ��(s)\n>');
  end

  %---------���㲢��ȡ����-----------
  %��ȡ����ϵĵ���ܶ�����
  sigma_e = get_shell_electric_density();
  %��ȡ����糡����������
  [set_U, set_E_r, set_E_theta, set_E_phi] = ...
      pack_up_electric_field_data(r_min,r_step,r_max);

  %--------���ݵĻ��������-----------
  %����������ݵ�ָ���ļ�
  save_data_to_file(sigma_e, ...
                    set_U, set_E_r, set_E_theta, set_E_phi,...
                    r_min, r_max, r_step);
               
  %ʹ�����ݻ���ͼ��
  plot_data(sigma_e, ...
            set_U, set_E_r, set_E_theta, set_E_phi,...
            r_min, r_max, r_step, plot_sleep_time);

end
%====================================================
%%%%%%%%%%%%%%%%%%���ܺ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%====================================================
function save_data_to_file(sigma_e, ...
                  set_U, set_E_r, set_E_theta, set_E_phi,...
                  r_min, r_max, r_step)
%���ڴ洢���ݵ��ļ��ĺ���
  global kTotalThetaSet;  

  data_file_handle = fopen('func_electron_image_method_data.txt', 'w');
  %�����ĵ���˵��
  fprintf(data_file_handle, '############˵����ʼ##############\n');
  fprintf(data_file_handle, '#���ĵ���func_electron_image_method.m���ɵ������ĵ�\n');
  fprintf(data_file_handle, '#�����Ը����Լ���������ȡ����\n');
  fprintf(data_file_handle, '#��ͬ���ݵı�ʶ��Ϊ <==:datakind:==> \n');
  fprintf(data_file_handle, '#��ͬ��datakind����ͬ������\n');
  fprintf(data_file_handle, '#��ά����ĺ����ǲ�ͬ��theta\n');
  fprintf(data_file_handle, '#��ά����������ǲ�ͬ��r\n');
  fprintf(data_file_handle, '#theta��r��ȡֵ�������ļ��и���\n');
  fprintf(data_file_handle, '#��ȡ���ļ�����������Ҫ���б�д��س���\n');
  fprintf(data_file_handle, '####copyleft---liyang---####\n');
  fprintf(data_file_handle, '############˵������##############\n\n\n\n');
  fprintf(data_file_handle, '>>>DATA_READ_POINT<<<\n');
  
  %���theta������
  fprintf(data_file_handle, '<==:theta:==>\n');
  fprintf(data_file_handle, '%d  ', kTotalThetaSet);
  
  %���r������
  fprintf(data_file_handle, '\n\n<==:r:==>\n');
  fprintf(data_file_handle, '%d\n', r_min:r_step:r_max);
  
  %���sigma_e������
  fprintf(data_file_handle, '\n\n<==:sigma_e:==>\n');
  fprintf(data_file_handle, '%d  ', sigma_e);
  
  %���U������
  fprintf(data_file_handle, '\n\n<==:U:==>\n');
  fprintf_2_d_array(data_file_handle, set_U);
  
  %���E_r������
  fprintf(data_file_handle, '\n\n<==:E_r:==>\n');
  fprintf_2_d_array(data_file_handle, set_E_r);
  
  %���E_theta������
  fprintf(data_file_handle, '\n\n<==:E_theta:==>\n');
  fprintf_2_d_array(data_file_handle, set_E_theta);
   
  %���set_E_phi������
  fprintf(data_file_handle, '\n\n<==:E_phi:==>\n');
  fprintf_2_d_array(data_file_handle, set_E_phi);
  
  fclose(data_file_handle);
end

%%%%%%%%%

function fprintf_2_d_array(data_file_handle, two_d_array)
%�������ļ��д洢2ά�������ݵĺ���
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
%���ڻ�������ͼ��ĺ���
  global kTotalThetaSet;  

  %���Ƶ���ܶȵ�ͼ��
  subplot(3,2,1);
  plot(kTotalThetaSet, sigma_e,...
       'LineWidth',2,...
       'color','red');
  axis([0,2*pi,-0.4,0]);
  set(gca,'XTick',0:pi/2:2*pi,'xtickLabel',{'0','��/2','��','3��/2','2��'});
  title('��_e(��) picture');

  scrsz = get(0,'ScreenSize'); %��ȡ��Ļ�ֱ���
  set(gcf,'position', scrsz);  %���û�ͼ����ȫ��
  
  %���Ƶ糡��ص�ͼ��
  r = r_min:r_step:r_max;
  for r_index = 1:length(r)
    %����U����theta��ͼ��
    subplot(3,2,3);
    plot(kTotalThetaSet, set_U(r_index,:),...
         'LineWidth',2,...
         'color','red');
    axis([0,2*pi,0,6*10^11]);
    set(gca,'XTick',0:pi/2:2*pi,'xtickLabel',{'0','��/2','��','3��/2','2��'});
    title(['U(��) picture, r = ',num2str(r(r_index)),'m']);

    %����E_r����theta��ͼ��
    subplot(3,2,4);
    plot(kTotalThetaSet, set_E_r(r_index,:),...
         'LineWidth',2,...
         'color','red');
    axis([0,2*pi,-5*10^15,0]);
    set(gca,'XTick',0:pi/2:2*pi,'xtickLabel',{'0','��/2','��','3��/2','2��'});
    title(['E_r(��) picture, r = ',num2str(r(r_index)),'m']);

    %����E_theta����theta��ͼ��
    subplot(3,2,5);
    plot(kTotalThetaSet, set_E_theta(r_index,:),...
         'LineWidth',2,...
         'color','red');
    axis([0,2*pi,-5*10^12,5*10^12]);
    set(gca,'XTick',0:pi/2:2*pi,'xtickLabel',{'0','��/2','��','3��/2','2��'});
    title(['E_��(��) picture, r = ',num2str(r(r_index)),'m']);

    %����E_phi����theta��ͼ��
    subplot(3,2,6);
    plot(kTotalThetaSet, set_E_phi(r_index,:),...
         'LineWidth',2,...
         'color','red');
    axis([0,2*pi,-0.1,0.1]);
    set(gca,'XTick',0:pi/2:2*pi,'xtickLabel',{'0','��/2','��','3��/2','2��'});
    title(['E_��(��) picture, r = ',num2str(r(r_index)),'m']);

    %һ֡ͼƬ��ͣ��ʱ������
    pause(plot_sleep_time);
  end
end

%%%%%%%%

function [set_U, set_E_r, set_E_theta, set_E_phi] = ...
         pack_up_electric_field_data(r_min,r_step,r_max)
%�糡�������ݴ���������(���ݷ�װ����)
  
  %�����������ݰ��ĳ�ֵ
  set_U = [];
  set_E_r = [];
  set_E_theta = [];
  set_E_phi = [];
  
  %������������
  for r = r_min:r_step:r_max
    [U, E_r, E_theta, E_phi] = get_electric_field_characteristic(r);
    
    %�ö�̬��������ķ���, �洢����
    %>>�˷���������������ٶȵĽ���<<
    set_U = [set_U; U];
    set_E_r = [set_E_r; E_r];
    set_E_theta= [set_E_theta; E_theta];
    set_E_phi = [set_E_phi; E_phi];
  end
end

%%%%%%%%

function [U, E_r, E_theta, E_phi] = get_electric_field_characteristic(r)
%��������뾶�ϵ�������theta�ķֲ�
%��������뾶�ϵĵ糡����������theta�ķֲ�
  
  global q epsi_0 d a kPickThetaNum kTotalThetaSet
  
  %�������ݴ洢�����С(����)
  %�Է�ֹ����ö�̬����
  %��̬���齫ÿһ�ֶ���������ĳ���
  %������ö�̬������������������ٶ�.
  U = zeros(1,kPickThetaNum);       
  E_r = zeros(1,kPickThetaNum);
  E_theta = zeros(1,kPickThetaNum);
  E_phi = zeros(1,kPickThetaNum);   
  
  for theta_step = 1:kPickThetaNum
    %ÿ��ȡһ��theta��������
    theta = kTotalThetaSet(theta_step); 
    
    %����ȷ���뾶�ϵĲ�ͬtheta���Ƴ�ֵ
    U(theta_step) = q/(4*pi*epsi_0) * ...
                      (1/sqrt(r^2 + d^2 - 2*r*d*cos(theta)) - ...
                       1/sqrt((d^2/a^2)*r^2 + a^2 - 2*r*d*cos(theta)));
    %����ȷ���뾶�ϵĲ�ͬtheta�ĵ糡r����ֵ
    E_r(theta_step) = q/(4*pi*epsi_0) * ...
                     (r-d*cos(theta)/(r^2 + d^2 - 2*r*d*cos(theta))^(1.5) - ...
                     (d^2/a^2)*r-d*cos(theta)/((d^2/a^2)*r^2 + a^2 - 2*r*d*cos(theta))^(1.5));
    %����ȷ���뾶�ϵĲ�ͬtheta�ĵ糡theta����ֵ
    E_theta(theta_step) = q*d*sin(theta)/(4*pi*epsi_0) * ...
                          (1/(r^2 + d^2 - 2*r*d*cos(theta))^(1.5) - ...
                          1/((d^2/a^2)*r^2 + a^2 - 2*r*d*cos(theta))^(1.5));
    %�糡phi����ֵ��Ϊ0
    E_phi(theta_step) = 0;
  end 
end

%%%%%%%%

function electric_density = get_shell_electric_density()
%��������ϵĵ���ܶȵķֲ�  
  global q  d a kPickThetaNum kTotalThetaSet
  
  %�������ݴ洢�����С(����)
  electric_density = zeros(1,kPickThetaNum);
  
  for theta_step = 1:kPickThetaNum
    theta = kTotalThetaSet(theta_step);
   
    %����r=aʱ�ĵ���ܶ�
    electric_density(theta_step) = -q/(4*pi) * (d^2 - a^2)/...
                        a * (a^2 + d^2 - -2*a*d*cos(theta))^(1.5);
  end
end
