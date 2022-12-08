%%%######################%%%
%%%#Copyleft 2017 Liyang#%%%
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
%%%��������: �⾮�ź�ģ��ϵͳ
%%%���ʱ��: 2017.11.4
%%%Ԥ����ʱ: 6h
%%%Ŀ��:�⾮����
%%%����:����matlab����⾮�ź�ǿ��
%%%ע��:����matlab�汾����, �е�matlab�汾���ܲ���֧�ֽ��ű���
%%%    ����д��һ��, ֻ��Ҫ���ű��ͺ����ֿ�����, ���ǽ��ű����
%%%    ��������.
%%%    ������ѡȡ�˵ڶ��ֽ����ʽ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function func_loggin_multilayer_GB2312()
  %====================================================
  %%%%%%%%%%%%%%%%%%����ǰ��׼��%%%%%%%%%%%%%%%%%%%%%%%%%
  %====================================================
  clear; clc;  
  format long;                              %���ü��㾫��
  
  disp('########################');
  disp('******�⾮���ۼ����ʽ*****');
  disp('########################');
  disp('Version 1.3.1');
  disp('last modified 2017-10-11')
  disp('Copyleft 2017 Liyang');
  disp(' ');
  %====================================================
  %%%%%%%%%%%%%%%%%%����ģʽ��ѡ��%%%%%%%%%%%%%%%%%%%%%%%%
  %====================================================
  %ѡ�����ģʽ����, 1:����, 0:�ر�
  global kOpenTestMode
  choice_cal_method = input('���Ƿ�Ҫ�������ģʽ(Y/N):', 's');
  if strcmpi('Y', choice_cal_method) || strcmpi('y', choice_cal_method)
    kOpenTestMode = 1;
    disp('��ѡ�������>>����ģʽ<<!');
  else
    kOpenTestMode = 0; 
    disp('��ѡ�������>>����ģʽ<<!');
  end
  disp(' ');
  
  %����ģʽ(���ýضϷ���ֱ����ⷨ��V)
  global kUseFMinusI
  choice_cal_method = input('���Ƿ�Ҫʹ��[F-I]��������(Y/N):', 's');
  if strcmpi('N', choice_cal_method) || strcmpi('n', choice_cal_method)
    kUseFMinusI = 0;
    disp('��ѡ�������>>[F]����ģʽ<<!');
  else
    kUseFMinusI = 1; 
    disp('��ѡ�������>>[F-I]����ģʽ<<!');
  end
  disp(' ');
  
  %ѡ��⾮ģ�͵Ĳ�������
  kTotalLayerNum = 3;                       
  data_list_char = ['1. һ��ģ��';'3. ����ģ��'; '5. ���ģ��'];
  data_list_number = [1, 3, 5];
  
  fprintf('�����Ե��õ����ݱ�ǩ��:\n');
  sizeof_data_list_char = size(data_list_char);
  for print_step = 1:sizeof_data_list_char(1)
    fprintf('-->');
    fprintf('%c', data_list_char(print_step,:));
    fprintf('\n');
  end
  fprintf('\n');
  while 1
    kTotalLayerNum = input('������Ҫ�����ģ�͵����ݱ�ǩ��Ӧ������:');
    
    if ismember(kTotalLayerNum, data_list_number)
      break;
    end
    
    disp('!!!��������ݲ����б���, ������ѡ��!!!');
  end
  
  tic;                                        %�����ʱ��, �����ʱ��ʼ
  disp('Timer Begin...');
  
  %====================================================
  %%%%%%%%%%%%%%%%%%������������Ķ���%%%%%%%%%%%%%%%%%%%%
  %====================================================
  %%%-----------�����Ĺ��ò���------------%%%
  global omega I_T z_R rho_T rho_R t W 
  global kPickIntervalNumVref kPickIntervalNumVdir kPickIntervalNumEphi
  f = 2000000;                                %�ź�Ƶ��,Hz
  omega = 2*pi*f;                             %�źŽ�Ƶ��
  I_T = 1;                                    %������Ȧ������С,A (T_T)
  z_R = 0.762;                                %����-̽����Ȧ����,m
  
  %%%ʹrho_T��rho_R����΢��ͬ, ��ʹE_phi��������%%%
  %%%����������ֱ�ӷ�(������[F-I]�滻��)��V%%%
  rho_T = 0.1143;                             %������Ȧ�İ뾶,m
  rho_R = rho_T+0.0001;                       %̽����Ȧ�İ뾶,m

  kPickIntervalNumVref = 1000;                %V_ref����ȡ�����
  kPickIntervalNumVdir = 300;                 %V_dir����ȡ�����
  kPickIntervalNumEphi = 1000;                %����Ephi��ͼ��ȡ����

  t = [-0.9061798459, -0.5384693101,...
        0, ...
        0.5384693101, 0.9061798459];          %L5��˹�����е����
  W = [0.2369268851, 0.4786286705,...
       0.568888889, ...
       0.4786286705, 0.2369268851];           %L5��˹�����еı���
  %%%+++++++++�����Ĺ��ò�������++++++++++%%%
  
  %%%--------�⾮��������ݿ�-------------%%%
  global mu rho_layer sigma epsilon_star k
  
  mu_0 = 1.25663706144 * 10^(-6);             %��ս��ʴŵ���
  epsilon_0 =  8.854187817*10^(-12);          %��ս�糣��
  
  %------5�����ݵĵ��ѧ���ʲ���---------
  if kTotalLayerNum == 5     
    mu_r = [1, 1, 1, 1, 1];                   %���ﵽ��������Դŵ���
    rho_layer = [0.1016, 0.1118, 0.1270, 0.3048];
                                              %����ķֽ�뾶ֵ 
    epsilon_r = [1, 1, 80, 1, 1];             %���ﵽ��������Խ�糣��
    sigma = [10^8, 10^(-5), 10, 1, 1];        %���ﵽ�����ĵ絼��, S/m
  
  %------3�����ݵĵ��ѧ���ʲ���---------
  elseif kTotalLayerNum == 3
    mu_r = [1, 1, 1];                         %���ﵽ��������Դŵ���
    rho_layer = [0.1016, 0.1270];             %����ķֽ�뾶ֵ
    epsilon_r = [1, 80, 1];                   %���ﵽ��������Խ�糣��
    sigma = [10^8, 10, 1];                    %���ﵽ�����ĵ絼��, S/m
  %------���Ƚ��ʵĵ�ѧ���ʲ���----------
  elseif kTotalLayerNum == 1
    mu_r = [1, 1, 1];                         %���ﵽ��������Դŵ���
    rho_layer = [0.1016, 0.1270];             %����ķֽ�뾶ֵ
    epsilon_r = [1, 1, 1];                    %���ﵽ��������Խ�糣��
    sigma = [1, 1, 1];                        %���ﵽ�����ĵ絼��, S/m
  %----------------------------------
  %----------------------------------
  end
  %+++++++++++����ѡ�����++++++++++++
  
  mu = mu_0 * mu_r;                           %���ﵽ�����ľ��Դŵ���
  epsilon = epsilon_0 * epsilon_r;            %���ﵽ�����ľ��Խ�糣��
  epsilon_star = epsilon+1i*(sigma/omega);    %���ﵽ�����ĸ���糣��
  k = sqrt(omega^2*mu.*epsilon_star);         %��Ų���ʸ��ֵ
  %%%++++++++�⾮��������ݿ����+++++++++++++%%%
  
  %%%------------������������--------------%%%
  global m total_layer_num
  total_layer_num = length(mu_r);             %��ȡ�ܹ��Ĳ���
  m = get_device_layer();                     %��ȡ̽��������һ����ʵ���
  %%%++++++++++����������������+++++++++++++%%%
  
  %====================================================
  %%%%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %====================================================
 
  if kOpenTestMode
    %%%�������ģ��
    clc;
    disp('================================');
    disp('!!!��Ŀǰ�������ݼ�����ȷ�Բ���ģʽ!!!');
    disp('================================');
    disp('���ڼ����ʼָ����sigma�µĵ�ѹֵ...');
    
    if kUseFMinusI                              %���ʹ����[F-I]���E_phi, ����Ҫ��V_dir
      test_V_all = cal_V_ref() + cal_V_dir(); 
    else                                        %���ʹ��[F]���E_phi����Ҫ��V_dir
      test_V_all = cal_V_ref();
    end
    disp('V_special = ');
    disp(test_V_all);
    disp('(V)');
    disp(' ');
    
    %------------------------����E_phi��ģ��--------------------
    disp('��������E_phi��ͼ��...');
    
    %��ʼ������k_z�����E_phi����
    k_z_test = get_interval_split(10^(-10), 10^15, kPickIntervalNumEphi-1);
    E_ref_array = zeros(1,kPickIntervalNumEphi);
    E_all_array = zeros(1,kPickIntervalNumEphi);
    
    %����E_phi����
    for interval_num = 1:kPickIntervalNumEphi
      %����(F-I)��Ӧ�ĵ糡ֵ
      E_ref_array(interval_num) = cal_E_phi_kz(k_z_test(interval_num), 1);
      %����(F)��Ӧ�ĵ糡ֵ
      E_all_array(interval_num) = cal_E_phi_kz(k_z_test(interval_num), 0);
    end
    
    %����ͼ��
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
                                              %���ô��ڴ�С
    legend([ref_imag_pic,ref_real_pic,all_imag_pic,all_real_pic],...
            '|Ref Image|',' |Ref Real|', '|Total Image|','|Total Real|');
                                              %����ͼ��
    xlabel('k_z');                     
    ylabel('E_{phi}(k)');
    title('Test for E_{phi}');
    saveas(gcf,['loggin_', num2str(kTotalLayerNum), '_layers_Ephi.jpg'])
                                              %�洢ͼ���ļ�
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
     disp('E_phi��ͼ���������!');
     disp(' ');
     disp('���ݲ��Խ���!'); 
     
  else 
    %%%������ʽ����ģ��
    %--------------��������ʾģ��.����(��).-----------------------------
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
    
    %---------------sigmaĩβԪ����ֵ�ĸı�------------------------------
    kPickSigmaFNum = 79;                          %�����sigma�仯ȡ�����(+1)
    kSigmaFLowerLimit = 10^(-3);                  %�����sigmaȡֵ����
    kSigmaFUpperLimit = 10;                       %�����sigmaȡֵ����

    sigma_end = get_interval_split...
               (kSigmaFLowerLimit, kSigmaFUpperLimit , kPickSigmaFNum);         
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    test_V_all = zeros(1,kPickSigmaFNum+1);                %��ʼ��V�洢����
    
    %%ѭ��������ɸ���ͼ��
    for step_sigma = 1:kPickSigmaFNum+1  
      %���¸�������
      sigma(end) = sigma_end(step_sigma);         %���ﵽ�����ĵ絼��, S/m
      epsilon_star = epsilon+1i*(sigma/omega);    %���ﵽ�����ĸ���糣��
      k = sqrt(omega^2*mu.*epsilon_star);         %��Ų���ʸ��ֵ
      
      %###############���ļ������#######################################
      if kUseFMinusI                              %���ʹ����[F-I]���E_phi, ����Ҫ��V_dir
        V_all = cal_V_ref() + cal_V_dir();
      else                                        %���ʹ��[F]���E_phi����Ҫ��V_dir
        V_all = cal_V_ref();
      end
      test_V_all(step_sigma) = V_all;             %���������Ȧ�ĵ���
      %################################################################
      
      %--------------��������ʾģ��.����(��).------------------------------
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

    disp('Calculation Complete!');                %���������ʾ               

    %-----------�������߻���ģ��--------------------------------------------
    imag_pic = semilogx(sigma_end, imag(test_V_all), '-r','LineWidth',2);
                                              %�����鲿
    hold on;   
    real_pic = semilogx(sigma_end, real(test_V_all), '-b','LineWidth',2);
                                              %����ʵ��
    grid on;
    xlabel('sigma(S/m)');                     
    ylabel('Voltage(V)');
    legend([imag_pic,real_pic],'Imag','Real');%����ͼ��
    title(['Loggin: ', num2str(kTotalLayerNum), ' Layers Model']);
                                              %�������߱���
    scrsz = get(0,'ScreenSize');              %��ȡ��Ļ�ֱ���
    set(gcf,'position', ...
        0.6*scrsz+0.1*[scrsz(3), scrsz(3),scrsz(4),scrsz(4)]);           
                                              %���û�ͼ����ȫ��
    saveas(gcf,['loggin_', num2str(kTotalLayerNum), '_layers_V.jpg'])
                                              %�洢ͼ���ļ�
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    %----------���ݴ洢�ļ�����---------------------------------------------
    char_file_name = ['loggin_', num2str(kTotalLayerNum), '_layers_data.txt'];
    data_file_handle = fopen(char_file_name,'w');

    fprintf(data_file_handle, '############˵����ʼ##############\n');
    fprintf(data_file_handle, '#���ĵ���func_loggin_multilayer.m���ɵ������ĵ�\n');
    fprintf(data_file_handle, '#�����Ը����Լ���������ȡ����\n');
    fprintf(data_file_handle, '#��ͬ���ݵı�ʶ��Ϊ <==:datakind:==> \n');
    fprintf(data_file_handle, '#��ͬ��datakind����ͬ������\n');
    fprintf(data_file_handle, '#����Ҫ�Լ������ȡ���ĵ��е�����\n');
    fprintf(data_file_handle, '####copyleft---liyang---2017####\n');
    fprintf(data_file_handle, '############˵������##############\n\n\n\n');
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
  
  use_time = toc;                            %�����ʱ��, ��ʱ����
  disp('Timer End.');
  disp(['Total Time: ', num2str(use_time),'s']);
  disp(' ');
  disp('Program Complete!');                 %���������ʾ
end%����������

%=======================================================
%%%%%%%%%%%%%%%%%%������������%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=======================================================
%%%----------------ֱ�Ӽ���V_dir����ֵ------------------%%%
function V_dir = cal_V_dir()
%�ɾ�ȷ��ʽ, ֱ�ӵõ�V_dirֵ�ĺ���
%����ʹ�õ��Ǽ򵥵ľ��λ��ַ�

  global omega mu m I_T rho_R rho_T z_R sigma kPickIntervalNumVdir
  
  phi = linspace(0, pi, kPickIntervalNumVdir+1);      %��ȡȡ����������
  
  first_phi = phi;
  last_phi = phi;
  first_phi(end) = [];
  last_phi(1) = [];
  delta_phi = last_phi - first_phi;                   %�����������֮��ľ���
  
  R = sqrt(rho_T^2 + rho_R^2 - 2*rho_R*rho_T .* cos(phi) + z_R^2);  
  k_dir = sqrt(1i * omega * mu(m) * sigma(m));
  f_phi = exp(1i*k_dir*R) ./ R .* cos(phi);           %������������ʽ
  
  first_f_phi = f_phi;
  last_f_phi = f_phi;
  first_f_phi(end) = [];
  last_f_phi(1) = [];
  average_f_phi = (first_f_phi + last_f_phi) / 2;     %�������������ֵ
  
  inte_f_phi = sum(delta_phi .* average_f_phi);       %�����������ڻ�
  
  V_dir = 1i * omega * mu(m) * I_T * rho_R * rho_T * inte_f_phi;  
end



%%%------���õ糡����kz�ĺ���,���ø�˹���ַ�����V_ref------%%%
function V_ref = cal_V_ref()
%����������
  global rho_R kPickIntervalNumVref
  k_z_min = 10^(-10);
  k_z_max = 10^7;
  k_z = get_interval_split(k_z_min, k_z_max , kPickIntervalNumVref);
                                           %����ȡ������
  
  inte_V = 0;                              %��ʼ���������ӵ�ֵ
  for interval_num = 1:kPickIntervalNumVref
    inte_V = inte_V + ... 
             cal_inte_cosE_kz(k_z(interval_num), k_z(interval_num+1));
                                           %�ۼӼ������                                       
  end
  
  V_ref = 4 * pi * rho_R * inte_V;
end

%%%%%%%%

function inte_cosE_kz = cal_inte_cosE_kz(k_z_j, k_z_jplus1)
%ʹ�ø�˹���õ¼���ÿһ��С����Ļ���ֵ
  global z_R t W;
  
  k_z = (k_z_jplus1 + k_z_j)/2 + ((k_z_jplus1 - k_z_j)/2) * t;  
                                           %�������t��Ӧ��k_z��ֵ
  g_kz = array_cal_E_phi_kz(k_z) .* cos(k_z*z_R);
                                           %����k_zֵ��Ӧ�ı�������ֵ
  inte_cosE_kz = (k_z_jplus1 - k_z_j)/2 * sum(g_kz .* W);    
                                           %����򻯺�Ļ���ֵ
end
  
%%%%%%%%

function interval_split = get_interval_split(lower, upper, split_num)
%��ȡ�����������ĺ���
  k_1 = lower;                             %������½�
  k_max = upper;                           %������Ͻ�
  delta_k = log(k_max/k_1)/split_num;      %�������������
  
  interval_split = zeros(1,split_num+1);   %��ʼ����������ݴ洢����
  for interval_split_index=1:split_num+1 
    interval_split(interval_split_index) = ...
          k_1*exp((interval_split_index-1)*delta_k);
  end 
end

%%%-----------------����phi����糡ֵ-----------------------%%%
function array_E_phi_kz = array_cal_E_phi_kz(array_kz)
%֧������k_z�������E����������ݴ��Ϊ�����ģ��
  global kUseFMinusI
  
  array_size = length(array_kz);          %��ȡ��������ĳ���, 5
  
  array_E_phi_kz = zeros(1, array_size);  %��ʼ�����ݴ洢����
  
  for kz_index = 1:array_size
    array_E_phi_kz(kz_index) = cal_E_phi_kz(array_kz(kz_index), kUseFMinusI);
  end                                     %�������E����ֵ
end

%%%%%%%%

function E_phi_kz = cal_E_phi_kz(k_z, use_F_minus_I)
%����糡��ֵ�ĺ���
  global rho_T rho_R I_T omega mu m
  
  gamma_m = cal_gamma_n(m, k_z);
         
  c_m = 1i * I_T * rho_T / 4.0 * gamma_m * [0;1];%�õ�c_m������  
  
  if rho_R >= rho_T 
    Fm_phi_plus = get_F_phi(use_F_minus_I, +1, gamma_m, k_z);
    E_phi_kz = [0, 1i*omega*mu(m)/gamma_m] * (Fm_phi_plus * c_m *...
                cal_multiply_besseljh1(1, 1, gamma_m*rho_T, gamma_m*rho_R));
  else
    Fm_phi_minus = get_F_phi(use_F_minus_I, -1, gamma_m, k_z);
    E_phi_kz = [0, 1i*omega*mu(m)/gamma_m] * (Fm_phi_minus * c_m *...
                cal_multiply_besseljh1(1, 1, gamma_m*rho_R, gamma_m*rho_T));
  end                                             %�ж�ʹ���ĸ���ʽ����糡E 
end

%%%%%%%%

function F_phi_pm = get_F_phi(use_F_minus_I, sign, gamma_m, k_z)
%����F����ĺ���
  global m rho_layer rho_R rho_T
  
  %ͨ�����ú���, �������͸�䷴��ϵ������
  %GAMMA_out�а�˳��洢��GAMMA_12, GAMMA_23,... ,GAMMA_(N-1)N
  %GAMMA_in�а�˳��洢��GAMMA21, GAMMA_32 ,... , GAMMA_N(N-1) 
  %T����
  [GAMMA_out, T_out, GAMMA_in, T_in] = cal_local_reflect_factor(k_z);
  
  %�õ�������������巴��ϵ����ֵ
  tilde_GAMMA_m_mplus1 = cal_general_reflect_factor...
                         (m, m+1, k_z, GAMMA_out, T_out, GAMMA_in, T_in);
  tilde_GAMMA_m_mminus1 = cal_general_reflect_factor...
                         (m, m-1, k_z, GAMMA_out, T_out, GAMMA_in, T_in); 
                                                  
  if sign == 1
    %%����Fphi+��ֵ  
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
    %%����Fphi-��ֵ         
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
    disp('F_phi_pluse �� F_phi_minus�����������,����ѡȡ����ֻ����+1��-1.');
    return;
  end
  
  if use_F_minus_I
    F_phi_pm = F_phi_pm - eye(2);           %��ȥ��λ�����������Ե�ѹ��С
  end
end

%%%%%%%%

function general_reflect_factor = cal_general_reflect_factor...
         (layer_from, layer_to, k_z, GAMMA_out, T_out, GAMMA_in, T_in)
%������巴��ϵ���ĺ���
  
  global rho_layer total_layer_num
  
  if layer_from>layer_to 
    general_reflect_factor = GAMMA_in{1};         %��ʼ���ݹ���ֵ
    
    for n = 2:layer_to
      gamma_n = cal_gamma_n(n, k_z); 
      
      mmminus1_frac_besselh1_multiply_frac_besselj = ...   
              cal_frac_besselh1(0, 0, gamma_n*rho_layer(n), ...
                                      gamma_n*rho_layer(n-1))...
              *...
              cal_frac_besselj(0, 0, gamma_n*rho_layer(n-1), ...
                                     gamma_n*rho_layer(n));
                                                  %�����ظ�����, �Ż������ٶ�
      
      %��һ��͸��ϵ��
      general_transmit_factor = ...
          (eye(2) - GAMMA_out{n} * general_reflect_factor *...
                    mmminus1_frac_besselh1_multiply_frac_besselj)\T_in{n};                       
      
      %��һ������ϵ��
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
      
      %��һ����͸��ϵ��
      general_transmit_factor = ...
          (eye(2) - GAMMA_in{n} * general_reflect_factor * ...
                    mmplus1_frac_besselh1_multiply_frac_besselj)\T_out{n};                    
      
      %��һ������ϵ��
      general_reflect_factor = ...
          T_in{n} * general_reflect_factor * general_transmit_factor *...
          mmplus1_frac_besselh1_multiply_frac_besselj + GAMMA_out{n};
    end
  end
end

%%%%%%%%
function [GAMMA_out, T_out, GAMMA_in, T_in] = cal_local_reflect_factor(k_z)
%��������͸��ϵ���ĺ���  
  global total_layer_num
  
  T_out = cell(1,total_layer_num-1);
  GAMMA_in = cell(1,total_layer_num-1);
  GAMMA_out = cell(1,total_layer_num-1);
  T_in = cell(1,total_layer_num-1);         %��ʼ�����巴��͸��ϵ�������cell����
  
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
%����j��h����ĺ���
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
%����gamma_n�ĺ���
  global k
  gamma_n = sqrt(k(n)^2 - k_z^2);
    if imag(gamma_n) < 0
      gamma_n = -gamma_n;
    end                                         %gamma_nֻȡ�鲿�������Ľ�
end

%%%%%%%%

function frac_besselj = cal_frac_besselj(up_n, down_n, up_x, down_x)
%�����һ��bessel�����ı�ֵ��ʽ��ֵ
%up_x, up_n�� ���� �ϵ�һ��bessel�����Ľ����ͺ����Ա���ֵ
%down_x, down_n�� ��ĸ �ϵ�һ��bessel�����Ľ����ͺ����Ա���ֵ

  if up_n~=1 && down_n~=1 && up_n~=0 && down_n~=0
      disp('cal_frac_besselj��֧��1,0֮���n������!');
      return;
  end
  
  %��������ؼ����鲿, ��600������ʵ�ʲ��Ե�.
  %�˴�ȡabs��ԭ�����, bessel�������Ա�����ģ��������ʱ����Ľ�����ʽ
  %��eָ���ϵ�ָ������600ʱ�ᵼ�µ��������,�Լ�bessel����ͬʱ�ڷ��ӷ�ĸ��
  %�ۺϷ����õ�, ��������
  if abs(imag(up_x))>600 || abs(imag(down_x))>600                                                
    frac_besselj = sqrt(down_x/up_x)*exp(1i*(down_x-up_x));  
    if down_n>up_n
      frac_besselj = -1i * frac_besselj;    
    elseif down_n<up_n
      frac_besselj = 1i * frac_besselj;   
    end                                     %���bessel�����Ա����鲿abs����600, 
                                            %��ʹ�ý��Ʒ�����
  else
    frac_besselj = besselj(up_n, up_x)/besselj(down_n, down_x);    
  end                                       %���bessel�����Ա����鲿absС��600, 
                                            %��ֱ�ӵ����ڲ���������
end

function frac_besselh1 = cal_frac_besselh1(up_n, down_n, up_x, down_x)
%���㺺�˶������ı�ֵ��ʽ��ֵ
  if up_n~=1 && down_n~=1 && up_n~=0 && down_n~=0
      disp('cal_frac_besselh��֧��1,0֮���n,m������!');
      return;
  end
  
  %��cal_frac_besselj��������·��ͬ
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
%���㺺�˶��������Ե�һ�౴����������ֵ
  if j_n~=1 && h1_n~=1 && j_n~=0 && h1_n~=0
      disp('cal_multiply_besseljh1��֧��1,0֮���n,m������!');
      return;
  end
  
  %��cal_frac_besselj��������·��ͬ
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

%%%-----------------������������-----------------------%%%
function device_layer = get_device_layer()
%����װ��λ����һ������еĺ���
  global rho_layer rho_R rho_T total_layer_num
  
  for search_index = 1:(total_layer_num-1)
    if rho_layer(search_index) > rho_R
      device_layer_R = search_index;        %̽��λ����һ��������
      break;
    else
      device_layer_R = total_layer_num;     %������ڷֽ������ݷ�Χ��, ��ôһ��λ�������
    end
  end
  
  for search_index = 1:(total_layer_num-1)
    if rho_layer(search_index) > rho_T
      device_layer_T = search_index;        %̽��λ����һ��������
      break;
    else
      device_layer_T = total_layer_num;     %������ڷֽ������ݷ�Χ��, ��ôһ��λ�������
    end
  end
  
  if device_layer_R ~= device_layer_T
    disp('̽����Ȧ�뾶�������Ȧ���򲻷�, �޷�����!');
    disp('��Ҫ�������, ���޸�get_device_layer()����');
    return;                                   %��ֹ��������
  end
  
  device_layer = device_layer_T;
end

%============================================================
%--------------------------������־---------------------------
%============================================================
%--2017.11.5
%  -1- �����д��ɲ����е�һ������. ���ִ�������.
%  -2- �����˳�����ֵ������﷨����, �������������, ��������ľ���
%
%--2017.11.6
%  -1- ʹ��bessel����������Զ���Ľ���, �����˾����������������
%  -2- ��������������������
%
%--2017.11.6
%  -1- ���ֶ������Ķ��������, û��ʹ����ν��[F-I]�㷨
%  -2- ����ÿ�ζ���ȫ��һ�����ʮ�ַ�ʱ, ��˴�����'���ݲ���ģʽ'
%  -3- ʵ��[F-I]�㷨��, �ɹ����������
%
%--2017.11.8
%  -1- ���㷨�����Ż�, ɾ������������㲽��, �ٶ������2������
%  -2- ����˻�ӭ��������ɸ�������(copyleft)
%  -3- ����processing bar
%  -4- ������ʱϵͳ
%  -5- 3�����ʱ��: 21.9822s; 5�����ʱ��: 43.3241s
%
%--2017.11.9
%  -1- Ϊ�˱�֤��ͼ��ƽ��������, ��40������������Ϊ80��
%  -2- Ϊ�˱�֤���������ȶ���, ��V_ref���ֲ�������550������Ϊ1000�� 
%  -3- 3�����ʱ��: 90.979s; 5�����ʱ��: 153.0233s 
%
%--2017.11.11
%  -1- ��ͼͨ����rho_T, rho_R��ֵ��һ��С�����ֶ�, ʹE_phi��������, 
%      ����ֱ����F����V_all
%  -2- ֱ�ӷ�������, ��ȻE��ͨ��������rho_Rʹ��Ѹ������, ��E_phi����ͼ��ȷ
%      ��������, ��ʹ��[F-I]���������E_phi���鲿ԶԶ����ʵ��(��Լ��������),
%      ��һֱ����ʵ��. ���,����ֵ���ɱ���Ļ���ֽϴ����. Ȼ������ʵ���ǲ���
%      ��!(��������?�����ֵ��������?)
%
