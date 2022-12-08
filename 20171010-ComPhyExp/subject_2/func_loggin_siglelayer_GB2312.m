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
%%%�������: ���Ƚ����е�Ų��⾮�ź�ģ��
%%%���ʱ��: 2017.10.21
%%%Ԥ��ʱ��: 5h
%%%Ŀ��: ʹ�õ��ѧ�����ķ���, ����Matlabģ��⾮�źŵĵ�ѹֵ
%%%����: ��ֵ��������Ȧ֮��Ĵų�
%%%ע��:����matlab�汾����, �е�matlab�汾���ܲ���֧�ֽ��ű���
%%%     ����д��һ��, ֻ��Ҫ���ű��ͺ����ֿ�����, ���ǽ��ű����
%%%     ����.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function func_loggin_siglelayer_GB2312()
  %#########################################
  %%%%%%%%%%%%����ǰ��׼��%%%%%%%%%%%%%%%%%%%%
  %#########################################
  clear;   
  format long; %���ü��㾫��

  %#########################################
  %%%%%%%%%%%%��������ѧ�����ļ���%%%%%%%%%%%%%%
  %#########################################
  global I_T rho S M L mu omega k t W 

  I_T = 1;                               %������Ȧ������С,A (T_T)
  rho = 10^(-3);   %rho_T = rho_R = rho  %����/̽����Ȧ�뾶,m
  S = pi * rho^2;                        %����/̽����Ȧ���
  M = S * I_T;                           %������Ȧ�ž�
 
  L = 0.8;                               %����-̽����Ȧ����,m

  sigma = 1;                             %���ʵ絼��,S/m
  mu_0 = 1.25663706144 * 10^(-6);        %��մŵ���
  mu_r = 1;                              %������Դŵ���
  mu = mu_r * mu_0;                      %���ʴŵ���
  f = 20000;                             %�ź�Ƶ��,Hz
  omega = 2*pi*f;                        %�źŽ�Ƶ��
  k = sqrt(1i * omega * mu * sigma);     %��Ų���ʸ

  t = [-0.9061798459, -0.5384693101,...
        0, ...
        0.5384693101, 0.9061798459];     %��˹�����е����
  W = [0.2369268851, 0.4786286705,...
       0.568888889, ...
       0.4786286705, 0.2369268851];      %��˹�����еı���

  %#########################################
  %%%%%%%%%%%%%%%%%������%%%%%%%%%%%%%%%%%%%%
  %#########################################

  %%�����ַ������ų�ǿ��
  approx_H = cal_approx_H();
  precise_H = cal_precise_H();

  %%�����ַ��������ƴ�С
  %ʹ�����ôų�������
  approx_V_with_approx_H = cal_approx_V(approx_H);
  approx_V_with_precise_H = cal_approx_V(precise_H);
  %ֱ�Ӿ�ȷ������
  precise_V = cal_precise_V();

  %%�������ý��Ϊ���Ǻͷ�ֵ����ʽ
  %��ֵ
  R_precise_V = abs(precise_V);
  R_approx_V_with_approx_H = abs(approx_V_with_approx_H);
  R_approx_V_with_precise_H = abs(approx_V_with_precise_H);
  %����
  theta_precise_V = angle(precise_V)/pi;
  theta_approx_V_with_approx_H = angle(approx_V_with_approx_H)/pi;
  theta_approx_V_with_precise_H = angle(approx_V_with_precise_H)/pi;
  %��ʽ����
  char_exp_precise_V = [num2str(R_precise_V),...
                       '*e^(i��*',num2str(theta_precise_V),')'];
  char_exp_approx_V_with_approx_H = ...
                       [num2str(R_approx_V_with_approx_H),...
                       '*e^(i��*',num2str(theta_approx_V_with_approx_H),')'];
  char_exp_approx_V_with_precise_H = ...
                       [num2str(R_approx_V_with_precise_H),...
                       '*e^(i��*',num2str(theta_approx_V_with_precise_H),')'];

  %%������ý��
  final_print_char = ...
       ['-----------------------------------------------\n',           ...
        '=============���Ƚ��ʲ⾮�ź�(V)==================\n',           ...
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
  %������������Ļ
  fprintf(final_print_char); 
  
  %�������������ݴ洢�ļ�
  loggin_data_file_handle = fopen('loggin_sigle_layers_data.txt','w');
  fprintf(loggin_data_file_handle, final_print_char);
  fclose(loggin_data_file_handle);
end 
%#########################################
%%%%%%%%%%%%�������ܺ���%%%%%%%%%%%%%%%%%%%%
%#########################################
function approx_V = cal_approx_V(H)
%������֪�Ĵų�����Ƶĵ��Ƶĺ���

  global omega mu S;
  approx_V = 1i * omega * mu * S * H; 
end

%%%%%%%%

function precise_V = cal_precise_V()
%�ɾ�ȷ��ʽ, �������ų���С, ֱ�ӵõ����Ƶĺ���
%����ʹ�õ��Ǽ򵥵ľ��λ��ַ�

  global omega mu I_T rho L k
  
  kPhiIntervalSplitNum = 300;                         %������phi����ɢȡֵ����
  
  phi = linspace(0, pi, kPhiIntervalSplitNum+1);      %��ȡȡ����������
  
  first_phi = phi;
  last_phi = phi;
  first_phi(end) = [];
  last_phi(1) = [];
  delta_phi = last_phi - first_phi;                   %�����������֮��ľ���
  
  R = sqrt(rho^2 + rho^2 - 2*rho*rho .* cos(phi) + L^2);  
  f_phi = exp(1i*k*R) ./ R .* cos(phi);               %������������ʽ
  
  first_f_phi = f_phi;
  last_f_phi = f_phi;
  first_f_phi(end) = [];
  last_f_phi(1) = [];
  average_f_phi = (first_f_phi + last_f_phi) / 2;     %�������������ֵ
  
  inte_f_phi = sum(delta_phi .* average_f_phi);       %�����������ڻ�
  
  precise_V = 1i * omega * mu * I_T * rho * rho * inte_f_phi;  
end

%%%%%%%%

function approx_H = cal_approx_H()
%������ƵĴų���С
  global M k L;
  approx_H = (M/(2*pi*L^3)) * (1-1i*k*L) * exp(1i*k*L);
end

%%%%%%%%

function precise_H = cal_precise_H()
%���㾫ȷ�Ĵų���С
  global M;
  
  kPickIntervalPoints = 300;                           %������, ������ĸ���
  k_rho = get_interval_split(kPickIntervalPoints);     %������ķֽ���������
  
  C = get_C(k_rho, kPickIntervalPoints);               %��ȡC��ֵ
  
  inte_H = 0;                                          %��ʼ������ֵ
  n = 3;                                               %����ʽ���е�nֵ
  for C_index = 1:4                                    %l = C_index - 1
    for interval_num = 1:kPickIntervalPoints
      inte_H = inte_H + ...
               C(C_index, interval_num) * ...
               cal_inte_J_0_kn(n, C_index-1, ...
                               k_rho(interval_num), ...
                               k_rho(interval_num+1) );   
                                                       %�ۼӼ�����ֵ�ֵ
    end
  end
  
  precise_H = 1i*M/(4*pi) * inte_H;
end

%%%%%%%%

function inte_J_0_kn = cal_inte_J_0_kn(n, l, k_rho_j, k_rho_jplus1)
%ʹ�ø�˹���õ¼���bessel���������ݴ�n�Ļ���
  global rho t W;
  
  x_j = k_rho_j * rho;
  x_jplus1 = k_rho_jplus1 * rho;                       %��k_rho_jת��Ϊx_j
  
  x = (x_jplus1 + x_j)/2 + ((x_jplus1 - x_j)/2) * t;   %�������t��Ӧ��x��ֵ
  g_x = besselj(0, x) .* (x.^(n+l));                     
  
  inte_J_0_xn = (x_jplus1 - x_j)/2 * sum(g_x .* W);    %����򻯺��bessel����
  
  inte_J_0_kn = (1/rho^(n+l+1)) * inte_J_0_xn;         %������ջ���ֵ
end
  
%%%%%%%%

function C = get_C(k_rho, all_interval_num)
%���ڼ���ϵ��C�ĺ���
  global k L;
  %k_z = sqrt(k^2 - k_rho.^2);
  f_k_rho = exp(1i*sqrt(k^2 - k_rho.^2)*L) ./ sqrt(k^2 - k_rho.^2);        
                                                       %����ֽ������ϵĺ���ֵ
  
  %ʹ�ö�������ͷȥβ�ķ�������
  %����M_j+1 - M_j����ֵ
  first_f_k_rho = f_k_rho;                              
  last_f_k_rho = f_k_rho;                              %�����鸴��Ϊ����
  first_f_k_rho(end) = [];                             %ȥ��β��                      
  last_f_k_rho(1) = [];                                %ȥ��ͷ��
  delta_f_k_rho = last_f_k_rho - first_f_k_rho;        %���������
  
  first_k_rho = k_rho;
  last_k_rho = k_rho;
  first_k_rho(end) = [];
  last_k_rho(1) = [];
  h = last_k_rho - first_k_rho;                         %ȡ����h
  
  d_numberator_split = delta_f_k_rho ./ h;              %d���ʽ�з����е�һ��
  
  first_d_numberator_split = d_numberator_split;
  last_d_numberator_split = d_numberator_split;
  first_d_numberator_split(end) = [];
  last_d_numberator_split(1) = [];
  d_numberator = 6 * (last_d_numberator_split - first_d_numberator_split); 
                                                        %d�ķ���
  
  first_h = h;
  last_h = h;
  first_h(end) = [];
  last_h(1) = [];
  h_add = first_h + last_h;                             %h��λ������õ�ֵ
  
  d = d_numberator ./ h_add;                              
  d = [d d(end)];
  d = [d(1) d];                                         %d��ֵ
  
  mu_Nplus1 = 1;
  mu = first_h ./ h_add;
  mu = [mu mu_Nplus1];                               %mu��ֵ
  
  lambda_1 = 1;
  lambda = last_h ./ h_add;
  lambda = [lambda_1 lambda];                           %lambda��ֵ
  
  M_matrix = diag(repmat(2,[1,all_interval_num+1])) + ...
                  diag(mu, -1) + diag(lambda, +1);     %Mϵ�������ֵ
  
  M_res = M_matrix \ d';
  M_res = M_res';                                       %M�Ľ�
  
  first_M_res = M_res;
  last_M_res = M_res;
  first_M_res(end) = [];
  last_M_res(1) = [];                                   %׼��M����ͷȥβ����
  
  %�������Cϵ����ֵ
  C_3 = (last_M_res - first_M_res) ./ (6 * h); 
  C_2 = (last_k_rho .* first_M_res - first_k_rho .* last_M_res) ./ (2 * h);
  C_1 = ((first_M_res - last_M_res)/6) .* h + ...
        (-last_k_rho.^2 .* first_M_res + first_k_rho.^2 .* last_M_res +...
        2*(last_f_k_rho - first_f_k_rho))./ (2 * h); 
  C_0 = (first_k_rho.*last_M_res - last_k_rho.*first_M_res)/6 .* h +...
        (last_k_rho.^3.*first_M_res - first_k_rho.^3.*last_M_res - ...
        6*(first_k_rho.*last_f_k_rho - last_k_rho.*first_f_k_rho)) ./ (6*h);
  %��Cϵ�����Ϊһ��2D���� 
  C = [];
  C = [C;C_0];
  C = [C;C_1];
  C = [C;C_2];
  C = [C;C_3];
end

%%%%%%%%

function interval_split = get_interval_split(split_num)
%��ȡ�����������ĺ���
  k_1 = 10^(-6);                           %������½�
  k_max = 10^8;                            %������Ͻ�
  delta_k = log(k_max/k_1)/split_num;      %�������������
  
  interval_split = zeros(1,split_num+1);   %��ʼ����������ݴ洢����
  for interval_split_index=1:split_num+1 
    interval_split(interval_split_index) = ...
          k_1*exp((interval_split_index-1)*delta_k);
  end 
end
