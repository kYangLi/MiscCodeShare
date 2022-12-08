m_set = -5:0.02:5;
C_set = []; 

fid = fopen('res_C.txt','w');
for m = m_set
    func = @(x,y) 1/(4*pi)*...
                  (m*cos(x).*cos(y)-cos(x)-cos(y))./ ...
                  ((2-2*m*(cos(x)+cos(y))+2*cos(x).*cos(y)+m^2).^(3/2));
    C_set(end+1) = dblquad(func,-pi,pi,-pi,pi);
    fprintf(fid,"%f    %f\n",m,C_set(end)); 
end    
fclose(fid);

plot(m_set,C_set,'LineWidth',3);
grid;
xlabel('m');
ylabel('C');