
%%SCRIPT HIPOSENTER UTS INVERSI
%ANGGOTA KELOMPOK :
%TEKNIK GEOFISIKA
%FAKULTAS TEKNIK SIPIL, LINGKUNGAN, DAN KEBUMIAN
%INSTITUT TEKNOLOGI SEPULUH NOPEMBER 

clear all
clc
close

%Parameter/ Data Stasiun Pengamatan (x,y,z)
x = [ 300 700 1000 200 3500 1800];
y = [ 1200 200 450 600 100 100];
z = rand(size(x)).*1000;
v_p = 10.9;

%NOISE
e =0;                        %Persen Error Maksimum
k =e*randn(1,1);

%Parameter Forward Modelling
x_hipo = 200;y_hipo=400;z_hipo=1000;to= 0;     %Asumsi Letak Pusat Sebenarnya

%Forward Modelling DATA ASUMSI ASLI
t_obs = zeros(length(x),1);
for i=1:length(x)        
    t_obs(i) = (to+(sqrt(((x_hipo-x(i))^2+(y_hipo-y(i))^2+(z_hipo-z(i))^2))/v_p))+(k*(to+(sqrt(((x_hipo-x(i))^2+(y_hipo-y(i))^2+(z_hipo-z(i))^2))/v_p))); %Pure
end

% INVERSI DATA AWAL/DUGAAN

g         = 100;
f         = 10^-5;
e_plot    = zeros(g,1);
iterasi   = zeros(g,1);
format long

for iter =1:g;
    
    %Kondisi pemakaian parameter model
    
    if iter == 1;
        x_model =1000;y_model = 300 ;z_model =150;to_model=0;
    else
        x_model  = xo_pertu;
        y_model  = yo_pertu;
        z_model  = zo_pertu;
    end
    
    % DATA arrival time secara teoritis 
    t_cal = zeros (length(x),1);
    for i = 1:length(x)
        t_cal(i) = to + sqrt(((x_model-x(i))^2+(y_model-y(i))^2+(z_model-z(i))^2))/v_p; %t_cal_inversi
    end
    
    % MISFIT & ERROR
    dt_misfit    = t_obs - t_cal;
    error        = std(abs(dt_misfit));
    e_plot(iter) = error;
    
    %ITERASI JACOBI
    if error >= f
        for i=1:length(x)
            deriv_x(i) = (x_model-x(i))/(v_p*sqrt((x_model-x(i))^2 + (y_model-y(i))^2 + (z_model-z(i))^2));
            deriv_y(i) = (y_model-y(i))/(v_p*sqrt((x_model-x(i))^2 + (y_model-y(i))^2 + (z_model-z(i))^2));
            deriv_z(i) = (z_model-z(i))/(v_p*sqrt((x_model-x(i))^2 + (y_model-y(i))^2 + (z_model-z(i))^2));
        end
        
        %Matrifikasi Turunan 
        J = zeros(length(x),3);
        J(:,1) = deriv_x';
        J(:,2) = deriv_y';
        J(:,3) = deriv_z';
      
        %PERTURBASI MODEL 
        dm_perturbasi =inv(J'*J)*J'*dt_misfit;
        xo_pertu = x_model + dm_perturbasi(1);
        yo_pertu = y_model + dm_perturbasi(2);
        zo_pertu = z_model + dm_perturbasi(3);
        iterasi (iter)  = iter;     
    end 
end

%INFORMASI INVERSI

Model_Asumsi         = [x_hipo;y_hipo;z_hipo]
Model_Forward        = [1000;300;150]
Model_Inversi        = [xo_pertu;yo_pertu;zo_pertu]
Model_Perturbasi     = dm_perturbasi
Error_Final          = e_plot(100)
Jumlah_Iterasi       = g
Matriks_Error        = [iterasi e_plot]

%PLOTTING KURVA

figure(1)
plot3(x,y,z,'rv','MarkerFaceColor','b','MarkerEdgeColor','g')
hold on
plot3(xo_pertu,yo_pertu, (-1*zo_pertu),'o', 'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
grid on; 
xlabel('Koord. X (m)'); ylabel('Koord. Y (m)');
zlabel('Kedalaman(m)'); title('Model Prediksi');
legend('Stasiun','Hiposenter');

figure(2)
plot(iterasi,e_plot)
xlabel('Iterasi');ylabel('Standar Deviasi'); title('Grafik Misfit');
