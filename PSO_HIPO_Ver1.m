%%SCRIPT HIPOSENTER UAS INVERSI
%ANGGOTA KELOMPOK :Ade, Yolanda, Rusba, Yogic 
%TEKNIK GEOFISIKA
%FAKULTAS TEKNIK SIPIL, LINGKUNGAN, DAN KEBUMIAN
%INSTITUT TEKNOLOGI SEPULUH NOPEMBER 

%%PEMODELAN PUSAT GEMPA BUMI PADA LAPISAN HOMOGEN ISOTROPIS

tic
clear all
clc
close

%Parameter/ Data Stasiun Pengamatan (x,y,z)
x = [ 300 700 1000 200 3500 1800 2000 50];
y = [ 1200 200 450 600 100 100 1000 200 ];
z = rand(size(x)).*1000;
v_p = 10.9;

%NOISE
e =0;                        %Persen Error Maksimum
k =e*randn(1,1);

%Parameter Forward Modelling
x_hipo = 200;y_hipo=400;z_hipo=1000;to= 0;     %Asumsi Letak Pusat Sebenarnya

%Forward Modelling DATA ASUMSI ASLI dengan Input Error
t_obs = zeros(length(x),1);
for i=1:length(x)        
    t_obs(i) = (to+(sqrt(((x_hipo-x(i))^2+(y_hipo-y(i))^2+(z_hipo-z(i))^2))/v_p))+...
        (k*(to+(sqrt(((x_hipo-x(i))^2+(y_hipo-y(i))^2+(z_hipo-z(i))^2))/v_p))); %Data Sintetis
end

%CONSTRAIN/BATASAN KOORDINAT PUSAT GEMPA
LB=[198 398 999] ;    %Batas Bawah 
UB=[201 400 1010];    %Batas Atas

%PSO Parameter
m=3  ; % jumlah koordinat gempa bumi
n=100; % swarm size
wmax=0.9; % inertia weight
wmin=0.4; % inertia weight
c1=2; % acceleration factor 1
c2=c1; % acceleration factor 2

%INISIASI SWARM
swarm = zeros (n,m);  %kolom pertama x , kolom kedua y, kolom ketiga z
for p=1:n;
    for t=1:m
        swarm(p,t)= LB(t)+(rand(1,1)*(UB(t)-LB(t)));
    end
end

%Pbest dan Gbest Awal
Val = zeros (n,1);
Val_update = Val;
t_cal = zeros(size(t_obs));
for p = 1:n;
for i=1:length(x)        
    t_cal = (to+(sqrt(((swarm(p,1)-x(i))^2+(swarm(p,2)-y(i))^2+(swarm(p,3)-z(i))^2))/v_p)); %Data Sintetis
end
Val(p,1)=std((t_obs-t_cal));
end
VALUE_INITIAL = Val
[fmin,index]  = min(Val);
Pbest_awal    = swarm;            %Lokal Minimum Awal
Gbest_awal    = swarm(index,1:m); %Global Minimum Awal

%Looping PSO
MAX_RUN = 20;     %MAIN LOOP 
max_iterasi = n;  %NEST LOOP OR SWARM SIZE CALCULATION
v             = zeros (max_iterasi,m);
error         = 10^-15;
gbest         = Gbest_awal;
pbest         = Pbest_awal;
BEST_VALUE    = zeros(MAX_RUN,m);
Matriks_Error =zeros(MAX_RUN,1);
Iterasi       = Matriks_Error;

for t = 1:MAX_RUN
    if gbest > error;
w=wmin+((wmax-wmin)*((MAX_RUN-t)/t)); % update inertial weight
Val = Val_update; 
Pbest=pbest;
Gbest=gbest
    for i= 1: max_iterasi;
       for g=1:m;
        v(i,g) =w*v(i,g)+c1*rand*(Pbest(i,g)-swarm(i,g))+c2*rand*(Gbest(g)-swarm(i,g));
       end
    end
    swarm = swarm + v;
    for p = 1:max_iterasi; 
      for i=1:length(x)        
          t_cal_update(i) = (to+(sqrt(((swarm(p,1)-x(i))^2+(swarm(p,2)-y(i))^2+(swarm(p,3)-z(i))^2))/v_p)); %Data Sintetis
      end
    Val_update(p)=std((t_obs'-t_cal_update));
    
     %UPDATE Pbest
       if Val_update(p) < Val(p);
          pbest([p]) = swarm(p);
          Val_update(p) = Val(p);
       end
    
        [Gbest_Val,index]           = min(Val_update); %Gbest Indexing
        gbest                       = swarm(index,:);
        Matriks_Error(t)            = min(Val_update);
        
    figure(1)
    clf
        if t<=(MAX_RUN-1);
            plot3(gbest(1),gbest(2),(-1)*gbest(3),'r.','Markersize',15)
            hold on
        else
            plot3(gbest(1),gbest(2),(-1)*gbest(3),'g.','Markersize',15)
            hold on
        end
    plot3(x,y,z,'bv')
    hold on
    title ('Koordinat Hiposenter')
    grid on
    legend ('Hiposenter','Stasiun')     
    end
    end 
    Iterasi(t) = t;
end
figure(2)

plot(Iterasi,Matriks_Error,'b--');
title ('PSO Convergence');

    
        

        
        
        
        
        
        
            
    



