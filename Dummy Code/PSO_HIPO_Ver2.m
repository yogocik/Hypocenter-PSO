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
LB=[198 398 999] ;    %Batas Bawah x,y,z
UB=[201 400 1010];    %Batas Atas x,y,z

% pso parameters values
m=3; % number of variables
n=100; % swarm size
wmax=0.9; % inertia weight
wmin=0.4; % inertia weight
c1=0.5+log(2); % acceleration factor 1
c2=c1; % acceleration factor 2

% pso main program start
maxite=1000; % set maximum number of iteration
maxrun=10; % set maximum number of runs need to be
for run=1:maxrun
% pso initialization start
for i=1:n
for j=1:m
x0(i,j)=round(LB(j)+rand()*(UB(j)-LB(j)));   %data dugaan hiposenter, distribusi ae
end
end

x=x0; % initial swarm
v=zeros(x0); % initial velocity

for i=1:n
f0(i,1)=(to+(sqrt(((x_hipo-x(i))^2+(y_hipo-y(i))^2+(z_hipo-z(i))^2))/v_p))+...
        (k*(to+(sqrt(((x_hipo-x(i))^2+(y_hipo-y(i))^2+(z_hipo-z(i))^2))/v_p))); %Waktu arrival distribusi hiposenter
end
[fmin0,index0]=min(f0);
pbest=x0; % initial pbest
gbest=x0(index0,:); % initial gbest

% pso initialization end

% pso algorithm start
ite=1;
tolerance=1;
while ite<=maxite && tolerance>10^-12
w=wmax-(wmax-wmin)*ite/maxite; % update inertial weight
% pso velocity updates
for i=1:n
for j=1:m
v(i,j)=w*v(i,j)+c1*rand()*(pbest(i,j)-x(i,j))...
+c2*rand()*(gbest(1,j)-x(i,j));
end
end
% pso position update
for i=1:n
for j=1:m
x(i,j)=x(i,j)+v(i,j);
end
end
% handling boundary violations
for i=1:n
for j=1:m
if x(i,j)<LB(j)
x(i,j)=LB(j);
elseif x(i,j)>UB(j)
x(i,j)=UB(j);
end
end
end
% evaluating fitness
for i=1:n
f(i,1)=(to+(sqrt(((x_hipo-x(i))^2+(y_hipo-y(i))^2+(z_hipo-z(i))^2))/v_p))+...
        (k*(to+(sqrt(((x_hipo-x(i))^2+(y_hipo-y(i))^2+(z_hipo-z(i))^2))/v_p)));
end
% updating pbest and fitness
for i=1:n
if f(i,1)<f0(i,1)
    pbest(i,:)=x(i,:);
f0(i,1)=f(i,1);
end
end
[fmin,index]=min(f0); % finding out the best particle
ffmin(ite,run)=fmin; % storing best fitness
ffite(run)=ite; % storing iteration count
% updating gbest and best fitness
if fmin<fmin0
gbest=pbest(index,:);
fmin0=fmin;
end
% calculating tolerance
if ite>100;
tolerance=abs(ffmin(ite-100,run)-fmin0);
end
% displaying iterative results
if ite==1
disp(sprintf('Iteration Best particle Objective fun'));
end
disp(sprintf('%8g %8g %8.4f',ite,index,fmin0));
ite=ite+1;
end
% pso algorithm-----------------------------------------------------end
gbest;
fvalue=10*(gbest(1)-1)^2+20*(gbest(2)-2)^2+30*(gbest(3)-3)^2;
fff(run)=fvalue;
rgbest(run,:)=gbest;
disp(sprintf('g_best'));
end

% pso main program------------------------------------------------------end
disp(sprintf('\n'));
disp(sprintf('PSOMainProgram'));
disp(sprintf('Final Results-----------------------------'));
[bestfun,bestrun]=min(fff)
best_variables=rgbest(bestrun,:)
disp(sprintf('Emboh'));
toc
% PSO convergence characteristic
plot(ffmin(1:ffite(bestrun),bestrun),'-k');
xlabel('Iteration');
ylabel('Fitness function value');
title('PSO convergence characteristic')
