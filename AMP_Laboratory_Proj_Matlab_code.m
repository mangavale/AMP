% Advanced Manufacturing Processes Laboratory
tic
clear all
close all
clc

set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')
%% Thermophysical properties

k=0.3; %Thermal conductivity (W/mK)
rho=1450; %Density (kg/m^3)
cp=2500; %Specific heat capacity at (J/kgK)
alpha=k/(rho*cp); %Thermal diffusivity (m^2/s)
Tvap=300; %Vaporization temperature (°C)
T0=25; % ambient temperature (°C)

%% Model parameters

eta0=0.18;%0.3; %efficiency coefficient (start with 1 or a reasonable estimated)
c_z=20;%10; % Geometrical correction factor
%% Process parameters

P=1.6; % Power (W)
v=800; % (mm/min) fade rate
v=v/60; %speed (mm/s)
v=v*1e-3; % (m/s)

%% Experimental data
Track_w=122; %(um) insert here track width experimentally measured value

%% Spatial domain
%xmin=-1000; %(um) here x refers to the relative spatial coordinate!!!
xmin=-230/2; %(mm)
%xmax=50; %(um)
xmax=230/2; %(mm)
dx=1;%(um) 10before
%ymax=200;%(um)
ymax=250/2;%(mm)
ymin=-ymax;%(um)
dy=dx;%(um)

zmin=0;
zmax=20;
dz=10;
% zmin=0;
% zmax=0;
% dz=0;

xmin=xmin*1e-6;
xmax=xmax*1e-6;
dx=dx*1e-6;
ymin=ymin*1e-6;
ymax=ymax*1e-6;
dy=dy*1e-6;
zmin=zmin*1e-6;
zmax=zmax*1e-6;
dz=dz*1e-6;

xv=[xmin:dx:xmax];
yv=[ymin:dy:ymax];
zv=[zmin:dz:zmax];
nx=length(xv);
ny=length(yv);
nz=length(zv);
if nz==0
   zv=[0]; nz=length(zv); 
end
if nx==0
   xv=[0]; nx=length(xv); 
end
if ny==0
   yv=[0]; ny=length(yv); 
end

%% Calculate

%Numerical
Tmat=zeros(nx,ny,nz); %intialize matrix
TmatLIM=Tmat;
for ii=1:nx
    for jj=1:ny
        for kk=1:nz
            r=sqrt(xv(ii)^2+yv(jj)^2+(c_z*zv(kk))^2); %calculate radius in specific position
            Tmat(ii,jj,kk)=T0+eta0*(P/(2*pi*k*r))*exp((v/(2*alpha))*(-xv(ii)-r)); % calculate solution               
            if Tmat(ii,jj,kk)>Tvap %Temperature limited solution
                TmatLIM(ii,jj,kk)=Tvap; 
            else
                TmatLIM(ii,jj,kk)= Tmat(ii,jj,kk);
            end
            %pause
        end
    end
    fprintf('%d',ii)
    if mod(ii,20)==0
       fprintf('\n') 
    end
end
%% Plot temperature field - Example of calculated temperature fields

xp=xv*1e6;
yp=yv*1e6;
zp=zv*1e6;
TP1=Tmat; % temperature limited solution
TP2=TmatLIM; %
TP_XY=TP1(:,:,1);
TP_XY2=TP2(:,:,1);
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
surf(yp,xp,TP_XY)
ylabel('\xi (\mum)')
xlabel('y (\mum)')
zlabel('T (°C)')
title('T(\xi,y)')
cbar=colorbar; cbar.Title.String='T (°C)';
%shading interp
axis tight
subplot(1,2,2)
surf(yp,xp,TP_XY2)
ylabel('\xi (\mum)')
xlabel('y (\mum)')
zlabel('T (°C)')
cbar=colorbar; cbar.Title.String='T (°C)';
title('T(\xi,y) - T limited solution')
%shading interp
axis tight

%% Plot temperature field - Calculate track width
FSIZE=7;% fontsize of contour labes
figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,3,1)
surf(yp,xp,TP_XY2)
ylabel('\xi (\mum)')
xlabel('y (\mum)')
zlabel('T (°C)')
cbar=colorbar; cbar.Title.String='T (°C)';
title('T(\xi,y) - T limited solution')
shading interp
axis image
view(0,90)

subplot(1,3,2)
[CV, HV]=contour(yp,xp,TP_XY2);
clabel(CV,HV,'FontSize',FSIZE)
hold on
tit_txt=sprintf('Contour plot');
title(tit_txt)
ylabel('\xi (\mum)')
xlabel('y (\mum)')
%ylim([-1000 xmax*1e6])
axis image

subplot(1,3,3)
[CV, HV]=contour(yp,xp,TP_XY2,[Tvap Tvap],'r');
clabel(CV,HV,'FontSize',FSIZE)
hold on
%Calculate VAP width
[WmaxV,Ind_max]=max(CV(1,2:end));
WmaxV=WmaxV*2;
wLineX=[-CV(1,Ind_max+1), CV(1,Ind_max+1)];
wLineY=[CV(2,Ind_max+1), CV(2,Ind_max+1)];
plot(wLineX,wLineY,'LineWidth',2)
tit_txt=sprintf('VAP width: %.0f \\mum',WmaxV);
title(tit_txt)
ylabel('\xi (\mum)')
xlabel('y (\mum)')
axis image

%% Plot -- Efficiency calibration
eta_test1=0.1;
eta_test2=0.9;
Tstar=(Tmat-T0)/eta0;
Tmat_eta1=Tstar*eta_test1+T0;
Tmat_eta2=Tstar*eta_test2+T0;

Teta1=squeeze(Tmat_eta1(:,:,1));
Teta2=squeeze(Tmat_eta2(:,:,1));

figure
subplot(2,2,1)
surf(yp,xp,Teta1)
ylabel('\xi (\mum)')
xlabel('y (\mum)')
zlabel('T (°C)')
cbar=colorbar; cbar.Title.String='T (°C)';
caxis([T0, Tvap]);
tittxt=sprintf('T_{lim}(\\xi,y)\n\\eta=%.2f',eta_test1);
title(tittxt)
shading interp
axis image
view(0,90)

subplot(2,2,2)
[CV, HV]=contour(yp,xp,Teta1,[Tvap Tvap],'r');
clabel(CV,HV,'FontSize',FSIZE)
hold on
%Calculate VAP width
[WmaxV,Ind_max]=max(CV(1,2:end));
WmaxV=WmaxV*2;
wLineX=[-CV(1,Ind_max+1), CV(1,Ind_max+1)];
wLineY=[CV(2,Ind_max+1), CV(2,Ind_max+1)];
plot(wLineX,wLineY,'LineWidth',2)
tit_txt=sprintf('VAP width: %.0f \\mum',WmaxV);
title(tit_txt)
ylabel('\xi (\mum)')
xlabel('y (\mum)')
axis image

subplot(2,2,3)
surf(yp,xp,Teta2)
ylabel('\xi (\mum)')
xlabel('y (\mum)')
zlabel('T (°C)')
cbar=colorbar; cbar.Title.String='T (°C)';
caxis([T0, Tvap]);
tittxt=sprintf('T_{lim}(\\xi,y)\n\\eta=%.2f',eta_test2);
title(tittxt)
shading interp
axis image
view(0,90)

subplot(2,2,4)
[CV, HV]=contour(yp,xp,Teta2,[Tvap Tvap],'r');
clabel(CV,HV,'FontSize',FSIZE)
hold on
%Calculate VAP width
[WmaxV,Ind_max]=max(CV(1,2:end));
WmaxV=WmaxV*2;
wLineX=[-CV(1,Ind_max+1), CV(1,Ind_max+1)];
wLineY=[CV(2,Ind_max+1), CV(2,Ind_max+1)];
plot(wLineX,wLineY,'LineWidth',2)
tit_txt=sprintf('VAP width: %.0f \\mum',WmaxV);
title(tit_txt)
ylabel('\xi (\mum)')
xlabel('y (\mum)')
axis image



%% Calibration of efficiency coefficient

etav=[0.01:0.01:1];
nv=length(etav);
errv=zeros(nv,1);
Wtmp=zeros(nv,1);
Tstar=(Tmat-T0)/eta0;
f=figure;
for kk=1:nv

    Tmat_tmp=Tstar*etav(kk)+T0;
    Tmat_tmp2=squeeze(Tmat_tmp(:,:,1));
   [CVAP, hVAP]=contour(yp,xp,Tmat_tmp2,[Tvap Tvap]); 
   [Wtmp(kk),Ind_max]=max(CVAP(1,2:end));
    Wtmp(kk)=Wtmp(kk)*2;
    errv(kk)=100*abs(Track_w-Wtmp(kk))/Track_w;
    tit_txt=sprintf('\\eta= %.2f',etav(kk));
    title(tit_txt)
    legend('Vaporisation iso-therm','Location','South')
    axis image
    
end

[minerr, idx_minerr]=min(errv);
Wtmp(idx_minerr);
Tmat_CAL=Tstar*etav(idx_minerr)+T0;
[CHZ, hHZ]=contour(yp,xp,squeeze(Tmat_CAL(:,:,1)),[Tvap Tvap]);
    axis image
    tit_txt=sprintf('\\eta= %.2f',etav(idx_minerr));
    title(tit_txt)
[WVAP,Ind_VAP]=max(CHZ(1,2:end));
WVAP=WVAP*2;
fprintf('\n\n')
fprintf('----------------------------------\n')
fprintf('Track Width (Efficiency=%.2f): %.3fum\n',etav(idx_minerr),Wtmp(idx_minerr))
fprintf('----------------------------------\n')
%% Plot efficiency coefficient calibration
figure
yyaxis left
plot(etav,errv)
ylabel('Relative error, {\iterr} (%)')
yyaxis right
plot(etav,Wtmp,[min(etav) max(etav)],[Track_w Track_w])
ylabel('Track width, {\itw} (\mum)')
xlabel('Efficiency, \eta')
set(gca,'FontName','Times New Roman','FontSize',12)
legend('{\iterr}','{\itw_{sim}}','{\itw_{meas}}')
set(gcf,'color','w');

%% Plot -- Penetration depth estimation


AXIM=0; %set to 1 to have correct scaling
FSIZE=7;% fontsize of contour labes
% nrow=2;
% ncol=3;

nrow=4;
ncol=4;
xp=xv*1e6;
yp=yv*1e6;
zp=zv*1e6;
%TP1=TmatLIM; % non-calibrated temperature limited solution
%TP1=Tmat; % non calibrated solution
Tstar=(Tmat-T0)/eta0; %calibrated solution
TP1=Tstar*etav(idx_minerr)+T0;
TP_XY=TP1(:,:,1);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(nrow,ncol,[1 5])
surf(yp,xp,TP_XY)
shading interp
view(0,90)
if AXIM==1
    axis image
else
    axis tight
end
ylabel('\xi (\mum)')
xlabel('y (\mum)')
zlabel('T (°C)')
cbar=colorbar; cbar.Title.String='T (°C)';
caxis([T0 Tvap])
title('T_{lim}(\xi,y)')

subplot(nrow,ncol,[9 13])
[CV, HV]=contour(yp,xp,TP_XY,[Tvap Tvap],'r');
clabel(CV,HV,'FontSize',FSIZE)
hold on
%Calculate VAP width
[WmaxV,Ind_max]=max(CV(1,2:end));
WmaxV=WmaxV*2;
wLineX=[-CV(1,Ind_max+1), CV(1,Ind_max+1)];
wLineY=[CV(2,Ind_max+1), CV(2,Ind_max+1)];
plot(wLineX,wLineY,'LineWidth',2)
tit_txt=sprintf('VAP width: %.0f \\mum',WmaxV);
title(tit_txt)
ylabel('\xi (\mum)')
xlabel('y (\mum)')
if AXIM==1
    axis image
else
    axis tight
end


pos1=ceil(ny/2);
TP_XZ=TP1(:,pos1,:);
%subplot(nrow,ncol,2)
subplot(nrow,ncol,[2 3 4])
surf(xp,zp,squeeze(TP_XZ)')
shading interp
view(0,-90)
if AXIM==1
    axis image
else
    axis tight
end
xlabel('\xi (\mum)')
ylabel('z (\mum)')
zlabel('T (°C)')
cbar=colorbar; cbar.Title.String='T (°C)';
caxis([T0 Tvap])
title('T_{lim}(\xi,z)')

%subplot(nrow,ncol,5)
subplot(nrow,ncol,[6 7 8])
[CV, HV]=contour(xp,-zp,squeeze(TP_XZ)',[Tvap Tvap],'r');
clabel(CV,HV,'FontSize',FSIZE)
hold on
%Calculate VAP width
[HmaxV,Ind_max]=min(CV(2,2:end));
wLineX=[CV(1,Ind_max+1), CV(1,Ind_max+1)];
wLineY=[CV(2,Ind_max+1), 0];
plot(wLineX,wLineY,'LineWidth',2)
tit_txt=sprintf('VAP depth: %.0f \\mum',abs(HmaxV));
title(tit_txt)
xlabel('\xi (\mum)')
ylabel('z (\mum)')
if AXIM==1
    axis image
else
    axis tight
end

pos0=find(xp==CV(1,Ind_max+1)); % extract width profile where highest depth is found
TP_YZ=TP1(pos0,:,:);
%subplot(nrow,ncol,3)
subplot(nrow,ncol,[10 11 12])
surf(yp,zp,squeeze(TP_YZ)')
shading interp
view(0,-90)
if AXIM==1
    axis image
else
    axis tight
end
xlabel('y (\mum)')
ylabel('z (\mum)')
zlabel('T (°C)')
cbar=colorbar; cbar.Title.String='T (°C)';
caxis([T0 Tvap])
title('T_{lim}(y,z)')

%subplot(nrow,ncol,6)
subplot(nrow,ncol,[14 15 16])
[CV2, HV2]=contour(yp,-zp,squeeze(TP_YZ)',[Tvap Tvap],'r');
clabel(CV2,HV2,'FontSize',FSIZE)
hold on
%Calculate VAP width
[HmaxV,Ind_max]=min(CV2(2,2:end));
wLineX=[CV2(1,Ind_max+1), CV2(1,Ind_max+1)];
wLineY=[CV2(2,Ind_max+1), 0];
plot(wLineX,wLineY,'LineWidth',2)
tit_txt=sprintf('VAP depth: %.0f \\mum',abs(HmaxV));
title(tit_txt)
xlabel('y (\mum)')
ylabel('z (\mum)')
yt = yticks;
if AXIM==1
    axis image
else
    axis tight
end

fprintf('\nPenetration depth h: %.2f',abs(HmaxV))
fprintf('\n----------------------------------\n')
