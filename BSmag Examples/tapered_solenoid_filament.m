%---------------------------------------------------
%  NAME:      tapered_solenoid_filament.m
%  WHAT:      Calculation of the magnetic field due to two tapered solenoids
%             on a volume (+ 3D plot).
%  REQUIRED:  BSmag Toolbox 20150407
%  AUTHOR:    20150407, L. Queval (loic.queval@gmail.com)
%  COPYRIGHT: 2015, Loic Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause).
%----------------------------------------------------

% Initialize
clear all,
BSmag = BSmag_init(); % Initialize BSmag analysis (first coil)
BSmag1 = BSmag_init(); % Initialize BSmag analysis (second coil)

% Source points (where there is a current source)
theta = linspace(-8*2*pi,8*2*pi,800);
z_min = 5e-3;
z_max = 25e-3;
z = linspace(z_min, z_max, 800);
r_min = 15e-3;
r_max = 55e-3;
r = linspace(r_min, r_max, 8*100);
Gamma = [r'.*cos(theta'),r'.*sin(theta'),z']; % x,y,z [m,m,m]
I = 1; % filament current [A]
dGamma = 1e9; % filament max discretization step [m]
[BSmag] = BSmag_add_filament(BSmag,Gamma,I,dGamma);


%source points for the second coil
z_min1 = -25e-3;
z_max1 = -5e-3;
z1 = linspace(z_max1, z_min1, 800);
r_min = 15e-3;
r_max = 55e-3;
r = linspace(r_min, r_max, 8*100);
Gamma1 = [r'.*cos(theta'),r'.*sin(theta'),z1']; % x,y,z [m,m,m]
% Note that if you change the sign of the cosine OR sine term then the 
% windings will be opposite. THis will change the field between a null
% in the center and a maximum in the center
I = 1; % filament current [A]
dGamma = 1e9; % filament max discretization step [m]
[BSmag1] = BSmag_add_filament(BSmag1,Gamma1,I,dGamma);

% Field points (where we want to calculate the field)
x_M = linspace(-3*r_max,3*r_max,21); % x [m]
y_M = linspace(-3*r_max,3*r_max,22); % y [m]
z_M = linspace(5*z_min1,5*z_max,23); % z [m]
[X_M,Y_M,Z_M]=meshgrid(x_M,y_M,z_M);
BSmag_plot_field_points(BSmag,X_M,Y_M,Z_M); % shows the field points volume

% Biot-Savart Integration
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M);
[BSmag1,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag1,X_M,Y_M,Z_M, BX, BY, BZ);
z_plot_min = min([z_min1, z_min, z_max, z_max1]);
z_plot_max = max([z_min1, z_min, z_max, z_max1]);

% Plot B/|B|
figure(2)
    plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
    hold on
    plot3(Gamma(1,1),Gamma(1,2),Gamma(1,3),'*g') % plot beginning
    plot3(Gamma(end,1),Gamma(end,2),Gamma(end,3),'*m') % plot end
    
    plot3(Gamma1(:,1),Gamma1(:,2),Gamma1(:,3),'.-b')
    plot3(Gamma1(1,1),Gamma1(1,2),Gamma1(1,3),'*g')
    plot3(Gamma1(end,1),Gamma1(end,2),Gamma1(end,3),'*m')
    
figure(3)
    plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
    hold on
    plot3(Gamma1(:,1),Gamma1(:,2),Gamma1(:,3),'.-r')
    normB=sqrt(BX.^2+BY.^2+BZ.^2);
    
    quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'b')
%axis tight

% Plot Bz on the volume
figure(4), hold on, box on, grid on
	plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
    plot3(Gamma1(:,1),Gamma1(:,2),Gamma1(:,3),'.-r')
	slice(X,Y,Z,BZ,[0],[],[z_min,0,z_max]), colorbar % plot Bz
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Bz [T]')
view(3), axis equal, axis tight
caxis([-0.5,0.5]*1e-5)

% Plot some flux tubes
figure(5), hold on, box on, grid on
	plot3(Gamma(:,1),Gamma(:,2),Gamma(:,3),'.-r') % plot filament
    plot3(Gamma1(:,1),Gamma1(:,2),Gamma1(:,3),'.-r')
	[X0,Y0,Z0] = ndgrid(-3*r_max:r_max:3*r_max,-3*r_max:r_max:3*r_max,5*z_plot_min); % define tubes starting point
    [X1,Y1,Z1] = ndgrid(-3*r_max:r_max:3*r_max,-3*r_max:r_max:3*r_max,5*z_plot_max); % define tubes starting point
	htubes = streamtube(stream3(X,Y,Z,BX,BY,BZ,X0,Y0,Z0), [0.2 10]);
    htubes = streamtube(stream3(X,Y,Z,BX,BY,BZ,X1,Y1,Z1), [0.2 10]);
xlabel ('x [m]'), ylabel ('y [m]'), zlabel ('z [m]'), title ('Some flux tubes')
view(3), axis equal, axis tight
set(htubes,'EdgeColor','none','FaceColor','c') % change tube color
camlight left % change tube light