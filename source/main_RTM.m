%% Seismic Reverse-time Migration (RTM) Example
clear;
close all;

method = 'adj_conv';
%method = 'recon_corr';
diagHess = 1;
bg_velocity = 'smooth';
%bg_velocity = 'true';

%% Read in velocity model and plot it

% two-layer model------
velocityModel=zeros(100,100);
velocityModel(:,:)=3.0;
velocityModel(51,:)=2.0; % 4.0
%----------------------

% fault model----------
% tmp=load('model.txt');
% nx1=600;
% nz1=300;
% tmp=reshape(tmp,nz1,nx1);
% velocityModel=tmp(51:251,400:600);
%----------------------

% salt-dome model
% tmp=load('model.txt');
% nx1=600;
% nz1=300;
% velocityModel=reshape(tmp,nz1,nx1);

[nz,nx]=size(velocityModel);

dx = 5;
dz = 5;
x = (1:nx)*dx;
z = (1:nz)*dz;

subplot(2,2,1);
imagesc(x,z,velocityModel);
xlabel('Distance (m)'); ylabel('Depth (m)');
title('Velocity Model');
hold on;
hshot = plot(x(1),z(1),'w*','LineWidth',3,'MarkerSize',12);
hold off;
colorbar;
colormap('jet');

%% Create shot gathers
% Use the velocity model to simulate wavefield. The 2nd-order acoustic wave equation
% is solved using finite differences for a defined initial wavefield.

% calculate time step dt from stability crierion for finite difference
% solution of the wave equation.
dt = 0.9*min(min(dz./velocityModel/sqrt(2)));

% determine time samples nt from wave travelime to depth and back to surface
vmin = min(velocityModel(:));
nt = round(sqrt((dx*nx)^2 + (dz*nx)^2)*2/vmin/(2*dt) + 1);
t  = (0:nt-1).*dt;

% add region around model for applying absorbing boundary conditions (20 nodes wide)
V = [repmat(velocityModel(:,1),1,20) velocityModel repmat(velocityModel(:,end),1,20)];
V(end+1:end+20,:) = repmat(V(end,:),20,1);

% Define frequency parameter for Ricker source wavelet
f  = 0.05;

%% Generate shots and save to file and video

data = zeros(size(nt,nx));
figure(gcf)

for ixs = 21:nx+20  % shot loop
    
    % initial wavefield
    rw = ricker(f,nz+40,dt,dt*ixs,0);
    rw = rw(1:nz+20,:);
    
    % plot initial wavefield
    set(hshot,'XData',x(ixs-20),'YData',z(1));
    subplot(2,2,2);
    imagesc(x,z,rw(1:end-20,21:end-20));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(['Shot ',num2str(ixs-20),' at ',num2str(x(ixs-20)),' m']);
    colormap('jet');
    
    % generate shot record
    tic
    [data, snapshot] = fm2d(V,rw,dx,nt,dt);
    toc
    save(['./snapshot/snapshot',num2str(ixs-20),'.mat'],'snapshot');
    save(['./data/shotfdm',num2str(ixs-20),'.mat'],'data');
    
    data = data(21:end-20,:)';
    
    if ismember(ixs-20,[1 nx/2 nx])
        start = 1;
    else
        start = nt;
    end
    
    for i = start:nt
        % plot shot record evolution
        ds = zeros(nt,nx);
        ds(1:i,:) = data(1:i,:);
        subplot(2,2,3);
        imagesc(x,t,ds);
        xlabel('Distance (m)'); ylabel('Time (s)');
        title('Shot Record');
        %caxis([-0.1 0.1])
        
        % plot wave propagation
        subplot(2,2,4);
        imagesc(x,z,snapshot(1:end-20,21:end-20,i));
        xlabel('Distance (m)'); ylabel('Depth (m)');
        title(['Wave Propagation t = ',num2str(t(i),'%10.3f')]);
        %caxis([-0.14 1])
        
        %writeVideo(vidObj,getframe(gcf));
        drawnow;
    end %shot loop
end


%% Process Shots - Reverse Time Migration

Stacked = zeros(nz+20,nx+40);
colormap('jet');

if ~strcmp(bg_velocity,'true')
    V = zeros(size(V))+V(1); % two-layer model
    %V = imgaussfilt(V,70); % fault model
end

for ixs = 1:nx  % shot loop
    
    load(['./data/shotfdm',num2str(ixs),'.mat']);
    shot = data(21:end-20,:)';
    
    tic
    if strcmp(method,'recon_corr')
        [~, rtmsnapshot] = reconstruct2d(V,data,dx,dt);      % reconstructed wavefield
    elseif strcmp(method,'adj_conv')
        [~, rtmsnapshot] = adjoint2d(V,fliplr(data),dx,dt);  % adjoint wavefield
    end
    toc
    %save(['faultModelData\rtmsnapshot',num2str(ixs),'.mat'],'rtmsnapshot');
    
    load(['./snapshot/snapshot',num2str(ixs),'.mat']);
    
    M = 0;
    s2 = zeros(size(V)); %s2 = 0;
    for i = 1:nt
        if strcmp(method,'adj_conv')
            M = snapshot(:,:,i) .* rtmsnapshot(:,:,nt-i+1) + M;  % conv imaging condition (adjoint wavefield)
        elseif strcmp(method,'recon_corr')
            M = snapshot(:,:,i) .* rtmsnapshot(:,:,i) + M;       % xcorr imaging condition (reconstructed wavefield)
        end
        
        if diagHess==1
            s2 = snapshot(:,:,i).^2 + s2; % on: approx diagonal Hessian
        else
            s2(:,:)=1.0;                  % off: approx diagonal Hessian
        end
        
        if ismember(ixs,[1 nx/2 nx])
            subplot(2,2,3);
            imagesc(x,z,snapshot(1:end-20,21:end-20,i));
            xlabel('Distance (m)'); ylabel('Depth (m)');
            title(['Forward Time Wave Propagation t = ',num2str(t(i),'%10.3f')]);
            %caxis([-0.14 1])
            
            subplot(2,2,4);
            imagesc(x,z,rtmsnapshot(1:end-20,21:end-20,nt-i+1));
            xlabel('Distance (m)'); ylabel('Depth (m)');
            title('Reverse Time Wave Propagation');
            %caxis([-0.14 1])
            
            subplot(2,2,2);
            imagesc(x,z,diff(M(1:end-20,21:end-20)./s2(1:end-20,21:end-20),2,1));
            xlabel('Distance (m)'); ylabel('Depth (m)');
            title(['Current Migrated Shot ',num2str(ixs)]);
            %caxis([-.05 .05])
            
            drawnow
            %writeVideo(vidObj,getframe(gcf));
        end 
    end
    
    if diagHess==1
        Stacked = Stacked + M./(s2+1e-5);  % on: illumination
    else
        Stacked = Stacked + M;             % off: illumination
    end
    
    subplot(2,2,2);
    imagesc(x,z,diff(Stacked(1:end-20,21:end-20),2,1));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title('Stacked Image');
    if diagHess==0
        clim([-200 200]);
    else
        clim([-30 30]);
    end

    subplot(2,2,3);
    imagesc(x,t,shot);
    xlabel('Distance (m)'); ylabel('Time (s)');
    title(['Current Shot ',num2str(ixs)]);
    clim([-0.1 0.1]);
    
    subplot(2,2,4);
    imagesc(x,z,diff(M(1:end-20,21:end-20)./s2(1:end-20,21:end-20),2,1));
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(['Current Migrated Shot ',num2str(ixs)]);
    %clim([-10 10]);
    
    set(hshot,'XData',x(ixs));
    drawnow
    %writeVideo(vidObj,getframe(gcf));
end
