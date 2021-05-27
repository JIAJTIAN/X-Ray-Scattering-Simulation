%% 3D-Image Simulation 4.0
% 1. Generate Spheres
% 2. FFT
% 3. Self-avoiding

clc,clear
close all

%% Setup
SaveFolder = 'E:\文档\Beamline\Simulation\Result\3DSimulation3\';

% Creat folder to store results
if ~exist(SaveFolder, 'dir') % Check if the folder is not exist
    mkdir(SaveFolder)        % If not, create one
end

Resolution = 1e2;   % Modeling resolution
NumPtcl = 10;       % Particle number
RadiusParticle = 15;% Particle radius
PtclIdx = 1;

%% Generate Real Space 3D-Image
SurfaceZoom = RadiusParticle; % A surface layer to avoid particle out of boundary

ParticleVoxels = zeros(Resolution,Resolution,Resolution); % Zero Matrix for storing binary data

figure;
while PtclIdx <= NumPtcl % Loop for all particles
        
    XC(PtclIdx) = randi([SurfaceZoom , (Resolution-SurfaceZoom)],1,1); % Generate particle center coordinates
    YC(PtclIdx) = randi([SurfaceZoom , (Resolution-SurfaceZoom)],1,1);
    ZC(PtclIdx) = randi([SurfaceZoom , (Resolution-SurfaceZoom)],1,1);
    
    [XInImage, YInImage, ZInImage] = meshgrid(1:Resolution,1:Resolution, 1:Resolution); % Generate meshgrid coordinate of the box
    SingleVoxels = (XInImage - XC(PtclIdx)).^2/RadiusParticle.^2 + ...
        (YInImage - YC(PtclIdx)).^2/RadiusParticle.^2 + (ZInImage - ZC(PtclIdx)).^2/RadiusParticle.^2 <= 1; % Calculate the voxels inside a sphere
    
   
    VoxelsTemp = ParticleVoxels + SingleVoxels;
    
    if sum(VoxelsTemp(VoxelsTemp>1)) == 0  % Self avoiding check: If overlap, the voxel > 1 
        ParticleVoxels = VoxelsTemp;
        
        
        % Generate spheres in a surface figure for ploting
        [x,y,z] = GenerateSphere(XC(PtclIdx),YC(PtclIdx),ZC(PtclIdx),RadiusParticle); 
        
        % Plot sphere
        SurfFig = surf(x,y,z);
        hold on
        set(SurfFig,'FaceLighting','phong','FaceColor','blue',...
            'EdgeColor','none','AmbientStrength',0.5)
        
        PtclIdx = PtclIdx + 1;
    end
    
end
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([0,Resolution])
ylim([0,Resolution])
zlim([0,Resolution])

light('Position',[1 0 0],'Style','infinite'); % Add light effect to 3d plot

FigureTitle = ['3DView_N=',num2str(NumPtcl)] ;
title(FigureTitle,'interpreter','none');
SaveTitle = [FigureTitle,'.png'];
SavePath = strcat(SaveFolder,SaveTitle);
saveas(gcf,SavePath);

%% Projection on Z-Axis
ProjVoxels = squeeze(sum(ParticleVoxels,3)); % Projection along Z-axis(3)

% Plot projection image
figure;
imagesc(ProjVoxels);
xlim([0,Resolution])
ylim([0,Resolution])
colormap cool

FigureTitle = ['ProjectionImage_N=',num2str(NumPtcl)] ;
title(FigureTitle,'interpreter','none');
SaveTitle = [FigureTitle,'.png'];
SavePath = strcat(SaveFolder,SaveTitle);
saveas(gcf,SavePath);

%% 2D-FFT

% Setup a range for intnsity image
MaxInt=log10(10.^5*NumPtcl); 
MinInt=log10(10.^1*NumPtcl);

% Fast Fourier transform
IntImage = abs(fftshift(fft2(ProjVoxels))).^2;

% Plot FFT image
figure;
colormap jet
xIntLim = linspace(-pi,pi,Resolution);
yIntLim = linspace(-pi,pi,Resolution);
imagesc(xIntLim,yIntLim,log10(IntImage),[MinInt,MaxInt]);

FigureTitle = ['ScatteringImage_N=',num2str(NumPtcl)];
title(FigureTitle,'interpreter','none');
SaveTitle = [FigureTitle,'.png'];
SavePath = strcat(SaveFolder,SaveTitle);
saveas(gcf,SavePath);

figure;
[vol_handle]=VoxelPlotter(ParticleVoxels,1); 
%visual effects (I recommend using the FigureRotator function from MATLAB
%Centeral
view(3);
daspect([1,1,1]);
set(gca,'xlim',[0 Resolution], 'ylim',[0 Resolution], 'zlim',[0 Resolution]);