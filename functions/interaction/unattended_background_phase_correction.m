function varargout = unattended_background_phase_correction(varargin)
% UNATTENDED_BACKGROUND_PHASE_CORRECTION unattended version for background_phase_correction.m
%      Performs manual background phase correction by fitting an nth-degree
%      polynomial to background phase in static tissue and outputs poly fits.
%
%      This script aims to generate an unattended version of the
%      correction, avoiding GUI input. This means that the inputs that are
%      manually selected in the previous script, now are passed as
%      arguments to the function.
%
%      BACKGROUND_PHASE_CORRECTION, by itself, creates a new BACKGROUND_PHASE_CORRECTION or raises the existing
%      singleton*.
%
%      H = BACKGROUND_PHASE_CORRECTION returns the handle to a new BACKGROUND_PHASE_CORRECTION or the handle to
%      the existing singleton*.
%
%      BACKGROUND_PHASE_CORRECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BACKGROUND_PHASE_CORRECTION.M with the given input arguments.
%
%      BACKGROUND_PHASE_CORRECTION('Property','Value',...) creates a new BACKGROUND_PHASE_CORRECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before background_phase_correction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to background_phase_correction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help background_phase_correction

% Last Modified by GUIDE v2.5 10-Sep-2014 15:46:51

%       Kevin Johnson, UW-Madison 2014
%       Used by: load_pcvipr.m
%       Dependencies: evaluate_poly.m


MAG= varargin{1};
VX = varargin{2};
VY = varargin{3};
VZ = varargin{4};

fit_order = varargin{5};
cd_thresh = varargin{6};
noise_thresh = varargin{7};

% polyfit_3D(MAG,VX,VY,VZ,fit_order,cd_thresh, noise_thresh)
% Creates a 3D polynomial fit over static tissue (thresholded areas)

Mask = int8(zeros(size(MAG)));
max_MAG = max(MAG(:));

% Create angiogram
for slice = 1:size(Mask,3)
    mag_slice= single(MAG(:,:,slice));
    vx_slice = single(VX(:,:,slice));
    vy_slice = single(VY(:,:,slice));
    vz_slice = single(VZ(:,:,slice));

    speed = sqrt(vx_slice.^2 + vy_slice.^2 + + vz_slice.^2);
    % Create magnitude-masked CD slice
    Mask(:,:,slice) = ( mag_slice > noise_thresh*max_MAG) .* ( speed < cd_thresh*1000);
end

%lots of memory problems w/ vectorization! solve Ax = B by (A^hA)x = (A^h*b)
if fit_order == 0
    poly_fitx.vals  = mean(VX(:));
    poly_fitx.px = 0;
    poly_fitx.py = 0;
    poly_fitx.pz = 0;
    
    poly_fity.vals  = mean(VY(:));
    poly_fity.px = 0;
    poly_fity.py = 0;
    poly_fity.pz = 0;
    
    poly_fitz.vals  = mean(VZ(:));
    poly_fitz.px = 0;
    poly_fitz.py = 0;
    poly_fitz.pz = 0;
else   
    [px,py,pz] = meshgrid(0:fit_order,0:fit_order,0:fit_order);
    idx2 = find( (px+py+pz) <= fit_order);
    px = px(idx2); % all linear combination of powers in x,y,z
    py = py(idx2);
    pz = pz(idx2);
    A = [px(:) py(:) pz(:)]; %form polynomial matrix operator
    
    N = size(A,1);
    
    AhA = zeros(N,N); %inverse matrix times matrix operator 
    AhBx = zeros(N,1);
    AhBy = zeros(N,1); 
    AhBz = zeros(N,1);
        
    xrange = single( linspace(-1,1,size(VX,1)) );
    yrange = single( linspace(-1,1,size(VY,2)) );
    zrange = single( linspace(-1,1,size(VZ,3)) );
    
    
    for slice = 1:numel(zrange)        
        vx_slice = single(VX(:,:,slice) ); %get velocity slices
        vy_slice = single(VY(:,:,slice) );
        vz_slice = single(VZ(:,:,slice) );
        
        [y,x,z] = meshgrid( yrange,xrange,zrange(slice) );  
        idx = find( Mask(:,:,slice) > 0);
        x = x(idx); %linear meshgrid. x,y,z are acting like variables
        y = y(idx);
        z = z(idx);
        vx_slice = vx_slice(idx);
        vy_slice = vy_slice(idx);
        vz_slice = vz_slice(idx);
        
        for i = 1:N
            for j = 1:N
                AhA(i,j) =  AhA(i,j) + sum( ( x.^px(i).*y.^py(i).*z.^pz(i)).*( x.^px(j).*y.^py(j).*z.^pz(j)) );
            end
        end
        
        for i = 1:N
            AhBx(i) = AhBx(i) + sum( vx_slice.* ( x.^px(i).*y.^py(i).*z.^pz(i)));
            AhBy(i) = AhBy(i) + sum( vy_slice.* ( x.^px(i).*y.^py(i).*z.^pz(i)));
            AhBz(i) = AhBz(i) + sum( vz_slice.* ( x.^px(i).*y.^py(i).*z.^pz(i)));
        end
    end
        
    poly_fitx.vals = linsolve(AhA,AhBx); %solve Ax = B by (A^hA)x = (A^h*b)
    poly_fitx.px = px;
    poly_fitx.py = py;
    poly_fitx.pz = pz;
        
    poly_fity.vals = linsolve(AhA,AhBy);
    poly_fity.px = px;
    poly_fity.py = py;
    poly_fity.pz = pz;
    
    poly_fitz.vals = linsolve(AhA,AhBz);
    poly_fitz.px = px;
    poly_fitz.py = py;
    poly_fitz.pz = pz;
end

varargout{1} = poly_fitx; %array of polynomial fits 
varargout{2} = poly_fity;
varargout{3} = poly_fitz;