%--------------------------------------------------------------------------
% Image shifting - Bi-linear and Windowed-Sinc
%--------------------------------------------------------------------------

close all;
clear all;

%--------------------------------------------------------------------------
% Editable variables
%--------------------------------------------------------------------------
% Change this to use different images
fileName = 'hawaii.tif';
% Change these values for different shifts
sig1 = 0;    % Horizontal shift
sig2 = -0.5;    % Vertical shift
% Change this value for different window sizes
N = 10;

%--------------------------------------------------------------------------
% Import image to Matlab
%--------------------------------------------------------------------------
imFile = imread(fileName);  % Reads the file into Matlab
[imRow,imCol,imComp] = size(imFile);    % Get dimensions of the image

%--------------------------------------------------------------------------
%   Bi-linear Shifting
%--------------------------------------------------------------------------
imProcBi = zeros(imRow+2, imCol+2, imComp);   % Create matrix with boundary for processing
[imBiRow,imBiCol,imBiComp] = size(imProcBi);   % Get dimensions of processing matrix
% Copy image matrix into processing matrix with boundary of 2
for i=1:imBiComp
    imProcBi(2:imRow+1, 2:imCol+1, i) = imFile(1:imRow,1:imCol,i);  % Add original image file into processing matrix
end

% Symmetric boundary extension for Bilinear interpolation
for i=1:imBiComp  % Set for each color component
    % Top extension
    imProcBi(1,1:imBiCol,i) = imProcBi(3,1:imBiCol,i);
    % Bottom extension
    imProcBi(imBiRow,1:imBiCol,i) = imProcBi(imBiRow-2,1:imBiCol,i);
    % Left extension
    imProcBi(1:imBiRow,1,i) = imProcBi(1:imBiRow,3,i);
    % Right extension
    imProcBi(1:imBiRow,imBiCol,i) = imProcBi(1:imBiRow,imBiCol-2,i);
end

% Bi-linear shifting operations - Horizontal
for i=1:imBiComp
   for k=2:imBiCol
       for j=2:imBiRow
           if sig1 > 0
                tempval = (1-sig1)*imProcBi(j+1,k,i) + sig1*imProcBi(j,k,i);
           elseif sig1 < 0
               tempval = (1-abs(sig1))*imProcBi(j-1,k,i) + abs(sig1)*imProcBi(j,k,i);
           elseif sig1 == 0
               tempval = imProcBi(j,k,i);
           end
           % Overflow check
           if tempval > 255
               tempval = 255;
           elseif tempval < 0
               tempval = 0;
           end
           imProcBi(j,k,i) = tempval;
       end
   end
end

% Bi-linear shifting - Vertical
for i=1:imBiComp
   for k=2:imBiRow
       for j=2:imBiCol
           if sig2 > 0
               tempval = (1-sig2)*imProcBi(k,j+1,i) + sig2*imProcBi(k,j,i);
           elseif sig2 < 0
               tempval = (1-abs(sig2))*imProcBi(k,j-1,i) + abs(sig2)*imProcBi(k,j,i);
           elseif sig2 == 0
               tempval = imProcBi(k,j,i);
           end
           % Overflow check
           if tempval > 255
               tempval = 255;
           elseif tempval < 0
               tempval = 0;
           end
           imProcBi(k,j,i) = tempval;
       end
   end
end

% Remove border from the processed image
imProcBi = imProcBi(2:imBiRow-1, 2:imBiCol-1,:);

%--------------------------------------------------------------------------
%   Windowed-Sinc Shifting
%--------------------------------------------------------------------------
imProcSinc = zeros(imRow+2*N, imCol+2*N, imComp);     % Create an image matrix for differing window sizes
[imSincRow, imSincCol, imSincComp] = size(imProcSinc);
for i=1:imComp
    imProcSinc(N+1:imRow+N, N+1:imCol+N, i) = imFile(1:imRow, 1:imCol, i);  % Add original image file
end

% Symmetric boundary extension
for i=1:imSincComp
    for j=0:N-1
        % Top extension
        imProcSinc(N-j,1:imSincCol,i) = imProcSinc(N+j+1,1:imSincCol,i);
        % Bottom extension
        imProcSinc(imSincRow-N+j+1,1:imSincCol,i) = imProcSinc(imSincRow-N-j,1:imSincCol,i);
        % Left extension
        imProcSinc(1:imSincRow,N-j,i) = imProcSinc(1:imSincRow,N+j+1,i);
        % Right extension
        imProcSinc(1:imSincRow,imSincCol-N+j+1,i) = imProcSinc(1:imSincRow,imSincCol-N-j,i);
    end
end

% Create windowed-sinc functions
t = -N:N;   % 2N+1 window taps
tau1 = N+sig1;  % Select tau value for the Hanning window
tau2 = N+sig2;  % Select tau for vertical Hanning window
% Horizontal window
sincHor = sin(pi*(t-sig1))./(pi*(t-sig1));   % Sinc function, must use ./ for array division
windHor = 0.5*(1 + cos(pi*(t-sig1)/tau1));   % Hanning window, must use .* for array multiplication
% Check for NaN in sincHor, meaning sig1 = 0 
if sig1 == 0
    index = find(isnan(sincHor));
    sincHor(index) = 1;
end
windSincHor = windHor.*sincHor;     % Windowed sinc
windSincHor = 1/sum(windSincHor) * windSincHor; % DC gain = 1
% Vertical window
sincVer = sin(pi*(t-sig2))./(pi*(t-sig2));   % Sinc function
windVer = 0.5*(1 + cos(pi*(t-sig2)/tau2));   % Hanning window
% Check for NaN in sincHor, meaning sig2 = 0
if sig2 == 0
    index = find(isnan(sincVer));
    sincVer(index) = 1;
end
windSincVer = windVer.*sincVer;     % Windowed sinc
windSincVer = 1/sum(windSincVer) * windSincVer; % DC gain = 1

% Horizontal
for i=1:imSincComp
    for j=1:imSincRow
        imProcSinc(j,:,i) = filter(windSincHor,1,imProcSinc(j,:,i));
    end
    for k=1:imSincCol
        imProcSinc(:,k,i) = filter(windSincVer,1,imProcSinc(:,k,i));
    end
end

% Remove border from processed image
imProcSinc = imProcSinc(2*N+1:imSincRow, 2*N+1:imSincCol,:);

% Plot images for comparison
figure('Name','Original Image');
imshow(uint8(imFile));
figure('Name','Bi-linear Shift');
imshow(uint8(imProcBi));
figure('Name','Windowed-Sinc');
imshow(uint8(imProcSinc));

% % Output images to the same directory
% imwrite(uint8(imProcBi),strcat('BiLin -',fileName));
% imwrite(uint8(imProcSinc),strcat('WinSinc -',fileName));
