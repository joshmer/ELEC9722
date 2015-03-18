%--------------------------------------------------------------------------
% Image Skewing - Bi-linear and Windowed-Sinc
%--------------------------------------------------------------------------

close all;
clear all;

%--------------------------------------------------------------------------
% Editable variables
%--------------------------------------------------------------------------
% Change this to use different images
fileName = 'hawaii.tif';
% Change these values for different shifts
sig1 = 0.3;    % Horizontal skew
sig2 = 0;    % Vertical skew
% Change this value for different window sizes
N = 5;
% Color offset
colorSet = 128;     % Set for desired color border, 0 = black, 128 = gray, 255=white

%--------------------------------------------------------------------------
% Import image and calculate skews
%--------------------------------------------------------------------------
% Import file into Matlab and get image dimensions
imFile = imread(fileName);
[imRow,imCol,imComp] = size(imFile);

% Calculate maximum skewing
maxSkewHor = ceil(abs(sig1*imCol));
borderHor = maxSkewHor;     % Extra border spacing
maxSkewVer = ceil(abs(sig2*imRow));
borderVer = maxSkewVer;

%--------------------------------------------------------------------------
% Bi-linear shifting
%--------------------------------------------------------------------------
imProcBi = colorSet*ones(imRow+2*borderVer, imCol+2*borderHor, imComp);     % Create a matrix for Bilinear interpolation
[imBiRow, imBiCol, imBiComp] = size(imProcBi);   % Get dimensions
% Copy the image matrix with an boundary of 2
for i=1:imBiComp
    imProcBi(borderVer+1:imRow+borderVer, borderHor+1:imCol+borderHor, i) = imFile(1:imRow, 1:imCol, i);
end

% Symmetric boundary extension - Bi-Linear interpolation adds 1 extra column
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

imProcBiTemp = colorSet*ones(imBiRow,imBiCol,imBiComp); %Temporary array

% Horizontal skew
if sig1 ~= 0    % Do not calculate if there is no skew
    for i=1:imBiComp
        for j=borderVer+1:imBiRow-borderVer
            % Calculate required shift/new column value
            tempSkew = sig1*(j-borderVer);
            tempCol = floor(tempSkew);
            tempSig = tempSkew-tempCol;
            for k=borderHor+1:imBiCol-borderHor
                % Bi-linear interpolation with new, calculated sigma
                tempval = (1-abs(tempSig))*imProcBi(j,k,i)+abs(tempSig)*imProcBi(j,k-1,i);
                imProcBi(j,k-1,i) = colorSet;    % Resets the previous value 
                if tempval > 255
                    tempval = 255;
                elseif tempval < 0
                    tempval = 0;
                end
                imProcBiTemp(j,k+tempCol-sign(sig1)*floor(borderHor/2),i) = tempval;
            end
            imProcBi(j,k,i) = colorSet; % Zero out the final column
        end
    end
    imProcBi = imProcBiTemp;
end

figure('Name','Bi-linear Skewing');
subplot(1,2,1);
imshow(uint8(imProcBi));

% Vertical skew
if sig2 ~= 0    % If there is no vertical skew, do not perform
    imProcBiTemp = colorSet*ones(imBiRow,imBiCol,imBiComp); %Reset array
    for i=1:imBiComp
        for k=1:imBiCol
            % Calculate required shift/new row value
            tempSkew = sig2*(k-borderHor);
            tempRow = floor(tempSkew);
            tempSig = tempSkew - tempRow;
            for j=borderVer+1:imBiRow-borderVer
                % Bi-linear interpolation with new, calculated sigma
                tempval = (1-abs(tempSig))*imProcBi(j,k,i) + abs(tempSig)*imProcBi(j-1,k,i);
                if tempval > 255
                    tempval = 255;
                elseif tempval < 0
                    tempval = 0;
                end
                imProcBiTemp(j+tempRow-sign(sig2)*floor(borderVer/2),k,i) = tempval;
            end
        end
    end
    imProcBi = imProcBiTemp;
end

subplot(1,2,2);
imshow(uint8(imProcBi));

%--------------------------------------------------------------------------
%   Windowed-Sinc Shifting
%--------------------------------------------------------------------------
imProcSinc = colorSet*ones(imRow + 2*N + 2*borderVer, imCol + 2*N + 2*borderHor, imComp);     % Create an image matrix for differing window sizes
[imSincRow, imSincCol, imSincComp] = size(imProcSinc);
for i=1:imComp
    imProcSinc(borderVer+N+1:imSincRow - borderVer - N, borderHor+N+1:imSincCol-N-borderHor, i) = imFile(1:imRow, 1:imCol, i);  % Add original image file
end

% Symmetric boundary extension
for i=1:imSincComp
    for j=0:N-1
        if sig1 ~= 0 % Don't extend if no horizontal skew
            % Left extension
            imProcSinc(1:imSincRow,borderHor+N-j,i) = imProcSinc(1:imSincRow,borderHor+N+j+1,i);
            % Right extension
            imProcSinc(1:imSincRow,imSincCol-borderHor-N+j+1,i) = imProcSinc(1:imSincRow,imSincCol-borderHor-N-j,i);
        elseif sig2 ~=0 % Don't extend if no vertical skew
            % Top extension
            imProcSinc(borderVer+N-j,1:imSincCol,i) = imProcSinc(borderVer+N+j+1,1:imSincCol,i);
            % Bottom extension
            imProcSinc(imSincRow-borderVer-N+j+1,1:imSincCol,i) = imProcSinc(imSincRow-borderVer-N-j,1:imSincCol,i);
        end
    end
end

% Horizontal interpolation
imProcSincTemp = colorSet*ones(imSincRow,imSincCol,imSincComp);  % Create border for the final image
t = -N:N;   % 2N+1 range, never changes
for i=1:imSincComp
    for j=borderVer+1:imSincRow-borderVer
        % Calculate required shift/new column value
        tempSkew = sig1*(j-borderVer);
        tempCol = floor(tempSkew);
        tempSig = tempSkew-tempCol;

        tau = N+tempSig;
        sincHor = sin(pi*(t-tempSig))./(pi*(t-tempSig));
        % If tempSig = 0, then will have sin(0)/0 = NaN
        if tempSig == 0
            index = find(isnan(sincHor));
            sincHor(index) = 1;     % sinc(0) = 1 normally
        end
        windHor = 0.5*(1 + cos(pi*(t-tempSig)/tau));
        windSincHor = windHor.*sincHor;
        windSincHor = 1/sum(windSincHor) * windSincHor;  % DC gain of 1
        
        % Clip the image
        tempval = filter(windSincHor,1,imProcSinc(j,1:imSincCol,i));   % 1D filter the entire row
        [row,col] = size(tempval);     % Get row dimensions
        % Keep vertical boundaries for next filtering operation
        imProcSincTemp(j,borderHor+N+1+tempCol-sign(sig1)*floor(borderHor/2):imSincCol-borderHor-N-sign(sig1)*floor(borderHor/2)+tempCol,i) = tempval(borderHor+2*N+1:col-borderHor);
    end
end

figure('Name','Windowed Sinc');
subplot(1,2,1);
imshow(uint8(imProcSincTemp));

% Vertical interpolation
imProcSinc = colorSet*ones(imSincRow,imSincCol,imSincComp);  % Reset for final image - removes stray artifacts
for i=1:imSincComp
    for k=1:imSincCol
        % Calculate required shift/new column value
        tempSkew = sig2*(k-borderHor);
        tempCol = floor(tempSkew);
        tempSig = tempSkew-tempCol;

        tau = N+tempSig;
        sincVer = sin(pi*(t-tempSig))./(pi*(t-tempSig));
        % If tempSig = 0, then will have sin(0)/0 = NaN
        if tempSig == 0
            index = find(isnan(sincVer));
            sincVer(index) = 1;     % sinc(0) = 1 normally
        end
        windVer = 0.5*(1 + cos(pi*(t-tempSig)/tau));
        windSincVer = windHor.*sincVer;
        windSincVer = 1/sum(windSincVer) * windSincVer;  % DC gain of 1
        
        % Clip the image
        tempval = filter(windSincVer,1,imProcSincTemp(1:imSincRow,k,i));   % 1D filter the entire row
        [row,col] = size(tempval);     % Get row dimensions
        % Keep just the image portion
        imProcSinc(borderVer+N+1+tempCol-sign(sig2)*floor(borderVer/2):imSincRow-borderVer-N-sign(sig2)*floor(borderVer/2)+tempCol,k,i) = tempval(borderVer+2*N+1:row-borderVer);
    end
end

subplot(1,2,2);
imshow(uint8(imProcSinc));

% % Output the images to the same directory
% imwrite(uint8(imProcBiTemp),strcat('BiLin_Hor-',fileName));
% imwrite(uint8(imProcBi),strcat('BiLin_HorVer-',fileName));
% imwrite(uint8(imProcSincTemp),strcat('WinSinc_Hor-',fileName));
% imwrite(uint8(imProcSinc),strcat('WinSinc_HorVer-',fileName));
