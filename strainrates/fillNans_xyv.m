function [filledvx,filledvy,filledvv] = fillNans_xyv(vx,vy,vv,limits)
% Open a new figure and display the image of interest
figure()
filledvv = vv;
filledvx = vx;
filledvy = vy;
imagesc(filledvv)

% Set the axes to be able to see the right area
if exist('limits','var')==1    
    axis(limits)
end
    
caxis([0 5])
resp = 'y';
while strcmpi(resp,'y')
    % Use the mouse to select a rectange around a NaN region
    rect = getrect();

    % Extract the portion of the grid in the rectangle
    xmin = floor(rect(1));
    xmax = floor(rect(1))+ceil(rect(3));
    ymin = floor(rect(2));
    ymax = floor(rect(2))+ceil(rect(4));
    subGridvv = filledvv(ymin:ymax,xmin:xmax);
    subGridvx = filledvx(ymin:ymax,xmin:xmax);
    subGridvy = filledvy(ymin:ymax,xmin:xmax);

    %Reshape into a vector
    [rows,cols] = size(subGridvv);
    dataVecvv = reshape(subGridvv',rows*cols,1);
    dataVecvx = reshape(subGridvx',rows*cols,1);
    dataVecvy = reshape(subGridvy',rows*cols,1);

    % Make vectors for x and y
    x = 1:cols;
    y = 1:rows;
    [X,Y] = meshgrid(x,y);
    X_vec = reshape(X',rows*cols,1);
    Y_vec = reshape(Y',rows*cols,1);

    % Assemble into single data frame
    Avv = [X_vec, Y_vec, dataVecvv];
    Avx = [X_vec, Y_vec, dataVecvx];
    Avy = [X_vec, Y_vec, dataVecvy];

    % Get rid of rows with NaNs
    Avv(isnan(Avv(:,3))==1,:) = [];
    Avx(isnan(Avx(:,3))==1,:) = [];
    Avy(isnan(Avy(:,3))==1,:) = [];

    % Regrid data, interpolating the missing NaNs
    Bvv = griddata(double(Avv(:,1)),double(Avv(:,2)),double(Avv(:,3)),double(X),double(Y));
    Bvx = griddata(double(Avx(:,1)),double(Avx(:,2)),double(Avx(:,3)),double(X),double(Y));
    Bvy = griddata(double(Avy(:,1)),double(Avy(:,2)),double(Avy(:,3)),double(X),double(Y));
    
    filledvv(ymin:ymax,xmin:xmax) = Bvv;
    filledvx(ymin:ymax,xmin:xmax) = Bvx;
    filledvy(ymin:ymax,xmin:xmax) = Bvy;
    imagesc(filledvv)
    if exist('limits','var')==1   
        axis(limits)
    end
    caxis([0 5])
    resp = input('Do you wish to define another rectangle? Y/N: ','s');
    if strcmpi(resp,'n')
        break
    end
end

