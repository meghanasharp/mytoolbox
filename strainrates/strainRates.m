function [elon,etrans,eshear,exGrid,eyGrid,exyGrid] = strainRates(vx,vy,pixel_size,varargin)
% strainRates: Computes flow-oriented components of the strain-rate tensor
% using the x- and y-components of a raster velocity field.

% Outputs: Longitudinal, transverse, and shear-strain grids

% Inputs:
    % vx: Grid of the x-component of velocity
    % vy: Grid of the y-component of velocity
    % pixel_size: Pixel size in map units

% Name-value options:
    % thick: Thickness grid, if using one.
    % length_scale: Actually half the length-scale over which strain rates
        % will be calculated. Defined in map units. Default is 3000.
    % maxR: Defines how far the "virtual stakes" are allowed to move before 
        % the script continues. If using a thickness grid, choose a 
        % reasonable estimate. Smaller maximum values speed up the script 
        % but cut off the length scale in areas with very thick ice. A 
        % maximum value is necessary for determining where the for-loop 
        % should start in the grid. Default is length_scale/pixel_size if
        % thick_grid_toggle = 0, and 8 if thick_grid_toggle = 1.
    % thick_grid_toggle: Set to 1 if using a thickness grid. Default is 0.
    % thick_multiplier: Set if using a thickness grid. This is used to set
        % the length scale during each iteration. Half-length-scale is set
        % during each iteration as thick_multiplier*thickness
    % tol: Tolerance for error in adaptive time-stepping scheme; value
        % is the percent difference between the two stake position estimates
        % divided by 100
    % vxNoData: Value representing no data in vx grid. Default is NaN.
    % vyNoData: Value representing no data in vy grid. Default is NaN.
    % thickNoData: Value representing no data in thickness grid. Default is NaN.
    % ydir: Defines which y-direction is positive. Set to 1 if positive
        % y-velocities go up, -1 if down. Default is 1.
    
p = inputParser;
addParameter(p,'length_scale',3000)
addParameter(p,'maxR',8)
addParameter(p,'thick_grid_toggle',0)
addParameter(p,'thick_multiplier',8)
addParameter(p,'tol',10^-4)
addParameter(p,'vxNoData',NaN)
addParameter(p,'vyNoData',NaN)
addParameter(p,'thickNoData',NaN)
addParameter(p,'ydir',1)

parse(p,varargin{:})

length_scale = p.Results.length_scale;
maxR = p.Results.maxR;
thick_grid_toggle = p.Results.thick_grid_toggle;
thick_multiplier = p.Results.thick_multiplier;
tol = p.Results.tol;
vxNoData = p.Results.vxNoData;
vyNoData = p.Results.vyNoData;
thickNoData = p.Results.thickNoData;
ydir = p.Results.ydir;

if thick_grid_toggle == 0
    maxR = ceil(length_scale/pixel_size);
end
    
%% Finish setting parameters

time_max = 0.1*maxR*pixel_size/0.01;

r = round(length_scale/pixel_size); % Finds the nearest number of pixels to the given length scale;
if r == 0 % If the length scale rounds to 0, set it to 1
    r = 1;
end
r = cast(r,'double');

% The actual length scale used will be r*pixel_size

% Remove erroneous thickness values
if thick_grid_toggle == 1 
    thick(thick<0)=0;
    thick(thick==thickNoData)=NaN;
end
 
 % Set no data values to NaN
 vx(vx==vxNoData)=NaN;
vy(vy==vyNoData)=NaN;
 vx(vx<-100)=NaN;
vy(vy<-100)=NaN;

% Local square dimensions 
locMult = 2; 
% This sets the dimensions of the local square extracted at each time step;
% the default dimensions are 2*locMult*length_scale

%% Initialize calculations

% Create arrays that the calculated values will be written into
centerAlphas = zeros(size(vx));
curXVels = zeros(1,5);
curYVels = zeros(1,5);
checkXVels = zeros(1,5);
checkYVels = zeros(1,5);
checkRowCoords = zeros(1,5);
checkColCoords = zeros(1,5);
errorCriteriaX = zeros(1,5);
errorCriteriaY = zeros(1,5);
newRowCoords = zeros(1,5);
newColCoords = zeros(1,5);
exGrid = zeros(size(vx)) + NaN;
exyGrid = zeros(size(vx)) + NaN;
eyGrid = zeros(size(vx)) + NaN;
rGrid = zeros(size(vx)) + NaN;

dHdx = zeros(size(vx)) + NaN;
dHdy = zeros(size(vx)) + NaN;

[rows,cols] = size(vx);

% Define the lengths of the segments at the beginning of each calculation
l0a1 = r;
l0a2 = r;
l0b1 = r*sqrt(2);
l0b2 = r*sqrt(2);
l0c1 = r;
l0c2 = r;
l0d1 = r*sqrt(2);
l0d2 = r*sqrt(2);

% Local square dimensions
locDim = (2*locMult*r)+1;
locCent = ceil(locDim/2);

% Assign local coordinates to the stakes around each strain square
rowCoords = [locCent,locCent-r,locCent,locCent+r,locCent];
colCoords = [locCent,locCent,locCent+r,locCent,locCent-r];
% (In order, these are the center point, the top point, the right-hand
% point, the bottom point, and the left-hand point. In other words, it
% starts from the center and then moves clockwise from the top point.)


%% Loop through the velocity grids
for i=locMult*maxR+1:rows-locMult*maxR
    for j=locMult*maxR+1:cols-locMult*maxR
    %% Assign a local length scale if using an ice thickness grid
    if thick_grid_toggle==1
        h = thick(i,j);
        length_scale = thick_multiplier*h;
        r = round(length_scale/pixel_size); % Finds the nearest number of pixels to the given length scale;
        if isnan(h)==1
            continue
        end
        if r == 0 % If the length scale rounds to 0, set it to 1
            r = 1;
        elseif r > maxR
            r = maxR;
        end
        r = cast(r,'double');
        rGrid(i,j) = r;
        
        % Define the lengths of the segments at the beginning of each calculation
        l0a1 = r;
        l0a2 = r;
        l0b1 = r*sqrt(2);
        l0b2 = r*sqrt(2);
        l0c1 = r;
        l0c2 = r;
        l0d1 = r*sqrt(2);
        l0d2 = r*sqrt(2);

        % Local square dimensions
        locDim = (2*locMult*r)+1;
        locCent = ceil(locDim/2);

        % Assign local coordinates to the stakes around each strain square
        rowCoords = [locCent,locCent-r,locCent,locCent+r,locCent];
        colCoords = [locCent,locCent,locCent+r,locCent,locCent-r];
    end
    %% Initialize calculations for each center pixel
        % Set the initial strain experienced by each strain segment to zero
        stota1 = 0;
        stota2 = 0;
        stotb1 = 0;
        stotb2 = 0;
        stotc1 = 0;
        stotc2 = 0;
        stotd1 = 0;
        stotd2 = 0;
        
        % Set the current length of each strain segment to the original
        % lengths comprising the strain square
        lLasta1 = l0a1;
        lLasta2 = l0a2;
        lLastb1 = l0b1;
        lLastb2 = l0b2;
        lLastc1 = l0c1;
        lLastc2 = l0c2;
        lLastd1 = l0d1;
        lLastd2 = l0d2;
        
        % Set the current rows and columns to the coordinates of the strain
        % square
        curRows = rowCoords;
        curCols = colCoords;
        
        % Extract an array around the center point (i,j) that represents
        % twice the dimensions of the strain square in order to make later
        % calculations. This is the "local grid"
        sqVx = vx((i-(locMult*r)):(i+(locMult*r)),(j-(locMult*r)):(j+(locMult*r)));
        sqVy = vy((i-(locMult*r)):(i+(locMult*r)),(j-(locMult*r)):(j+(locMult*r)));
        if thick_grid_toggle == 1
            sqThick = thick((i-(locMult*r)):(i+(locMult*r)),(j-(locMult*r)):(j+(locMult*r)));
        end
        % Extract an array around the center point (i,j) that represents
        % just the strain square in order to calculate the average velocity
        % at the center point and determine a reasonable time interval
        sqVxmean = vx((i-r):(i+r),(j-r):(j+r));
        sqVymean = vy((i-r):(i+r),(j-r):(j+r));
        [sqRows,sqCols] = size(sqVx);
        % Calculate mean velocity
        meanX = mean(mean(sqVxmean, "omitnan"), "omitnan");
        meanY = mean(mean(sqVymean, "omitnan"), "omitnan");
        meanVel = sqrt(meanX^2+meanY^2);
        
        % Let the stakes move by approximately one tenth of the length scale
        time = 0.1*r*pixel_size/meanVel;
        time = min(time, time_max);

        dtOrig = pixel_size/meanVel*.05; % Initialize time step as the time it takes to move 
        % a twentieth of the pixel length, according to the average velocity
        dt = dtOrig;
        t = 0; % Initialize time tracker
        
        % Extract the x- and y-velocities at the original stake positions
        for k = 1:5
            curXVels(k) = sqVx(rowCoords(k),colCoords(k));
            curYVels(k) = sqVy(rowCoords(k),colCoords(k));
        end
        if any(isnan([curXVels,curYVels]))==1
            continue
            % Go to the next pixel in the for-loop if any of the current
            % velocities are NaNs
        end
        
        %% Let the stakes move and measure strain rates
            while t<time
                % Move the stakes according to the x- and y-velocities and the
                % time step. Calculate the new column and row positions.
                % Also calculate the column and row positions according to
                % the improved Euler method to check for accuracy.
                t = t+dt;
                for k = 1:5
                    newRowCoords(k) = curRows(k) - ydir*curYVels(k)*dt/pixel_size;
                    newColCoords(k) = curCols(k) + curXVels(k)*dt/pixel_size;
                    [checkXVels(k),checkYVels(k)] = locInterp2(curRows(k),curCols(k),sqVx,sqVy);
                    checkRowCoords(k) = curRows(k) - ydir*0.5*(curYVels(k)+checkYVels(k))*dt/pixel_size;
                    checkColCoords(k) = curCols(k) + 0.5*(curXVels(k)+checkXVels(k))*dt/pixel_size;
                    errorCriteriaY(k) = abs((checkRowCoords(k)-newRowCoords(k))/checkRowCoords(k));
                    errorCriteriaX(k) = abs((checkColCoords(k)-newColCoords(k))/checkColCoords(k));
                end
                if any(isnan([checkXVels,checkYVels]))==1
                    break
                end
                %% Adaptive time stepping
                if any([errorCriteriaY,errorCriteriaX] >= tol)==1
                    t = t-dt; % Reverse the time to what it was before
                    dt = dt/2; % Make a smaller time step
                    % We leave the current rows, columns, and velocities the
                    % same    
                else
                %% Make final calculations
                % Check to be sure that the stakes haven't moved outside of the
                % local grid
                    if max(newRowCoords) <= sqRows && max(newColCoords) <= sqCols && min(newRowCoords)>=1 && min(newColCoords)>=1
                        % Calculate the current length 
                        lfa1 = sqrt((newRowCoords(1) - newRowCoords(2))^2+(newColCoords(1)-newColCoords(2))^2);
                        lfa2 = sqrt((newRowCoords(1) - newRowCoords(4))^2+(newColCoords(1)-newColCoords(4))^2);
                        lfb1 = sqrt((newRowCoords(2) - newRowCoords(5))^2+(newColCoords(2)-newColCoords(5))^2);
                        lfb2 = sqrt((newRowCoords(3) - newRowCoords(4))^2+(newColCoords(3)-newColCoords(4))^2);
                        lfc1 = sqrt((newRowCoords(1) - newRowCoords(5))^2+(newColCoords(1)-newColCoords(5))^2);
                        lfc2 = sqrt((newRowCoords(1) - newRowCoords(3))^2+(newColCoords(1)-newColCoords(3))^2);
                        lfd1 = sqrt((newRowCoords(2) - newRowCoords(3))^2+(newColCoords(2)-newColCoords(3))^2);
                        lfd2 = sqrt((newRowCoords(4) - newRowCoords(5))^2+(newColCoords(4)-newColCoords(5))^2);

                        % Calculate the current strains and strain rates
                        stra1 = log(lfa1/lLasta1);
                        stra2 = log(lfa2/lLasta2);
                        strb1 = log(lfb1/lLastb1);
                        strb2 = log(lfb2/lLastb2);
                        strc1 = log(lfc1/lLastc1);
                        strc2 = log(lfc2/lLastc2);
                        strd1 = log(lfd1/lLastd1);
                        strd2 = log(lfd2/lLastd2);

                        ea1 = stra1/dt;
                        ea2 = stra2/dt;
                        eb1 = strb1/dt;
                        eb2 = strb2/dt;
                        ec1 = strc1/dt;
                        ec2 = strc2/dt;
                        ed1 = strd1/dt;
                        ed2 = strd2/dt;
                  
                        % Update the new rows and columns as current
                        curRows = newRowCoords;
                        curCols = newColCoords;
                        
                        % Update the current lengths as the previous
                        % lengths
                        lLasta1 = lfa1;
                        lLasta2 = lfa2;
                        lLastb1 = lfb1;
                        lLastb2 = lfb2;
                        lLastc1 = lfc1;
                        lLastc2 = lfc2;
                        lLastd1 = lfd1;
                        lLastd2 = lfd2;
                        
                        % Update the running total of strain for each
                        % segment
                        stota1 = stota1 + stra1;
                        stota2 = stota2 + stra2;
                        stotb1 = stotb1 + strb1;
                        stotb2 = stotb2 + strb2;
                        stotc1 = stotc1 + strc1;
                        stotc2 = stotc2 + strc2;
                        stotd1 = stotd1 + strd1;
                        stotd2 = stotd2 + strd2;
                        
                        % Calculate final strain rate components
                        ea1f = stota1/t;
                        ea2f = stota2/t;
                        eb1f = stotb1/t;
                        eb2f = stotb2/t;
                        ec1f = stotc1/t;
                        ec2f = stotc2/t;
                        ed1f = stotd1/t;
                        ed2f = stotd2/t;

                        % Reset the time step
                        dt = dtOrig;

                        % Set the current velocities to those calculated at
                        % the end of the time step
                        curXVels = checkXVels;
                        curYVels = checkYVels;

             
                    else
                        break
                    end
                    % If a stake has moved outside the strain square, leave the
                    % current rows and columns and current velocities the same
                    % as the previous time step, and simply move on to
                    % calculating the strain rate components. This will be the
                    % final value for this strain square. Change the time to
                    % kick it out of the while loop
                end
            end
            if t<.5*time % Don't calculate values if the stakes have been 
                % allowed to move for less than half of the designated time
                % increment
                exGrid(i,j) = nan;
                eyGrid(i,j) = nan;
                exyGrid(i,j) = nan;
                centerAlphas(i,j) = nan;
            else
            %% Create strain rate grids

            % Average the strain rate components
            ea = (ea1f + ea2f)/2;
            eb = (eb1f + eb2f)/2;
            ec = (ec1f + ec2f)/2;
            ed = (ed1f + ed2f)/2;
            
            % Calculate coordinate-oriented strain values
            exGrid(i,j) = .25*(eb + ed - ea) + .75*ec;
            exyGrid(i,j) = .5*eb - .5*ed;
            eyGrid(i,j) = .75*ea + .25*(eb + ed - ec);
    
            
            %% Calculate flow orientation
            % Calculate a grid of flow directions so that the grid-oriented
            % strain rates can be rotated outside of the for-loop to align with
            % local flow directions
            centerVelX = vx(i,j);
            centerVelY = vy(i,j);
            if centerVelX>0 && centerVelY>0
                centerAlphas(i,j) = atand(centerVelY/centerVelX);
            elseif centerVelX<0 && centerVelY>0
                centerAlphas(i,j) = atand(centerVelY/centerVelX)+180;
            elseif centerVelX<0 && centerVelY<0
                centerAlphas(i,j) = atand(centerVelY/centerVelX)+180;
            elseif centerVelX>0 && centerVelY<0
                centerAlphas(i,j) = atand(centerVelY/centerVelX)+360;
            elseif centerVelX>0 && centerVelY==0
                centerAlphas(i,j) = 0;
            elseif centerVelX==0 && centerVelY>0
                centerAlphas(i,j) = 90;
            elseif centerVelX<0 && centerVelY==0
                centerAlphas(i,j) = 180;
            elseif centerVelX==0 && centerVelY<0
                centerAlphas(i,j) = 270;
            end
            end
    end  
    
end

%% Rotate strains and finalize grids
elon = exGrid.*cosd(centerAlphas).^2 + 2*exyGrid.*sind(centerAlphas).*cosd(centerAlphas) + eyGrid.*sind(centerAlphas).^2;
etrans = exGrid.*sind(centerAlphas).^2 - 2*exyGrid.*sind(centerAlphas).*cosd(centerAlphas)+eyGrid.*cosd(centerAlphas).^2;
eshear = (eyGrid-exGrid).*sind(centerAlphas).*cosd(centerAlphas) + exyGrid.*(cosd(centerAlphas).^2-sind(centerAlphas).^2);

eEff = sqrt(abs(exGrid.*eyGrid-exyGrid.^2));
ez = -elon-etrans;

