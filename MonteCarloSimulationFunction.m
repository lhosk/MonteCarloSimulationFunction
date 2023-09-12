clf; close all; clear all; clc;
%%                               lhosk, Computational Physics: Project 7
%Example of MonCamSimFunction to visualize area between circles or spheres

Radius1 = input('Input Radius 1: ');
Radius2 = input('Input Radius 2: ');
if Radius1 <= 0 || Radius2 <= 0
    error('Radius 1 and 2 must be positive!')
end
HoleArea = input(['If you would like circle/sphere 2 to be a hole, ' ...
                  'type "1". ' ...
                  'If you would like it to occupy area, type "0". ']);
Array = [Radius1 0 0 0 0
         Radius2 0 0 0 HoleArea
               4 3 2 1 1];
SampleSize = input('Input sample size for each loop. (ex. 10000): ');
NumberOfLoops = input(['Input Number of Loops (positive int' ...
                       'eger)(ex. 25): ']);
SampleSize = SampleSize*ones(1,NumberOfLoops);
d = input('Input the max size of d (ex. 2): ');
% ------------------------------------------------------------------------
MonCamSim(Array, 3, SampleSize, NumberOfLoops, d)

%%                                                       MonCamSimFunction

function MonCamSim(Array,Dim,SampleSize,NumberOfLoops,d)
%{
% Monte Carlo Simulation: lhosk


%%                                                            Introduction
% ------------------------------------------------------------------------


% The purpose of the MonCamSim function is to estimate the combined area 
% and total area of circles in 2-dimensional space and spheres in 
% 3-dimensional space. Any circle/sphere can either be a hole or occupy 
% space. The math behind calculating both types of occupied areas use a 
% Monte Carlo Simulation. This causes dots to be randomized between the 
% limits of the figure. If these dots land inside the circles/sphere, they 
% account for the area over the figure size. The radius for n number of 
% circles and spheres can have a length of any size. Also, depending on 
% their dimensions, they can be positioned anywhere in the x, y, and z 
% direction. To better understand this function, figures are shown for 
% each example. Lastly, keep in mind that some figures may appear on top
% of each other.


%%                                                          Current Issues
% ------------------------------------------------------------------------


% 1. When using MLH (The Minimized Limit Holder) 3 Dimensionally with 2 
%    Figures, where part of a hole resides both within and outside the 
%    limits of the figure, there is a visual misinterpreation in which the 
%    hole allows area given from non-hole spheres. Matlab currently does
%    not have a slice/volume fill function to fix this. Although, this
%    only change the visuals and not the estimations behind MonCamSim.
%
% 2. When using 3 Dimensional figures, the xlabel does not show 'X axis'
%    due to the fact that there is a comment showing the statistics of the
%    figure shown. 
%
% 3. For figures in which the area occupied by more than one circle/sphere
%    is a large percent of the total Area of Rectangle/Cube (RAD), it take
%    a increased amount of time to plot each of the plot points onto the
%    figure. To limit the run time of each program will keeping the
%    visuals intact, the local function CalculatePlotSize is used to
%    minimize the allowed point plot size. This function is still a bit
%    messy in terms of visuals. Although, it also does not change the math,
%    just the visuals.
%
% 4. As the size of the rectangle/cube of the figure increases, so does 
%    the error percentage, due to the fact that the plot point always stay
%    the same size. For example, if 2 spheres have a radius of 1 million 
%    each, the sample size will need to be extremely large, nearing 1 
%    billion to maintain accuracy. This is the opposite for radius that 
%    are decreasing in size, in which they might also be inaccurate, but 
%    instead, due to the fact that the plot points may be too large for 
%    the size of the circles/ spheres. I have not yet input something to 
%    counteract this, but am continuing to work on it.
%
% 5. Both problems 3 and 4 are interchangable in a way. Although for the
%    sake of this project, the each part of the project should work just
%    fine.


%%                                                         Important Notes
% ------------------------------------------------------------------------


% Array should be set up like below (C1 = Circle 1, C2 = Circle 2, so on)
% Array = [ RADIUS(C1) Xcenter(C1) CENTERY(C1) Zcenter(C1) 1(HOLE)
%           RADIUS(C2) Xcenter(C2) CENTERY(C2) Zcenter(C2) 0(AREA)
%           RADIUS(C3) Xcenter(C3) CENTERY(C3) Zcenter(C3) 1(HOLE) ]
%
% All arrays should have the 5 columns that are described as seen above. 
% If a circle is a hole, then it should be "1", but if it takes up space, 
% it should be a "0".

% Dimensions should be either "2" or "3"

% SampleSize should be a positive integer. Any sample size from 1 to 10
% million should run in less than a minute.

% The NumberOfLoops should be a positive integer. When NumberOfLoops is
% more than one, the figures only plot the last loop to maximize speed.

% If you would like the SampleSize to change along with the number of 
% runs, input an Array like [1000 2000 3000], in which the length of the
% Array is the same size as the NumberOfLoops.

% ALL OF THE PROJECT PARTS WITH FOR LOOPS ARE CURRENT WORKS IN PROGRESS
% THAT MIGHT HAVE TO BE EDITED FOR A DIFFERENT NUMBER OF CIRCLES/SPHERES

%}
%%               Running a For Loop Through the Function for Project Parts
% ------------------------------------------------------------------------
clf;close all;

% Everything that is correlated to the for loop that is running through 
% the whole function will be outline (similar to MLH later) like
% ------------------------------------------------ Function For Loop Start
if isempty(NumberOfLoops)
    NumberOfLoops = max(size(SampleSize));
end

Dimensions = Dim;

if isempty(d)
    d = 1;
    Array(2,2) = d*Array(2,2);
else
    din = zeros(1,NumberOfLoops);
    for i = 1:(NumberOfLoops-1)
    din(1,i+1) = (d/(NumberOfLoops-1))*i;
    end
end

x_plot = din;
y_plot = zeros(NumberOfLoops, 25);
for NOR2 = 1:25
for NOR = 1:NumberOfLoops
    Array(2,2) = din(NOR);

            
% and will end with
% -------------------------------------------------- Function For Loop End


%%                                                          Project Part 1


%%                          Making Sure Inputs of the Function Are Correct
% ------------------------------------------------------------------------


% Finding the dimensions of the array input
SA = size(Array);

% Making sure that the dimensions of the array input are correct
if SA(2) ~= 5
    error(['The length of the Array must be equal to 5, ' ...
           'thus following the correspoonding data format.'])
end

% Making sure that the dimensions input is correct
if Dimensions ~= 2 && Dimensions ~= 3
    error(['The Diemnsions must be "3"(3-Dimensional) or ' ...
           '"2" (2-Dimensional).'])
end

% Making sure that all plots in the Hole/Area Column are either 0 or 1,
% with holes being 1 and areas being 0.
for i = 1:SA(1)
    if Array(i,5) ~= 0 && Array(i,5) ~= 1
        error('All numbers in column 5 must be equal to 0 or 1')
    end
end

% Making sure the radius of each circle is positive
for i = 1:SA(1)
    if Array(i,1) < 0
        error(['The radius of circle ', num2str(i), ' must be positive.'])
    end
end

% Making sure that the sample size is a positive integer
SS = SampleSize(NOR);
% ------------------------------------------------ Function For Loop Start
if ismatrix(SampleSize)
    SSSA = size(SampleSize); % SSSA = Sample Size Size Array
    if max(SSSA) ~= NumberOfLoops
        error(['The size of SampleSize input must be the same size as' ...
               ' NumberOfLoops input.'])
    end
% -------------------------------------------------- Function For Loop End
elseif round(SS) ~= SS || SS < 1
    error('The sample size must be a positive integer.')
end



%%                                                           Limits Holder
% ------------------------------------------------------------------------


% Preallocating a cell array that holds individual limits of circles due 
% to change in dx, dy, and dz. Mainly used for visuals.
NOC = SA(1); % NOC: Number of Circles
LH = zeros((NOC), 10); % LH: Limits Holder
% Representation of the limits holder
%[ CntrX1 CntrY1 CntrZ1 X1min X1max Y1min Y1max Z1min Z1max Hole/Area]
%[ CntrX2 CntrY2 CntrZ2 X2min X2max Y2min Y2max Z2min Z2max Hole/Area]
%[ AND SO ON                                                         ]

% For loop that automates the x,y, and or z center/limits
for i = 1:NOC  
    LH(i,1) = (              Array(i,2) );
    LH(i,2) = (              Array(i,3) );
    LH(i,3) = (              Array(i,4) );
    LH(i,4) = (-Array(i,1) + Array(i,2) );
    LH(i,5) = ( Array(i,1) + Array(i,2) );
    LH(i,6) = (-Array(i,1) + Array(i,3) );
    LH(i,7) = ( Array(i,1) + Array(i,3) );
    LH(i,8) = (-Array(i,1) + Array(i,4) );
    LH(i,9) = ( Array(i,1) + Array(i,4) );

    if Array(i,5) == 1 % Marking the 10th column as a hole or not
        LH(i,10) = 1;
    end

end

% Setting Variables of the total X, Y, and Z, min's and max's.
% This will be used later to set the visual limits and the limits in which
% the dots can be randomized.
% T = total (xmin, xmax, ymin, ymax, zmin, zmax)
TXmin = min(LH(:,4)); 
TXmax = max(LH(:,5));
TYmin = min(LH(:,6));
TYmax = max(LH(:,7));

if Dimensions == 3
    TZmin = min(LH(:,8));
    TZmax = max(LH(:,9));
end


% Minimized Limits Holder -------------------------------------- Start MLH
% In cases where a hole extends the size of a figure, we can use the
% upcoming code to minimize the limits. If applicable this can be seen in 
% figure(2). If this is not applicable to the information given, this will
% not run, as it does not change the limits. The minimized version is 
% useful to increase accuracy of finding area by decreasing the range in 
% which the plots are randomized, thus increasing the chance of plots 
% landing in areas, rather than holes. Sets of code that that use MLH are
% commented at the start and end of the line with "Start MLH" and 
% "End MLH" repectively, for an easier read of the code.

% Finding the number of non-hole circles/spheres
NHC = NOC - sum(LH(:,10));

% Preallocating an array that holds limits of non-hole circles/spheres
MLH = zeros((NHC), 9); % MLH: Minimized Limits Holder

% Plotting only the non-hole limits in the minimized limits holder
k = 1; % k is used to add limits to the MLH
for i = 1:NOC
    if LH(i,10) == 0 % This checks if hole is true or not
        MLH(k,1:9) = LH(i, 1:9);

        k = k + 1;
    end
end

% Finding the x and y minimized and maximized limits
MXmin = min(MLH(:,4));
MXmax = max(MLH(:,5));
MYmin = min(MLH(:,6));
MYmax = max(MLH(:,7));

% Finding the z minimized and maximized limits if in 3-D
if Dimensions == 3
    MZmin = min(MLH(:,8));
    MZmax = max(MLH(:,9));
end


% ---------------------------------------------------------------- End MLH


%%                                            Deciding Which Limits to Use
% ------------------------------------------------------------------------


% This is used in later lines to decide whether to use T or M limits (2D)
if TXmin ~= MXmin || TXmax ~= MXmax || TYmin ~= MYmin || TYmax ~= MYmax 
    MvsT = 1; % Marks as the use of minimized limits
    UXmin = MXmin;
    UXmax = MXmax;
    UYmin = MYmin;
    UYmax = MYmax;
else
    MvsT = 0; % Marks as the use of total limits
    UXmin = TXmin;
    UXmax = TXmax;
    UYmin = TYmin;
    UYmax = TYmax;
end

% This is used in later lines to decide whether to use T or M limits (3D)
if Dimensions == 3
    if TZmin ~= MZmin || TZmax ~= MZmax || MvsT == 1
        MvsT = 1;
        UZmin = MZmin;
        UZmax = MZmax;
    else
        UZmin = TZmin;
        UZmax = TZmax;
    end
end


%%                                       Plotting Visuals of Circles (2-D) 
% ------------------------------------------------------------------------


% ------------------------------------------------ Function For Loop Start
if NOR2 == 25 && NOR == NumberOfLoops % Only visualizing the last loop  
% -------------------------------------------------- Function For Loop End


% This is used so the circles/sphere are the same colors in figures 1 and 2
ColorHolder = zeros(NOC,3);

figure(1)
% Plotting (2-D) circles as visuals
if Dimensions == 2
    % Plotting the circles that cover area
    for i = 1:NOC            
        if Array(i,5) == 0 % Making sure the circle is not a hole
    
            % Visually pltos the circles given
            circle = nsidedpoly(1000, 'Center', ...
                [Array(i,2) Array(i,3)], 'Radius', Array(i,1));
                
            % Randomizing colors of the cirlces into the Color Holder
            ColorHolder(i,:) = rand(1,3);
    
            % Plots the circle as a random color that is semi-transparent
            plot(circle, 'FaceColor', ColorHolder(i,:), ...
                'FaceAlpha', 0.75)
    
        end
        hold on
    end
    
    figure(1)
    % Only plotting the outline of the circles that cover area to see the
    % overlap on the original figure to show resemblance to figure(2)
    for i = 1:NOC            
        if Array(i,5) == 0  % Making sure the circle is not a hole
    
            % Visually plots the circles given
            circle = nsidedpoly(1000, 'Center', ...
                [Array(i,2) Array(i,3)], 'Radius', Array(i,1));
    
            % Plots the circle as a random color that is semi-transparent
            plot(circle, 'FaceColor', ColorHolder(i,:), 'FaceAlpha', 0)
        
        end
        hold on
    end

    % Plotting the circles that are holes last so that anything under
    % them are not visible
    for i = 1:NOC
        if Array(i,5) == 1 % Making sure the circle is a hole
            
            % Visually plots the circles given
            circle = nsidedpoly(1000, 'Center', ...
                [Array(i,2) Array(i,3)], 'Radius', Array(i,1));
    
            % Plots the hole as white with a black outline
            plot(circle, 'FaceColor', 'white', 'FaceAlpha', 1)
        
        end
        hold on
    end
    
    % Applying the limits of the figure visually
    xlim([(TXmin)   (TXmax)])
    ylim([(TYmin)   (TYmax)])
    axis equal
    
    % Applying the limits of the plot points (in red)
    line(([TXmin TXmin]), ([TYmin,TYmax]), 'Color', 'Red')
    line(([TXmax TXmax]), ([TYmin,TYmax]), 'Color', 'Red')
    line(([TXmin TXmax]), ([TYmin,TYmin]), 'Color', 'Red')
    line(([TXmin TXmax]), ([TYmax,TYmax]), 'Color', 'Red')
        

    % MLH Plotting another figure expressing MLH---------------- Start MLH
    
    
    if MvsT == 1 % Using MLH limits
        
    figure(2)
    % Plotting the circles that cover area
    for i = 1:NOC            
        if Array(i,5) == 0  % Making sure the circle is not a hole
    
            % Visually plots the circles given
            circle = nsidedpoly(1000, 'Center', ...
                [Array(i,2) Array(i,3)], 'Radius', Array(i,1));
    
            % Plots the circle as a random color that is semi-transparent
            plot(circle, 'FaceColor', ColorHolder(i,:), 'FaceAlpha', 1)
        
        end
        hold on
    end
        
    figure(2)
    % Plotting the circles that are holes last so that anything under
    % them are not visible
    for i = 1:NOC
        if Array(i,5) == 1
            
            % Visually plots the circles(holes) given
            circle = nsidedpoly(1000, 'Center', ...
                [Array(i,2) Array(i,3)], 'Radius', Array(i,1));
    
            % Plots the hole as white with a black outline
            plot(circle, 'FaceColor', 'white', 'FaceAlpha', 1)
            
        end
        hold on
    end

        % Applying the limits of the figure visually
        axis equal
        xlim([(MXmin)   (MXmax)])
        ylim([(MYmin)   (MYmax)])
        
        % Applying the limits of the plot points (in red)
        line(([MXmin MXmin]), ([MYmin,MYmax]), 'Color', 'Red')
        line(([MXmax MXmax]), ([MYmin,MYmax]), 'Color', 'Red')
        line(([MXmin MXmax]), ([MYmin,MYmin]), 'Color', 'Red')
        line(([MXmin MXmax]), ([MYmax,MYmax]), 'Color', 'Red')
    end


    % ------------------------------------------------------------ End MLH


end


%%                                       Plotting Visuals of Spheres (3-D)
% -----------------------------------------------------------------------


% Plotting (3-D) spheres for visuals
if Dimensions == 3

    % Plotting the spheres that cover area
    for i = 1:NOC            
        if Array(i,5) == 0  % Making sure the sphere is not a hole
    
            [X,Y,Z] = sphere(20); % 20 faces on the given spheres
            X = X*Array(i,1); % Defining X values
            Y = Y*Array(i,1); % Defining Y values
            Z = Z*Array(i,1); % Defining Z values
    
            % Randomizing colors of the spheres into the Color Holder
            ColorHolder(i,:) = rand(1,3);
    
            % Plotting the sphere as a random color that is
            % semi-transparent
            surf(X+Array(i,2), Y+Array(i,3), Z+Array(i,4), ...
                'FaceColor', ColorHolder(i,:), 'FaceAlpha', 0.35, ...
                'EdgeAlpha', 0.35)
            xlabel('X axis')
            ylabel('Y axis')
            zlabel('Z axis')
        
        end
        hold on
    end
       
    % Plotting the spheres that are holes last so that anything under
    % them are not visible
    for i = 1:NOC
        if Array(i,5) == 1
            
            [X,Y,Z] = sphere(20); % 20 faces on the sphere
            X = X*Array(i,1); % Defining X values
            Y = Y*Array(i,1); % Defining Y values
            Z = Z*Array(i,1); % Defining Z values
    
            % Plotting the hole as black
            surf(X+Array(i,2), Y+Array(i,3), Z+Array(i,4), ...
                'FaceColor', 'k', 'EdgeColor', 'w', 'EdgeAlpha', 0.2)
            
        end
        hold on
    end

    % Applying the limits of the figure visually
    xlim([(TXmin)   (TXmax)])
    ylim([(TYmin)   (TYmax)])
    zlim([(TZmin)   (TZmax)])
    axis equal
    
    
    % Applying the limits of the plot points (in red) (3D)
    line(([TXmin TXmin]), ([TYmin,TYmax]), ([TZmin TZmin]), 'Color', 'Red')
    line(([TXmax TXmax]), ([TYmin,TYmax]), ([TZmin TZmin]), 'Color', 'Red')
    line(([TXmin TXmax]), ([TYmin,TYmin]), ([TZmin TZmin]), 'Color', 'Red')
    line(([TXmin TXmax]), ([TYmax,TYmax]), ([TZmin TZmin]), 'Color', 'Red')
    line(([TXmin TXmin]), ([TYmin,TYmax]), ([TZmax TZmax]), 'Color', 'Red')
    line(([TXmax TXmax]), ([TYmin,TYmax]), ([TZmax TZmax]), 'Color', 'Red')
    line(([TXmin TXmax]), ([TYmin,TYmin]), ([TZmax TZmax]), 'Color', 'Red')
    line(([TXmin TXmax]), ([TYmax,TYmax]), ([TZmax TZmax]), 'Color', 'Red')
    line(([TXmin TXmin]), ([TYmax,TYmax]), ([TZmin TZmax]), 'Color', 'Red')
    line(([TXmax TXmax]), ([TYmax,TYmax]), ([TZmin TZmax]), 'Color', 'Red')
    line(([TXmin TXmin]), ([TYmin,TYmin]), ([TZmin TZmax]), 'Color', 'Red')
    line(([TXmax TXmax]), ([TYmin,TYmin]), ([TZmin TZmax]), 'Color', 'Red')

    
    % MLH Plotting another figure expressing MLH---------------- Start MLH


    if MvsT == 1 % Using MLH limits                  

        figure(2)
        % Plotting the spheres that cover area
        for i = 1:NOC            
            if Array(i,5) == 0  % Making sure the circle is not a hole

                [X,Y,Z] = sphere(20); % 20 faces on the sphere	
                X = X*Array(i,1); % Defining X values	
                Y = Y*Array(i,1); % Defining Y values	
                Z = Z*Array(i,1); % Defining Z values	

                % Plotting the sphere as a random color	
                surf(X+Array(i,2), Y+Array(i,3), Z+Array(i,4), ...	
                    'FaceColor', ColorHolder(i,:), 'FaceAlpha', 0.35, ...	
                    'EdgeAlpha', 0.35)
                xlabel('X axis')
                ylabel('Y axis')
                zlabel('Z axis')

            end
            hold on
        end
        
        % Plotting the spheres that are holes last so that anything under
        % them are not visible
        for i = 1:NOC
            if Array(i,5) == 1
                
                [X,Y,Z] = sphere(20); % 20 faces on the sphere	
                X = X*Array(i,1); % Defining X values	
                Y = Y*Array(i,1); % Defining Y values	
                Z = Z*Array(i,1); % Defining Z values	
                
                % Plotting the hole as black	
                surf(X+Array(i,2), Y+Array(i,3), Z+Array(i,4), ...	
                    'FaceColor', 'k', 'EdgeColor', 'w', 'EdgeAlpha', 0.2)
                
                
            end
            hold on
        end
            
    % Applying the limits of the figure visually
    axis equal
    xlim([(MXmin)   (MXmax)])
    ylim([(MYmin)   (MYmax)])
    zlim([(MZmin)   (MZmax)])
    
    % Applying the limits of the plot points (in red)(3D)
    line(([MXmin MXmin]), ([MYmin,MYmax]), ([MZmin MZmin]), 'Color', 'Red')
    line(([MXmax MXmax]), ([MYmin,MYmax]), ([MZmin MZmin]), 'Color', 'Red')
    line(([MXmin MXmax]), ([MYmin,MYmin]), ([MZmin MZmin]), 'Color', 'Red')
    line(([MXmin MXmax]), ([MYmax,MYmax]), ([MZmin MZmin]), 'Color', 'Red')
    line(([MXmin MXmin]), ([MYmin,MYmax]), ([MZmax MZmax]), 'Color', 'Red')
    line(([MXmax MXmax]), ([MYmin,MYmax]), ([MZmax MZmax]), 'Color', 'Red')
    line(([MXmin MXmax]), ([MYmin,MYmin]), ([MZmax MZmax]), 'Color', 'Red')
    line(([MXmin MXmax]), ([MYmax,MYmax]), ([MZmax MZmax]), 'Color', 'Red')
    line(([MXmin MXmin]), ([MYmax,MYmax]), ([MZmin MZmax]), 'Color', 'Red')
    line(([MXmax MXmax]), ([MYmax,MYmax]), ([MZmin MZmax]), 'Color', 'Red')
    line(([MXmin MXmin]), ([MYmin,MYmin]), ([MZmin MZmax]), 'Color', 'Red')
    line(([MXmax MXmax]), ([MYmin,MYmin]), ([MZmin MZmax]), 'Color', 'Red')
    end


    % ------------------------------------------------------------ End MLH


end


% ------------------------------------------------ Function For Loop Start
end % Only visualizing the last loop of the for loop
% -------------------------------------------------- Function For Loop End


%%                                                  Making the Plot Points
% ------------------------------------------------------------------------


% If the limits use the MLH, then we make the RAD, defined below, based on
% the MLH, and along with that, we randomize plots within the MLH, rather
% than the original limits. We have already defined whether we will use
% the total (T) or minimized (M) limits, by using U.
if Dimensions == 2
    % RAD = Rectangle Area based on the limits of the circle(s)
    RAD = (UXmax - UXmin) * (UYmax - UYmin);
else % Dimensions = 3
    % RAD is now a cube based the on limits of the sphere(s)
    RAD = (UXmax - UXmin) * (UYmax - UYmin) * (UZmax - UZmin);
end

    
% Preallocating array of every X, Y, and Z Dot for randomization
EXYZD = zeros(SS, NOC+5); % EXYZD = Every X, Y, and Z Dot
% EXYZD LOOKS LIKE (C = Circle)
%                          1: In Any C    1: In +1 C   1: In Any C
%                          0: No C        0: No C      0: No C
% [Xplot   Yplot   Zplot   NOC(i)(1/0)^   (1/0)^       (1/0)^     ]

    % Randomizing the dots from minimum limits to maximum limits
    EXYZD(:,1) = (UXmin) + (UXmax-UXmin) .* rand(SS,1); % X plot
    EXYZD(:,2) = (UYmin) + (UYmax-UYmin) .* rand(SS,1); % Y plot
if Dimensions == 3
    EXYZD(:,3) = (UZmin) + (UZmax-UZmin) .* rand(SS,1); % z plot
end


%%             Determining Whether the Plots Reside in Circle(s)/Sphere(s)
% ------------------------------------------------------------------------


for ii = 1:SS     % Goes through each of the dots  
    
    for i = 1:NOC % Goes through each of the circles/spheres to check if 
                  % the dot resides in said circle/sphere
    
        DotX = EXYZD(ii,1) - LH(i,1); % Finding x distance from center
        DotY = EXYZD(ii,2) - LH(i,2); % Finding y distance from center
    if Dimensions == 3
        DotZ = EXYZD(ii,3) - LH(i,3); % Finding z distance from center
    end

    
    if Dimensions == 2
        % Finding distance to center of circle(i) using x and y coords
        DotXYZ = sqrt(DotX^2 + DotY^2);
    else
        % Finding distance to center of circle(i) using x, y, and z coords
        DotXYZ = sqrt(DotX^2 + DotY^2 + DotZ^2);
    end


    % Decides if the dot is within the circles radius and assigns 1 as
    % "true" in its own circle/sphere column of EXYZD
    if DotXYZ < Array(i,1) 
        EXYZD(ii, i+3) = 1;
    
        % Shows that this plot resides in area. This is later used for
        % calculating the area of which at least one circle/sphere resides
        if EXYZD(ii, end) == 0
            EXYZD(ii, end) = 1;
        end
        
        % All "true" statements from other if statements are ignored and 
        % will continue to be ignored due to fact that the plot resides in 
        % a hole. This also breaks the for loop to examine the next plot. 
        if Array(i,5) == 1
            EXYZD(ii, (4:end)) = 0;
        break
        end  
    
    end


    end % end of 2nd for loop

end % end of 1st for loop


%%        Finding # of Plots and the Area of the Shared and Singular Space
% ------------------------------------------------------------------------
        

% Finding the percentage of dots that resides in more than 1 circle of 
% which none of them are holes
for i = 1:SS
    
    % Assigning one or "true" to the second to last column if a dot 
    % resides in more than 1 circle, of which none of them are holes.
    if sum(EXYZD(i, 4:(NOC+3))) >= 2
        EXYZD(i, end-1) = 1; 
    end                     

end

NODIMTOC  = sum(EXYZD(:,end-1)); % Number of dots in atleast two circles    import
NODIMTOCA = NODIMTOC/SS * RAD  ; % Area of dots in at least two circles
% disp(['The area of which at least two circles reside in is roughly ', ...
%       num2str(NODIMTOCA), ' units.'])

NODIAON   = sum(EXYZD(:,end  )); % Number of dots in atleast one circle
NODIAONA  = NODIAON /SS * RAD  ; % Area of dots in at least one circle
% disp(['The area of which at least one circle resides in is roughly ', ...
%       num2str(NODIAONA), ' units.'])

NodRad = NODIAON / RAD; % This is used along with the function below to 
% minimize the number of plot points that are being plotted, based on the 
% size of the figure and the size of the area in which more than one
% circles/spheres reside. This is done the minimize run speed, while 
% keeping visuals intact. [PlotSize] is the variable that holds that 
% number of points that will be plotted.  
[PlotSize] = CalculatePlotSize(SS, NodRad, NODIAON);


%%                                     Graphing Dots On Top Of Other Graph
% ------------------------------------------------------------------------


% ------------------------------------------------ Function For Loop Start
if NOR2 == 25 && NOR == NumberOfLoops % Only visualizing the last loop 
% -------------------------------------------------- Function For Loop End


% Used to limit the number of plotted points to minimized run time
PlotSizeNow = 0;

% Plotting the coordinates of the plots that reside in more than 1 circles
if MvsT == 0 % Making sure 
    figure(1)
    if Dimensions == 2


        for i = 1:SS
            if EXYZD(i,end-1) == 1 % Making sure the plot resides in more 
                                   % than one circle

                % Plotting the 2D points in red
                plot(EXYZD(i,1), EXYZD(i,2), 'ro-')
                

                PlotSizeNow = PlotSizeNow + 1;
                if PlotSizeNow == PlotSize
                    break % Breaking the for loop the PlotSize
                end
                hold on

            end
        end


    else % Dimensions = 3


        for i = 1:SS
            if EXYZD(i,end-1) == 1 % Making sure the plot resides in more 
                                   % than one sphere
                
                % Plotting the 3-D points in red
                plot3(EXYZD(i,1), EXYZD(i,2), EXYZD(i,3), 'ro-')

                PlotSizeNow = PlotSizeNow + 1;
                if PlotSizeNow == PlotSize
                    break % Breaking the for loop at PlotSize
                end
                hold on

            end
        end


    end


    % Defining titles for respective dimensions
    if Dimensions == 2
        title('Areas Of and Between Circles (Affected by Holes)')
    else
        title('Areas Of and Between Spheres (Affected by Holes)')
    end


    % Shows text on the x label which shows the areas
    if Dimensions == 2
        AAA = ['In red, the area of which more than one circle re' ...
                'sides in is roughly '];
        CCC = ['In randomized colors, the areas of which at least ' ...
                'one circle resides in is roughly '];
    else
        AAA = ['In red, the area of which more than one sphere re' ...
                'sides in is roughly '];
        CCC = ['In randomized colors, the areas of which at least ' ...
                'one sphere resides in is roughly '];
    end

    BBB = num2str(NODIMTOCA);
    DDD = num2str(NODIAONA);
    pos = [0 0 475 0];
    h = uicontrol('Style','Text','Position',pos);
    text = {[AAA BBB ' units.']; [CCC DDD ' units.']};
    xlabel(textwrap(h, text), 'FontSize', 10)


else % -------------------------------------------------------- Start MLH
    

    % There is now a normal limits figure and a minimized limits figure

    figure(1) % figure 1 title
    if Dimensions == 2
        title('Areas Of and Between Circles (Affected by Holes)')
    else
        title('Areas Of and Between Spheres (Affected by Holes)')
    end

    figure(2)
    if Dimensions == 2


        for i = 1:SS
            if EXYZD(i,end-1) == 1 % Making sure the plot resides in more
                                   % than one circle

                % Plotting the 2-D point in red
                plot(EXYZD(i,1), EXYZD(i,2), 'ro-')

                PlotSizeNow = PlotSizeNow + 1;
                if PlotSizeNow == PlotSize
                    break % Breaking the for loop at PlotSize
                end
                hold on

            end
        end

        
    else % Dimensions = 3


        for i = 1:SS
            if EXYZD(i,end-1) == 1 % Making sure the plot resides in more
                                   % than one sphere

                % Plotting the 3-D point in red
                plot3(EXYZD(i,1), EXYZD(i,2), EXYZD(i,3), 'ro-')

                PlotSizeNow = PlotSizeNow + 1;
                if PlotSizeNow == PlotSize
                    break % Breaking the for loop at PlotSize
                end
                hold on

            end
        end


    end
    

    % Defining titles for respective dimensions (figure(2))
    if Dimensions == 2
        title('Areas Of and Between Circles (Affected by Holes)')
    else
        title('Areas Of and Between Spheres (Affected by Holes)')
    end
    

    % Shows text on the x label which shows the areas
    if Dimensions == 2
        AAA = ['In red, the area of which more than one circle re' ...
                'sides in is roughly '];
        CCC = ['In randomized colors, the areas of which at least ' ...
                'one circle resides in is roughly '];
    else
        AAA = ['In red, the area of which more than one sphere re' ...
                'sides in is roughly '];
        CCC = ['In randomized colors, the areas of which at least ' ...
                'one sphere resides in is roughly '];
    end

    BBB = num2str(NODIMTOCA);
    DDD = num2str(NODIAONA);
    pos = [0 0 475 0];
    h = uicontrol('Style','Text','Position',pos);
    text = {[AAA BBB ' units.']; [CCC DDD ' units.']};
    xlabel(textwrap(h, text), 'FontSize', 10)


    % ------------------------------------------------------------ End MLH


end


% ------------------------------------------------ Function For Loop Start
end % Only visualizing the last loop of the for loop

y_plot(NOR,NOR2) = NODIAONA;

end % End of for loop
% figure(3)
% plot( x_plot, y_plot(:,NOR2), '|', 'Color', 'Blue', 'LineWidth', 1)
% hold on
end % End of for loop
% ymean = mean(y_plot,2);
% plot( x_plot, ymean, 'r', 'linewidth', 0.1)
% -------------------------------------------- Project Part 1 and 2 Graphs



               


%%         Function to Minimize Number of Points Plotted, Reducing Run Time

function [PlotSize] = CalculatePlotSize(SS, NodRad, NODIAON)
% This is a function that counts the number of plots that would originally
% be plotted and minimizes it in response to the magnitude of the sample
% size, thus minimizing the speed taken to run the program, while keeping
% the visuals intact.


    if SS >= 0 && SS < 10000
            PlotSize = NODIAON;
%     elseif SS >= 1000      && SS < 10000
%         if NodRad <= 0 && NodRad < 0.25
%             PlotSize = NODIAON;
%         else
%             PlotSize = NODIAON;
%         end
    elseif SS >= 10000     && SS < 100000
        if NodRad <= 0 && NodRad < 0.25
            PlotSize = NODIAON;
        else
            PlotSize = NODIAON/10;
        end
    elseif SS >= 100000    && SS < 1000000
        if NodRad <= 0 && NodRad < 0.25
            PlotSize = NODIAON/10;
        else
            PlotSize = NODIAON/100;
        end
    elseif SS >= 1000000   && SS < 10000000
        if NodRad <= 0 && NodRad < 0.25
            PlotSize = NODIAON/100;
        else
            PlotSize = NODIAON/1000;
        end
    elseif SS >= 10000000  && SS < 100000000
        if NodRad <= 0 && NodRad < 0.25
            PlotSize = NODIAON/1000;
        else
            PlotSize = NODIAON/10000;
        end
    elseif SS >= 100000000 && SS < 1000000000
        if NodRad <= 0 && NodRad < 0.25
            PlotSize = NODIAON/10000;
        else
            PlotSize = NODIAON/100000;
        end
    end
    PlotSize = round(PlotSize); % Rounding to a hole number so the for 
                                % loops have a number to break at, 
                                % considering the for loops are going from 
                                % 1 to (Positive Integer).
end % End of local function
end % End of function