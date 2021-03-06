%Robert MacCurdy - rbm7@cornell.edu - September 2010

%Newton-Raphson simulation of location-finding for TOA arrays based on 
%methods similar to GPS.
clear

%%%%%%%%%%%%%USER-DEFINED CONSTANTS%%%%%%%%%%%

%define the Tower positions; cols are X Y, rows are tower #; (position 
%defined in meters). Add more rows if more receivers are desired.
TwrPos=[-500 2500;
      -2000 -300;
      900 -2800;
      1600 500];
  
%define tag position for sim  
TagPos=[5550 3000];

%Define actual transmit time, relative to some universal reference. 
%A NOTE ABOUT TIME CALCULATIONS BELOW: though time is
%expressed here in seconds, it is converted into distance in all
%calculations below, so the units are in meters. This enables matrix
%computations with uniform units. 
TrueTransmitTime=-2; %the actual tag transmit time (seconds)

%A NOTE ABOUT COORDINATES AND TIME: all coordinates (X, Y, and time which is
%expressed as a range in meters) are assumed to be local. The algorithm's
%performance could deteriorate if very large values for any of the 
%user-supplied inputs are supplied. This is a matter of numerical 
%precision, particularly within the matrix operations. Therefore, if this 
%code is used outside this example, on real data expressed in UTM or 
%seconds from the beginning of the week (as GPS uses), a wrapper routine 
%should translate those absolute units into local positions and times.

%define the Initial Guess for position and transmit time. 
TagPosGuess = [0,0]; %Solution converges if init guess is inside array, but
                     %don't choose an init position too close to a
                     %receiver, or the matrix G will be ill-conditioned.
TransmitTimeGuess=0; %the time that the towers think the tag sent the 
%message (seconds)

ViewRegionSize=3; %scaling factor to show more or less of region outside array


%%%%%%NO USER-MODIFIABLE CONSTANTS BELOW%%%%%%%%
c=299792458; %speed of propagation is 299792458 m/sec
numTwrs=size(TwrPos,1);
guess=[TagPosGuess TransmitTimeGuess.*c]; %this seems to need to be INSIDE the array for best results (sometimes works outside array)

%define the search area (used mainly for display purposes); col's are 
%X Y, rows are min to max (meters)
Area=[min(TwrPos(:,1))*ViewRegionSize min(TwrPos(:,2))*ViewRegionSize;...
      max(TwrPos(:,1))*ViewRegionSize max(TwrPos(:,2))*ViewRegionSize];

%show a plot of the scenario
figure1 = figure;
axes1 = axes('Parent',figure1);
xlim(axes1,[Area(1,1) Area(2,1)]);
ylim(axes1,[Area(1,2) Area(2,2)]);
box(axes1,'on');
hold(axes1,'all');
%plot Tower Positions (in XY plane)
for i=1:size(TwrPos,1)
    plot(TwrPos(i,1),TwrPos(i,2),'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],...
        'MarkerSize',8,...
        'Marker','diamond',...
        'LineStyle','none');
end
%show tag position
plot(TagPos(1),TagPos(2),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
    'MarkerSize',12,...
    'Marker','o',...
    'LineStyle','none');
    
%get true range and signal flight time from each tower to tag
trueRange=sqrt(sum(((ones(numTwrs,1)*TagPos)-TwrPos).^2,2)); 
TrueTravelTime=trueRange./c;

%this is the time that each tower will detect a transmission, based on the
%true position of the tag.
MeasuredArrivalTime=TrueTravelTime+TrueTransmitTime;

%plot the initial guess
plot(guess(1),guess(2),'MarkerFaceColor',[1 .7 0],'MarkerEdgeColor',[1 .7 0],...
    'MarkerSize',8,...
    'Marker','diamond',...
    'LineStyle','none');
Soln_Track(1,1)=guess(1);
Soln_Track(2,1)=guess(2);


range = zeros(numTwrs,1);
stop = 0;

while (stop == 0)
    % create guess matrix
    guessmatrix = [];
    for i = 1:numTwrs
        guessmatrix = [ guessmatrix; guess ];
    end
    
    %find the vestor difference between the assumed tag position and each
    %of the towers.
    TagtoTwrVec = TwrPos - guessmatrix(:,1:2);
    
    %find the ranges from the assumed transmit position to each tower
    guessrange = sqrt(sum((TagtoTwrVec.^2)')');
    
    %find the difference between the measured ranges (based on guessed Tx
    %time), and the assumed ranges (based on assumed position). If we have
    %guessed correctly, the difference should be zero.
    del_p = (MeasuredArrivalTime.*c - guessrange - guess(3)); %the pseudorange difference
    
    pmatrix = guessrange * ones(1,2);
    G = TagtoTwrVec ./ pmatrix;
    G = [ -G ones(numTwrs,1) ];
    % solve for  'deltaPos' which contains dx, dy, and  dt
    deltaPOS = G \ del_p;
    % calculate 'obsPos' by adding 'deltaPos' to the current guess
    obsPos = guess + deltaPOS';
    % check to see if the initial guess and the computed result is
    % "close enough"; if it is, then stop the iteration by setting the
    % 'stop' flag to 1; else, set the new initial guess equal to the
    % last computed value, and iterate again
    
    %show computed intermediate position guess
    plot(obsPos(1),obsPos(2),'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0],...
        'MarkerSize',8,...
        'Marker','diamond',...
        'LineStyle','none');
    Soln_Track(1,end+1)=obsPos(1);
    Soln_Track(2,end)=obsPos(2);
    
    %the stop condition is based on position accuracy
    if ((abs(obsPos(1:2) - guess(1:2))) < 1e-6)
        stop = 1;
    end
    %update the last guess position & time estimate with the current one
    guess = obsPos; %obsPos is the observed position (and transmit time). 
                    %The first two entries are X and Y; the last is the
                    %transmit time expressed in meters. Divide by the speed
                    %of light to get the transmit time in seconds.
end

%plot the trajectory of the intermediate solutions
plot(Soln_Track(1,:),Soln_Track(2,:), 'LineStyle','--');

axis(axes1,'square');

%Print the final position estimate
obsPos(1:2)
%Print the final absolute tag transmit time.
obsPos(3)/c


%Now look at PDOP (position dilution of precision) for the system above. Plot it on a
%figure with the tower positions
x=linspace(Area(1,1),Area(2,1),100);
y=linspace(Area(1,2),Area(2,2),100);
z=zeros(size(x,2),size(y,2));

for i=1:size(x,2)
    for j=1:size(y,2)
        guess=[x(i) y(j) 0];
        guessmatrix = ones(numTwrs,1)*guess;
        TagtoTwrVec = TwrPos - guessmatrix(:,1:2);
        guessrange = sqrt(sum((TagtoTwrVec.^2)')');
        pmatrix = guessrange * ones(1,2);
        G = TagtoTwrVec ./ pmatrix;
        G = [ -G ones(numTwrs,1) ];
        H=inv(G'*G);
        z(j,i)=sqrt(H(1,1)+H(2,2));
        if z(j,i)>=50
            %z(j,i)=50; %scale data to fit nicely in color map by compressing values over 50
            z(j,i)=50-log2(50)+log2(z(j,i)); %scale data to fit nicely in color map by compressing values over 50
        end
    end
end

v = linspace(.5,max(max(z)),20);
%v = [.9 1.5 2 2.5 3 3.5 5 7 8 13];

%show a plot of the scenario
figure2 = figure;
axes2 = axes('Parent',figure2,'PlotBoxAspectRatio',[1 1 1],...
    'FontWeight','bold',...
    'FontSize',14);
xlim(axes2,[Area(1,1) Area(2,1)]);
ylim(axes2,[Area(1,2) Area(2,2)]);
box(axes2,'on');
hold(axes2,'all');
%plot PDOP
contourf(x,y,z,v); %only use as many contour lines as necessary
colormap(gray); %display in greyscale colors
%plot Tower Positions (in XY plane)
for i=1:size(TwrPos,1)
    plot(TwrPos(i,1),TwrPos(i,2),'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 1 1],...
        'MarkerSize',13,...
        'Marker','^',...
        'LineStyle','none');
end
% Create xlabel
xlabel('Easting (m)','FontWeight','bold','FontSize',14);

% Create ylabel
ylabel('Northing (m)','FontWeight','bold','FontSize',14);

% Create colorbar
colorbar('peer',axes2);

axis(axes2,'square');
