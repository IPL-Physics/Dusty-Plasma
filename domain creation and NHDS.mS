%Nonhomogeneous Drift Subtraction
%Breaking into Subsections
tracks = ma.tracks; %for 12 domains

xdomain=0:3750:15000;
ydomain=0:700:2100;
dom11tracks=[];
dom12tracks=[];
dom13tracks=[];
dom14tracks=[];

dom21tracks=[];
dom22tracks=[];
dom23tracks=[];
dom24tracks=[];

dom31tracks=[];
dom32tracks=[];
dom33tracks=[];
dom34tracks=[];


for i=1:numel(tracks)
    if xdomain(1)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(2) && ydomain(1)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(2)
        dom11tracks{i}=tracks{i};
    elseif xdomain(2)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(3) && ydomain(1)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(2)
        dom12tracks{i}=tracks{i};
    elseif xdomain(3)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(4) && ydomain(1)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(2)
        dom13tracks{i}=tracks{i};
    elseif xdomain(4)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(5) && ydomain(1)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(2)
        dom14tracks{i}=tracks{i};

    elseif xdomain(1)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(2) && ydomain(2)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(3)
        dom21tracks{i}=tracks{i};
    elseif xdomain(2)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(3) && ydomain(2)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(3)
        dom22tracks{i}=tracks{i};
    elseif xdomain(3)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(4) && ydomain(2)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(3)
        dom23tracks{i}=tracks{i};
    elseif xdomain(4)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(5) && ydomain(2)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(3)
        dom24tracks{i}=tracks{i};
    
    elseif xdomain(1)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(2) && ydomain(3)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(4)
        dom31tracks{i}=tracks{i};
    elseif xdomain(2)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(3) && ydomain(3)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(4)
        dom32tracks{i}=tracks{i};
    elseif xdomain(3)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(4) && ydomain(3)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(4)
        dom33tracks{i}=tracks{i};
    elseif xdomain(4)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(5) && ydomain(3)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(4)
        dom34tracks{i}=tracks{i};

    end
end


alldomtracks = {dom11tracks,dom12tracks,dom13tracks,dom14tracks,dom21tracks,dom22tracks,dom23tracks,dom24tracks,dom31tracks,dom32tracks,dom33tracks,dom34tracks};

AllDomTracks=[]; %remove empty doubles from each cell and make vertical

for i=1:numel(alldomtracks);
    currentTracks = alldomtracks{i};
    currentTracks = currentTracks(~cellfun('isempty', currentTracks));
    currentTracks = reshape(currentTracks, [], 1);
    AllDomTracks{i}=currentTracks;
end


tracks = ma.tracks; %for 20 domains

xdomain=0:3000:15000; % 5 x domains
ydomain=0:525:2100; % 4 y domains
dom11tracks=[];
dom12tracks=[];
dom13tracks=[];
dom14tracks=[];
dom15tracks=[];

dom21tracks=[];
dom22tracks=[];
dom23tracks=[];
dom24tracks=[];
dom25tracks=[];

dom31tracks=[];
dom32tracks=[];
dom33tracks=[];
dom34tracks=[];
dom35tracks=[];

dom41tracks=[];
dom42tracks=[];
dom43tracks=[];
dom44tracks=[];
dom45tracks=[];


for i=1:numel(tracks)
    if xdomain(1)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(2) && ydomain(1)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(2)
        dom11tracks{i}=tracks{i};
    elseif xdomain(2)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(3) && ydomain(1)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(2)
        dom12tracks{i}=tracks{i};
    elseif xdomain(3)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(4) && ydomain(1)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(2)
        dom13tracks{i}=tracks{i};
    elseif xdomain(4)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(5) && ydomain(1)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(2)
        dom14tracks{i}=tracks{i};
    elseif xdomain(5)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(6) && ydomain(1)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(2)
        dom15tracks{i}=tracks{i};

    elseif xdomain(1)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(2) && ydomain(2)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(3)
        dom21tracks{i}=tracks{i};
    elseif xdomain(2)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(3) && ydomain(2)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(3)
        dom22tracks{i}=tracks{i};
    elseif xdomain(3)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(4) && ydomain(2)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(3)
        dom23tracks{i}=tracks{i};
    elseif xdomain(4)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(5) && ydomain(2)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(3)
        dom24tracks{i}=tracks{i};
    elseif xdomain(5)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(6) && ydomain(2)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(3)
        dom25tracks{i}=tracks{i};
    
    elseif xdomain(1)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(2) && ydomain(3)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(4)
        dom31tracks{i}=tracks{i};
    elseif xdomain(2)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(3) && ydomain(3)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(4)
        dom32tracks{i}=tracks{i};
    elseif xdomain(3)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(4) && ydomain(3)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(4)
        dom33tracks{i}=tracks{i};
    elseif xdomain(4)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(5) && ydomain(3)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(4)
        dom34tracks{i}=tracks{i};
    elseif xdomain(5)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(6) && ydomain(3)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(4)
        dom35tracks{i}=tracks{i};
    
    elseif xdomain(1)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(2) && ydomain(4)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(5)
        dom41tracks{i}=tracks{i};
    elseif xdomain(2)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(3) && ydomain(4)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(5)
        dom42tracks{i}=tracks{i};
    elseif xdomain(3)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(4) && ydomain(4)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(5)
        dom43tracks{i}=tracks{i};
    elseif xdomain(4)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(5) && ydomain(4)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(5)
        dom44tracks{i}=tracks{i};
    elseif xdomain(5)<=tracks{i}(1,2) && tracks{i}(1,2)<=xdomain(6) && ydomain(4)<=tracks{i}(1,3) && tracks{i}(1,3)<=ydomain(5)
        dom45tracks{i}=tracks{i};

    end
end


alldomtracks = {dom11tracks,dom12tracks,dom13tracks,dom14tracks,dom15tracks,dom21tracks,dom22tracks,dom23tracks,dom24tracks,dom25tracks,dom31tracks,dom32tracks,dom33tracks,dom34tracks,dom35tracks,dom41tracks,dom42tracks,dom43tracks,dom44tracks,dom45tracks};

AllDomTracks=[]; %remove empty doubles from each cell and make vertical

for i=1:numel(alldomtracks);
    currentTracks = alldomtracks{i};
    currentTracks = currentTracks(~cellfun('isempty', currentTracks));
    currentTracks = reshape(currentTracks, [], 1);
    AllDomTracks{i}=currentTracks;
end


%Plot a certain amount of tracks of whatever subsection for viewing
%purposes
figure
axis([0 15000 -1000 4000])
N=100; %number of particles in subsection


for i=1:numel(AllDomTracks)
    currentdomtracks=AllDomTracks{i};
    %currentma=Allma{1,i};
    %currentdomtracks=currentma.tracks;
    for j = 1:min(N,numel(currentdomtracks))
        hold on
        plot(currentdomtracks{j}(:, 2), currentdomtracks{j}(:, 3));
    end
   
end

%Tells you how many particle trajectories you now have in AllDomTracks 
%and how many you lost from the original tracks
tracksize=0;
for i=1:numel(AllDomTracks)
    current=AllDomTracks{i};
    n=numel(current);
    tracksize=tracksize+n;
end
allDomainTrackSize = tracksize;
allDomainTrackSize
totalTrackSize = numel(tracks)
LostTracks = totalTrackSize - allDomainTrackSize
%less than 500 is good

SPACE_UNITS = 'um'; 
TIME_UNITS = 's';
%%Find each ma object info for each domain%%
Allma=[];

for i=1:numel(AllDomTracks);
    ctracks=AllDomTracks{:,i};
    currentma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS); %dimensionality and Specify units for the msdanalyzer
    currentma = currentma.addAll(ctracks); % Add the sample tracks 
    currentma = currentma.computeDrift('velocity');
    currentma = currentma.computeVCorr;
    currentma = currentma.computeMSD;% Compute MSD for each particle.
    Allma{i}=currentma;
end

save(['30Pa_0_7mA_20domains_ma_Cam1_'],'Allma') % Save msd for specific case into current path

figure
title('Combined Drift over all 12 Domains 46Pa 0.7mA');
ma.plotDrift
AllDomDrift=[];

for domain=1:numel(Allma);
    tma=Allma{domain};
    currentdrift=tma.drift;
    AllDomDrift{domain}=currentdrift(:,2:3);
end

figure
title('Separate Drift of All 12 Domains in 46Pa 0.7mA');

blueShades = [0 0 1; 0 0.5 1; 0 0.8 0.9; 0 0.65 0.98];     % Shades of blue
redShades = [1 0 0; 1 0.1 0.7; 1 0.7 0.7; 1 0.2 0.6];   % Shades of red
greenShades = [0 1 0; 0.1 1 0.7; 0.7 1 0.7; 0.6 0.95 0.6]; % Shades of green

for domain=1:12

    if domain <= 4
        colorSet = blueShades(domain,:);
    elseif domain <= 8
        colorSet = redShades(domain-4,:);
    else
        colorSet = greenShades(domain-8,:);
    end
    
    lineColor = colorSet(rem(domain-1, size(colorSet, 1)) + 1, :);

    hold on
    plot(AllDomDrift{domain}(:, 1), AllDomDrift{domain}(:, 2),'Color',lineColor,'LineWidth',1.5);
end

    
legendLabels = {
    'dom11', 'dom12', 'dom13','dom14', ...
    'dom21', 'dom22', 'dom23','dom24', ...
    'dom31', 'dom32', 'dom33','dom34'
};

% Add a legend
legend(legendLabels);

clf
figure
title('VCorrelation over all 70Pa 0.35mA');
ma.plotMeanVCorr


AllDomVCorr=[];

for domain=1:numel(Allma);
    tma=Allma{domain};
    currentmeanvcorr=tma.getMeanVCorr;
    AllDomMeanVCorr{domain}=currentmeanvcorr;
end

figure
axis([0 20 -2 4])
title('Mean VCorr of All 12 Domains in 46Pa 0.7mA');

blueShades = [0 0 1; 0 0.5 1; 0 0.8 0.9; 0 0.65 0.98];     % Shades of blue
redShades = [1 0 0; 1 0.1 0.7; 1 0.7 0.7; 1 0.2 0.6];   % Shades of red
greenShades = [0 1 0; 0.1 1 0.7; 0.7 1 0.7; 0.6 0.95 0.6]; % Shades of green

for domain=1:12
    hold on

    if domain <= 4
        colorSet = blueShades(domain,:);
    elseif domain <= 8
        colorSet = redShades(domain-4,:);
    else
        colorSet = greenShades(domain-8,:);
    end
    
    lineColor = colorSet(rem(domain-1, size(colorSet, 1)) + 1, :);

    plot(AllDomMeanVCorr{domain}(:, 1), AllDomMeanVCorr{domain}(:, 2),'Color',lineColor);
end

legendLabels = {
    'dom11', 'dom12', 'dom13','dom14', ...
    'dom21', 'dom22', 'dom23','dom24', ...
    'dom31', 'dom32', 'dom33','dom34'
};
legend(legendLabels);




%%Subtract drift from all domains
AllDomtracks_nodrift=[];

for k=1:numel(Allma)
    currentma=Allma{k};
    tracks = currentma.tracks;
    Total_drift = currentma.drift;
    roundedTotal_drift=round(Total_drift(:,1),4);
    number_particles=numel(tracks);

    track_nodrift=[];
    
    for i=1:number_particles;
        particle=tracks{i,1};
        time = particle(:,1);
        framecount = numel(particle)/3; %not an actuall timejust a counter for each particle. 
        % It is divided by 3 because numel(particle) counts amount of time
        % data, x data, and y data, thus divide by 3 since the length of all of
        % these is the same
        % Extract data from Total_drift array based on matching values in ith particle time array
        roundedtime = round(time,4);
        matchingRows = ismember(roundedTotal_drift(:, 1), roundedtime);
        drift = Total_drift(matchingRows, 2:3);
        
        particle_nodrift=zeros(framecount,1);
        for j=1:framecount;
            particle_nodrift(j,1:2) = particle(j,2:3)-drift(j,1:2);
        end
        total_particle = [time,particle_nodrift];
        track_nodrift{i,1}=total_particle;
       
                
    end
    
    AllDomtracks_nodrift{k,1}=track_nodrift;
    
end

%Plot a certain amount of tracks of whatever subsection for viewing
%purposes
hold on
figure
axis([0 15000 -1000 4000])

for i=1:12
    currentdomtracks=AllDomtracks_nodrift{i};
    for j = 1:min(100,numel(currentdomtracks))
        hold on
        plot(currentdomtracks{j}(:, 2), currentdomtracks{j}(:, 3));
    end
   
end

%The order of the cells with the domains matter for plotting
cell12=AllDomtracks_nodrift{1};
cell10=AllDomtracks_nodrift{2};
cell9=AllDomtracks_nodrift{3};
cell4=AllDomtracks_nodrift{8};
cell5=AllDomtracks_nodrift{5};
cell6=AllDomtracks_nodrift{6};
cell7=AllDomtracks_nodrift{7};
cell8=AllDomtracks_nodrift{4};
cell3=AllDomtracks_nodrift{9};
cell2=AllDomtracks_nodrift{10};
cell11=AllDomtracks_nodrift{11};
cell1=AllDomtracks_nodrift{12};


allNoDrift=[cell1;cell2;cell3;cell4;cell5;cell6;cell7;cell8;cell9;cell10;cell11;cell12];

%Cell ordering matters for this plotting as improper cell order will stack cells of particles
%on top of other cells when plotting. Unclear why and what the correct
%pattern is however the combined tracks from the cells will calculate an
%independent MSD ma object which has no prior knowledge of this cell order
%and is therefore fine
figure
axis([0 15000 -1000 4000])
hold on

N = 2000; % Set the desired number of tracks to see. for all set N=numel(track_nodrift)
M = numel(allNoDrift);

for i = 1:min(N, numel(allNoDrift))
    plot(allNoDrift{i}(:, 2), allNoDrift{i}(:, 3));
end
title(['Nonhomogeneous Drift subtraction 70Pa 07mA ' num2str(N) 'Tracks']);

NoDriftTracks = allNoDrift;
save([num2str(numel(allNoDrift)) '_C1_46Pa_07mA_msd_data_NHDS'],'NoDriftTracks')
