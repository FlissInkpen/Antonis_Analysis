function [output_analsysis_data]=Antonis_new_script_ver_1
clear all
close all
clc
%Antonis analysis code
%Written by Felicity Inkpen, 21st May 2019

%% Information to start

%Each folder is a specific day with a specific animal.
%There are multiple folders for different animals and different days.

%Code to take the sorted data, clustered using KlustaKwick, to analyse a
%bunch of parameters, including firing rate, spatial information, amplitude
%of the spikes, etc.

%Calculate how stable the firing is between different session
%Correlate activity over time and spatial bins
%Create and plot figures
%Export clusters from the spiking data, export a .csv file, to include all
%the parameters and metadata

%8 tetrodes in total
%Code should be transferrable for different number of tetrodes

%The code asks for some starting inputs:
% How many sessions you had (could be taken from the .tseg file returned
% from the KlustaKwik analysis)
% If there are no inputs, the script will run the analysis for all
% tetrodes, additional inputs allow you to specify which tetrode is of
% interest.

%Raw data files of all data formats have the number of the sessions run as
%the last digit

%the clustering script needs to be run first; this is a separate script.
%This script analyses the data already preprocessed by KlustaKwik (+manual curration).


%KlustaKwick uses principal component analysis on the raw spiking data,
%from all the data formats.
%It returns a file: merge.tseg, which contains the time segments of the
%sessions.  This information is necessary to separate out the different sessions in later analysis.
%each session should be about 10 minutes, but not precisely.

%Ultimately create an output excel file of all the rats

%The experimental protocol takes a snapshot of the data once the voltage
%passes a certain threshold on one of the four wires of the tetrode.

%% Basic Inputs to the function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the figure defaults
startup
%% Mannual inputs
% Manual inputs through dialogue box. This asks for four variables, identified through experimental setup and manual measurement:
% The pixel ratio, the number of pixels per metre. The secondary pixel ratio, assumed to be not applicable.
% The number of tetrodes used in experiment. The number of sessions.
%defaults of 297, 8, NaN and 6, respectively, are given
% user inputs the combinations of paradigms and animals (uigetfile, % questdlg, etc)
%long term goal to automate, not a priority right now.
directory = '\\mvm-sbms-130383.bms.ed.ac.uk\frax\Fliss\Antonis\';           % define the directory where all your data and analysis code is stored

[data_folder, date, pixel_ratio, pixel_ratio_2, tetrodes, arenas, sessions]...    % retrieve your starting information
    =starting_information(directory);
cd(data_folder);                                                            % change directory to the folder containing the relevant data
disp(['Cluanalysis will now run on: ', pwd]);                               % Shows which directory Matlab is running in
disp('-----------------------------------------------------------------');


%session_output makes a cell to collect the information from each session:
%use filenames, electrodes, clusters, etc. 

%want to make a structure of cells for each session
%this wants to later be populated with all the information from all the
%respective sessions. 
session_output=[];
for i=1:sessions
    session_output.(['session_' (num2str(i))])={};
end

%% Find and load the set files
[Set_Files]=set_files;

%add the set files to the relevant sub structures of the sessions structure
for i=1:sessions
session_output.(['session_' num2str(i)]){1,1}=Set_Files.(['file_' num2str(date) num2str(i)]);

end

%% Find the value of the gain at each electrode in each session. 
for i=1:sessions
    gain_mat=gain_analysis(i, session_output);
    session_output.(['session_' num2str(i)]){2,1}=gain_mat;
end

%% Check that you have the correct number of arenas and sessions
if arenas>1
    if exist('merge.goal')>0
        load('merge.goal')
        Arenas=merge;
        if length(Arenas)~=arenas
            error('Please check the number of arenas used in your experiment and run the code again')
        else
            fprintf('Analysis done on data from %d different arenas.  ', max(Arenas))
        end
    else
        error('No goal file present, please run the script again and check the number of arenas')
    end
else
    fprintf('Analysis done on data from %d arena.  ', arenas)
end



tic;                                                                        % meaure how long it takes the function to run
electrodes = length(dir('*.clu*'));                                         %define the number of electrodes from the number of items in the folder

%% Inputs from question dialogue boxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Is there a sleep session?
sleep=sleep_question;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the time and space coordinates
[merge_coords, time_stamp]=animal_coordinates;                              % function to find the time and space coordinates, and plot these


% want to be able to look at a specific tetrode, and a specific cell.,
%These would want to be deifned in the inputs of the function, with
%defaults that ran all tetrodesand all cells.






%limits of the field coordinates
limsX= [min(merge_coords(:,1)), max(merge_coords(:,1))];                    % the maximum and minimum values of the x-column of the position matrix
limsY= [min(merge_coords(:,2)), max(merge_coords(:,2))];                    % the minimum and maximum values of the y column of the position matrix

%variables for making a smoothed 2 dimensional histograph / placefield
%diagram
sigma = 15;                                                                                       % sigma (gaussian standard deviation) to be used for rate and position map smoothing
min_dwell = 0.0001;                                                                      % total number of seconds that rat has to be in a bin for it to count
bin_size = 2.5; 												                        % (cm), for calculating the rate map.
min_dwell_distance = 5; 										                                % (cm) the distance from a point to determine minimum dwell time
min_dwell_time = 100; 											                                % (ms) minimum dwell time in ms for plotting rate map (should be in the multiple of 20)
dt_position = 20; 											                                % sampling interval of position data (ms)
smooth = 3; 												                                % smooth factor for guassian smoothing AKA sigma. In each script it is converted from cm to pixels.




%within each electrode, there will be a number of clusters

for i=1:electrodes
    eno=i;                                                                      % electrode number, this will be looped through the number of electrodes
    filename = 'merge';                                                         % the files we are interested in are those that have been processed and merged by KlustaKwik
    spkfile = [filename '.spk.' num2str(eno)];
    fetfile = [filename '.fet.' num2str(eno)];                                  % there are three file types, the spiking data, the features (metadata) and the cluster data
    clufile = [filename '.clu.' num2str(eno)];                                  % define filenames by their extensions
    
    %[filenames, posx, posy, hd, post, timeseg, seg_size] = read_mypos(filename);
    %Load the appropriate files
    cd(data_folder)                                                             %change directory to the data folder
    
    %Import the spiking data
    fid = fopen(spkfile);                                                       % create a file identifier - an fid.
    if (fid == -1)                                                              % if the file cannot be opened
        disp(sprintf('\tERROR: File not found: %s',spkfile));                   % display an error message
    else
        waves = fread(fid,Inf,'int16=>int8');                                   % otherwise, read the file to create a variable called waves
        waves = reshape(waves, 4, 50, []);                                      % reshape this
    end
    fclose(fid);                                                                % close the file identifier
    
    %Import the feature metadata
    fid = fopen(fetfile);                                                       % define the file identifier as the feature file variable
    if (fid == -1)                                                              % if this filetype can't be opened
        disp(sprintf('\tERROR: File not found: %s',fetfile));                   % display an error message
    else
        nfeatures = sscanf(fgetl(fid),'%d');                                    % Otherwise, define the number of features by reading the file identifier string as formatted data
        features = fscanf(fid,'%d',[nfeatures Inf]);                            % define the features by scanning the file identifier; there are 24 features for each
        fclose(fid);                                                            % close the file identifier
    end
    fet_per_tet = (nfeatures -4) / 4;                                           %calculate the features per tetrode - the number of features analysed across the four tetrodes, minus the time and position coordinates.
    % Typically there are 24 features, minus the four for time and position, divided by four for the four tetrodes = 5 features per tetrode.
    num_spikes = length(features(1,:));                                         %The number of spike is equal to the number of columns of data
    disp(sprintf('\t...done'));                                                 % display 'done' - indented in the command window
    
    
    %Import the corresponding clusters data
    clusters = load(clufile);                                                   % Load the corresponding clusters data
    
    % define the cluster index
    nclu= max(clusters);                                                        % the number of clusters is defined by the largest number in the clu file
    temp=[];                                                                    % create a temp matrix
    for k=2:nclu                                                      % each data point is asigned a cluster, (1, 2, 3,4,..etc) find the number of clusters
        temp(k)=length(clusters(clusters==k));                                  % define the temp variable as the length of the longest cluster
    end
    index=NaN(max(temp),nclu);                                         % make an empty matrix for the index, with the dimensions of the length of the longest cluster, and the width of the number of clusters
    for j=2:nclu                                                      % go through each of the clusters
        temp=(find(clusters==j));                                               % find the indicies with that cluster number
        index(1:length(temp),j)=temp;                                           % in the index matrix, input the indicies into the corresponding column
    end
    figure
    for j=2:nclu                                                     % for 1: the number of clusters (1,2,3,4... etc)
        index_temp=index(:,j);                                                  % make a temporary version of the index, with just the indicies from that cluster
        index_temp(isnan(index_temp))=[];
        if max(index_temp)>length(features)                                     % if you index is longer than your content, remove the last entry in the index
            index_temp=index_temp(1:end-1);                                     % im not sure why this has happened, so this is a quick bug fix
        end
        features_temp=features(:,index_temp);                                   % make a temporary version of the features, with only the features with the corresponding index
        waves_temp=waves(:,:,index_temp);                                       % likewise, make a temporary version of the wave data, with only the features from the corresponding index
        fprintf('There are %d contributions to cluster %d from electrode %d. ', size(features_temp,2),j, i);             % display the sizes of each identified cluster
        %Whatever cluster analysis needs to be done should go here.
        
        %The positions recorded in the feature matrix
        n_features=size(features_temp,1);
        hds = features_temp(nfeatures-3,:);                                         % head direction information - irrelevant for Antonis
        posxs = features_temp(nfeatures-2,:);                                       % The x coordinates from the feature matrix
        posys = features_temp(nfeatures-1,:);                                       % The y coordinates from the feature matrix
        %remove any position data that is outside the preset limits
        posxs(find(posxs<limsX(1)|posxs>limsX(2)))=NaN;                             % remove x coordinates that are outside of the maximum and minimum limits defined
        posys(find(posys<limsY(1)|posys>limsY(2)))=NaN;                             % remove x coordinates that are outside of the maximum and minimum limits defined
        
        %Scale the bins to the pixel ratio
        smooth_scaled = smooth /100 * pixel_ratio;                                  % scale the smoothing factor to the pixel ratio
        binsize = bin_size /100 * pixel_ratio;                                      % scale the bin size to the pixel ratio
        min_dwell_dist = min_dwell_distance / 100 * pixel_ratio;                    % scale the min dwell distance to the pixel ratio
        
        limsx=([min(features_temp((size(features_temp,1)-2),:)),...
            max(features_temp((size(features_temp,1)-2),:))]);
        limsy=([min(features_temp((size(features_temp,1)-1),:)),...
            max(features_temp((size(features_temp,1)-1),:))]);
        % Define the number of bins
        xbins=ceil((limsx(2)-limsx(1))/binsize);
        ybins=ceil((limsy(2)-limsy(1))/binsize);
        
        test_mat_structure=zeros(ybins,xbins);
        
        for k=1:length(features_temp)
            x_coord=round((features_temp((size(features_temp,1)-2),k)-limsx(1))./((limsx(2)-limsx(1))/xbins));
            
            y_coord=round((features_temp((size(features_temp,1)-1),k)-limsy(1))./((limsy(2)-limsy(1))/ybins));
            if x_coord==0
                x_coord=1;
            end
            if y_coord ==0
                y_coord =1;
            end
            test_mat_structure(y_coord,x_coord)=test_mat_structure(y_coord,x_coord)+1;
        end
        
        output_structure=interp2( test_mat_structure, smooth);
        subplot(2, ceil((nclu-1)/2),j-1)
        contourf(output_structure,20, 'LineColor','none');
        colormap('jet');
        c=colorbar;
        c.Label.String='Spiking';
        title(['Place Field Map, Cluster ',num2str(j)])
        
        
    end
    suptitle(['Electrode ' num2str(i)])
    suplabel('Arena X-coordinate (cm)','x');
    suplabel('Arena Y-coordinate (cm)','y');
    
end

%

%% Spiking as a function of position to generate place fields

%The clustered feature data contains information about each spike that was
%fired.
%This data is only collected in the time period surrounding the spike, to
%ensure computation power is sufficient
% - you dont want to record all the time or that would take up all your RAM
% Therefore, the presence of each column of the feature matrix corresponds
% to one spike.
%We want to take the position of each spike, and create a place field
%diagram.
%this will be smoothed with a gaussian filter, the working parameters of
%which are outlined in the previous code.

%inputs - the feature matrix for each session
% this needs to be organised by electrode
% You are trying to locate a place field, find out how stable it is, how
% big it is, etc.

%outputs -a smoothed placefield matrix for each cluster across each session
%figures of the place fields, in usual conventions

%the crude position coordinates are taken at a sampling frequency of 50Hz,
%i.e one measurement every 20ms.

% the spiking data and metadata have a much higher sampling rate,
% which only kicks in when the threshold for spiking is passed,
%  this means that they just have a time stamp

%ultimately, what you want is a place field for each electrode and each
%cluster
%if you have 6 different sessions, you want 6 different place fields for
%each cluster
%if you have <100 spikes per session, the cluster will exist, but there
%wont be enough information to calculate the place field.

%the function that performs this analysis currently is called
%pfield_s3_posT


% 




%% Outputs,
%Master excel file, to contain inputs, in the form of .csv files, for each
%rat, to include all data and metadata

% Create data folders for later and save inputs %%
% CreateFolders
% disp(sprintf('\t...saving Clua_inputs'));
% save('Inputs/Clua_inputs.mat','date','elecs','clus','colourbar_set')
% disp('------------------------------------');


%Desired figures: both separate figures, saved as .svg and .png, and
%compiled data types and figures in master figures, saved as png.




toc

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   N E S T E D   I N P U T    F U N C T I O N S     %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function startup
%Figure defaults
%axes
set(0, 'DefaultTextFontSize', 14);
set(0, 'DefaultAxesFontSize', 16);
set(0, 'DefaultUicontrolFontSize', 8);
set(0,'DefaultAxesTickDir', 'out')
set(0,'DefaultAxesLineWidth', 1.5)

%data
set(0,'DefaultLineLineWidth',1)
co = [0    0    1
    1    0    0
    0.3  0.3  0.3
    0.2    1  0.2
    0      1    1
    0      0    0];
set(0,'DefaultAxesColorOrder',co)
set(0,'DefaultAxesAmbientLightColor',[1 1 1])
set(0, 'DefaultAxesBox', 'off')
end
function [data_folder, date, pixel_ratio, pixel_ratio_2, tetrodes, arenas, sessions]=starting_information(directory)
prompt = {'Date of experiment (year, month, day)',...                       % define a prompt in the form of an input window
    'Number of pixels per metre:','Secondary pixel ratio',...               % include all the necessary variables 
    'Number of tetrodes:','Number of differently shaped arenas:', ...
    'Number of recording sessions'};
definput = {'20160516','297','NaN','8', '1','6'};                           % set the default inputs
answer = inputdlg(prompt,'Inputs',[1 44],definput);                         % the output returned is a cell structure containing the answers
%define the individual outputs
date= str2double(answer{1,1}); pixel_ratio=str2double(answer{2,1});         % define numerically all the variables from the output cell structure
pixel_ratio_2=str2double(answer{3,1}); tetrodes=str2double(answer{4,1});
arenas = str2double(answer{5,1}); sessions=str2double(answer{6,1});

data_folder = ([directory  num2str(date)]);                                 %define the data folder

end
function CreateFolders                                                                                                    % by Roddy %%                                                                                                                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       This script makes folders for images etc in the current directory
%       Any folder names contained in fold will be created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Making data folders...')
fold = {'Figures','Data','Inputs'};                                                                                                     % a list of folders to make
num = 0;                                                                                                                                % start a counter
for file = 1:length(fold)                                                                                                               % for every folder name listed in fold
    name = fold{file};                                                                                                              % get the filename
    if ~exist(name, 'dir')	                                                                                                        % if a folder with that name does not exist in the current directory
        disp(sprintf('\t...making %s folder',name));                                                                 % display a message
        mkdir(name);                                                                                                            % make a new folder with the current folder name
        num = num + 1;	                                                                                                        % increment counter
    end % if exist(a, 'file')
end % for file = 1:length(ext)
disp(sprintf('\t...done'));
end
function [ax,h]=suplabel(text,whichLabel,supAxes)
% PLaces text as a title, xlabel, or ylabel on a group of subplots.
% Returns a handle to the label and a handle to the axis.
%  [ax,h]=suplabel(text,whichLabel,supAxes)
% returns handles to both the axis and the label.
%  ax=suplabel(text,whichLabel,supAxes)
% returns a handle to the axis only.
%  suplabel(text) with one input argument assumes whichLabel='x'
%
% whichLabel is any of 'x', 'y', 'yy', or 't', specifying whether the
% text is to be the xlable, ylabel, right side y-label,
% or title respectively.
%
% supAxes is an optional argument specifying the Position of the
%  "super" axes surrounding the subplots.
%  supAxes defaults to [.08 .08 .84 .84]
%  specify supAxes if labels get chopped or overlay subplots
%
% EXAMPLE:
%  subplot(2,2,1);ylabel('ylabel1');title('title1')
%  subplot(2,2,2);ylabel('ylabel2');title('title2')
%  subplot(2,2,3);ylabel('ylabel3');xlabel('xlabel3')
%  subplot(2,2,4);ylabel('ylabel4');xlabel('xlabel4')
%  [ax1,h1]=suplabel('super X label');
%  [ax2,h2]=suplabel('super Y label','y');
%  [ax3,h2]=suplabel('super Y label (right)','yy');
%  [ax4,h3]=suplabel('super Title'  ,'t');
%  set(h3,'FontSize',30)
%
% SEE ALSO: text, title, xlabel, ylabel, zlabel, subplot,
%           suptitle (Matlab Central)
% Author: Ben Barrowes <barrowes@alum.mit.edu>
%modified 3/16/2010 by IJW to make axis behavior re "zoom" on exit same as
%at beginning. Requires adding tag to the invisible axes
%modified 8/8/2018 to allow cells as text for multiline capability
currax=findobj(gcf,'type','axes','-not','tag','suplabel');
if nargin < 3
    supAxes=[.08 .08 .84 .84];
    ah=findall(gcf,'type','axes');
    if ~isempty(ah)
        supAxes=[inf,inf,0,0];
        leftMin=inf;  bottomMin=inf;  leftMax=0;  bottomMax=0;
        axBuf=.04;
        set(ah,'units','normalized')
        ah=findall(gcf,'type','axes');
        for ii=1:length(ah)
            if strcmp(get(ah(ii),'Visible'),'on')
                thisPos=get(ah(ii),'Position');
                leftMin=min(leftMin,thisPos(1));
                bottomMin=min(bottomMin,thisPos(2));
                leftMax=max(leftMax,thisPos(1)+thisPos(3));
                bottomMax=max(bottomMax,thisPos(2)+thisPos(4));
            end
        end
        supAxes=[leftMin-axBuf,bottomMin-axBuf,leftMax-leftMin+axBuf*2,bottomMax-bottomMin+axBuf*2];
    end
end
if nargin < 2, whichLabel = 'x';  end
if nargin < 1, help(mfilename); return; end
if (~isstr(text) & ~iscellstr(text)) | ~isstr(whichLabel)
    error('text and whichLabel must be strings')
end
whichLabel=lower(whichLabel);
ax=axes('Units','Normal','Position',supAxes,'Visible','off','tag','suplabel');
if strcmp('t',whichLabel)
    set(get(ax,'Title'),'Visible','on')
    title(text);
elseif strcmp('x',whichLabel)
    set(get(ax,'XLabel'),'Visible','on')
    xlabel(text);
elseif strcmp('y',whichLabel)
    set(get(ax,'YLabel'),'Visible','on')
    ylabel(text);
elseif strcmp('yy',whichLabel)
    set(get(ax,'YLabel'),'Visible','on')
    ylabel(text);
    set(ax,'YAxisLocation','right')
end
%for k=1:length(currax), axes(currax(k));end % restore all other axes
for k=1:length(currax), set(gcf,'CurrentAxes',currax(k));end % restore all other axes
if (nargout < 2)
    return
end
if strcmp('t',whichLabel)
    h=get(ax,'Title');
    set(h,'VerticalAlignment','middle')
elseif strcmp('x',whichLabel)
    h=get(ax,'XLabel');
elseif strcmp('y',whichLabel) | strcmp('yy',whichLabel)
    h=get(ax,'YLabel');
end
end
function [merge_coords, time_stamp]=animal_coordinates(sessions)
merge_coords=load('merge.mypos');                                           % load the coordinates (in x,y,z and time), as a matrix, merge_coords

% remove any coordinates that are outside of the arena.
four_stdevs=4.*std(merge_coords);                                           % find any data that is four standard deviations away from the mean
acceptable_range=mean(merge_coords)-four_stdevs;                            % define an acceptable range to be anything larger than four standard deviation less than the mean
for i=1:length(merge_coords)                                                % go through the coordinates
    for j=1:2
        if merge_coords(i,j)<acceptable_range(1,1:2)                        % is the x or y coordinates are smaller than the limit of the acceptable range
            merge_coords(i,j)=NaN;                                          % replace those with not-a-number
        end
    end
end


%Create a figure of the rat's trajectory over time
figure
colour_map=(round((merge_coords(:,4)./60),3));                              % create a colour vector from the time component of the merge.position data
scatter(merge_coords(:,1),merge_coords(:,2),8, colour_map,'filled');        % scatter the rat's position, with colour showing time
h=colorbar; h.Limits=[0 60];                                                % specify the limits of the colourbar
h.Label.String = 'Time (minutes)';                                          % label the colourbar as 'time'
xlabel('Arena X-coordinate (cm)'); ylabel('Arena Y-coordinate (cm)');       % label the axes
hold on


arena = questdlg('Define arena shape', ...                                  % Create a question dialogue box
    'arena', 'circular','square','other','circular');                       % Default is circular arena
if strcmp(arena, 'circular')==1
    centre_coords=((max(merge_coords(:,1:2))-min(merge_coords(:,1:2)))./2)...   %define the centre of the arena
        +min(merge_coords(:,1:2));
    x=centre_coords(1,1); y=centre_coords(1,2);                                 %define the coordinates of the centre of the arena
    r= max((max(merge_coords(:,1:2))-min(merge_coords(:,1:2)))./2)+1;           %define the radius of the arena
    th = 0:pi/50:2*pi;                                                          % define theta (angle)
    xunit = r * cos(th) + x;                                                    % define the factors of a circl
    yunit = r * sin(th) + y;
    plot(xunit,yunit,'Color',[0.5, 0.5, 0.5],'LineWidth',6);                   % plot a circle around the data, to indicate the boundary of the arena
else
    disp('arena is square or otherwise');                                       % Otherwise, display 'no sleep data to analyse
end




%load the time stamps
%The time segments (timestamps) merge.tseg
time_stamp=load('merge.tseg');                                              % Load the time segments file
if sessions~=length(time_stamp)
    fprintf('Please examine the time stamp file. \nThe number of time segments does not match up to the number of sessions.\n')
else
    return
end

sessions=length(time_stamp);                                                % the number of sessions is the length of the time stamps file
session_length=length(merge_coords)/sessions;                               % the session length is the length of the coordinates file, divided by the number of the sessions
individual_session=merge_coords;
figure
for i =1:sessions
    session_1=individual_session(1:session_length,:);
    subplot(3, 2, i);
    colour_map=(round((session_1(:,4)./60),3));                              % create a colour vector from the time component of the merge.position data
    scatter(session_1(:,1),session_1(:,2),8, colour_map,'filled');
    hold on
    if strcmp(arena, 'circular')==1
        plot(xunit,yunit,'Color',[0.5, 0.5, 0.5],'LineWidth',6);                   % plot a circle around the data, to indicate the boundary of the arena
    else
        disp('arena is square or otherwise');                                       % Otherwise, display 'no sleep data to analyse
    end
    
    h=colorbar;                                              % specify the limits of the colourbar
    h.Label.String = 'Time (minutes)';
    individual_session=individual_session(session_length+1:end,:);
end
suplabel('Arena X-coordinate (cm)','x');
suplabel('Arena Y-coordinate (cm)','y');

end
% timeseg is the start time of each trial
% post_all is the times the locations were sampled
function [filenames, posx_all, posy_all, hd_all, post_all, timeseg, seg_size] = read_mypos(out_name)                    %%% Steven Huang<s.huang@ed.ac.uk> %%%

posfn = [out_name '.mypos'];
disp(sprintf(['Reading positions from ',posfn]));

fid = fopen(posfn, 'r');
if (fid == -1)
    error(sprintf('\t...file not found: %s',out_name));
end % if (fid == -1)
[pos_all, count] = fscanf(fid,'%f', [4 inf]);
fclose(fid);

posx_all = pos_all(1,:)';
posy_all = pos_all(2,:)';
hd_all = pos_all(3,:)';
post_all = pos_all(4,:)';
disp(sprintf('\t...done'));                                                                                                             % display 'done' - indented in the command window

segfn = [out_name '.tseg'];
disp(sprintf('Reading time segment information from %s...',segfn));

fid = fopen(segfn,'r');
if (fid == -1)
    error(sprintf('\t...failed to open the time segment file.'));
end % if (fid == -1)
s = textscan(fid, '%s %s');
fclose(fid);
timeseg = str2double(s{2});
filenames = s{1};
time_index = zeros(1, length(timeseg));
seg_size = zeros(1, length(timeseg));
save('temp','post_all','timeseg');
find(post_all==timeseg(1));
for ii=1:length(time_index)
    if (size(find(post_all==timeseg(ii)),1)==0)
        error(sprintf('\tUhoh, the start time (%0.2f) of trial %i is outside the recording times.',timeseg(ii),ii));
        [temp,index] = min(abs(post_all-timeseg(ii)));
        error(sprintf('\tNearest recording time is.at index %i at time %0.2f.',index,post_all(index)));
        error(sprintf('\tThe recording times before and after are at: %0.2f %0.2f.',post_all(index-1),post_all(index+1)));
        error('\tSetting time segment value equal to nearest recorded time sampled.');
        timeseg(ii) = post_all(index);
    end % if (size(find(post_all==timeseg(ii)),1)==0)
    time_index(ii) = find(post_all==timeseg(ii));
end % for ii=1:length(time_index)

for ii=1:length(seg_size)
    if(ii<length(seg_size))
        seg_size(ii) = time_index(ii+1) - time_index(ii);
    else
        seg_size(ii) = length(post_all) - time_index(ii) + 1;
    end % if(ii<length(seg_size))
end % for ii=1:length(seg_size)

end
%Is there a sleep session?
function sleep=sleep_question

%Is there a sleep session?
sleep = questdlg('Is a sleep session included?', ...                        % Create a question dialogue box
    'Sleep', '1 sleep session','Multiple sleep sessions','No','No');        % Default is no sleep sessions
if strcmp(sleep, '1 sleep session')==1
    disp('One Sleep session included, more information needed');            % If there is a sleep session, display this, and perform relevant analysis
elseif strcmp(sleep, 'Multiple sleep sessions')==1
    disp('Multiple Sleep sessions included, more information needed');      % If there are multiple sleep sessions, display this and perform relevant analysis
else
    disp('No sleep data to analyse');                                       % Otherwise, display 'no sleep data to analyse
end
end



function [Set_Files]=set_files
setfiles=dir('*.set');
Set_Files=[];
for i=1:length(setfiles)
    setfile= setfiles(i).name;
    delimiter = ' ';
%% Format for each line of text:
%   column1: text (%q)
%	column2: text (%q)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(setfile,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Create output variable
Set_Files.(['file_' (num2str(setfile(1:9)))]) = [dataArray{1:end-1}];

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

end
    
end



%Function to find the gain applied to each channel in a session 
%need the input of the session_output variable, and the number of the
%session being studied - can be part of a loop. 

function gain_mat=gain_analysis(k, session_output)
                                                                         
temp=strfind(session_output.(['session_' num2str(k)]){1,1}, ...             % create a temporary variable the highlights the parts of the session output data that contains the word 'gain'
    'gain','ForceCellOutput',true);
    gain_index=[];                                                          % create a gain index matrix
    for i=1:length(temp)                                                    % go through the temporary cell structure
        if temp{i,1}==1
            gain_index=[gain_index;i];                                      % add the location of the gain values, one by one
        else
            continue
        end
    end
    
    gain= (session_output.(['session_' num2str(k)]){1,1}(gain_index,2));    % use the previously defined gain index to find the gain at each part of the session output cell structure
    gain_mat=NaN(length(gain),1);                                           % you want a matrix of the different gain values in the form of a double, not a characater, which you can now find from the index
    for i=1:length(gain)                                                    % go through the gain structure
        gain_mat(i,1)=str2double(gain{i,1});                                % do the string to double transformation
    end
    
    p1=bar(gain_mat, 0.8,  'FaceColor', [0.3 0.5 1]);                       % plot a bar graph of the gain values across the channels
    hold on
    p2=plot([32.5 32.5],[0 max(gain_mat)+1000], 'r--', 'LineWidth', 3);     % plot a line to show the point of 32 channels
    ylabel('Gain');    xlabel('Channels');                                  % add axis labels
    xlim([0 length(gain)+1]);                                               % scale the x axis
    legend([p2],{'32 Channels'}); legend boxoff;                            % add a legend to show the point of 32 channels
 
end










