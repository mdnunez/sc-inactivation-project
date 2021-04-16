function SCMap(RFx, RFy, rfindatapre, rfindatapost,filename)
%This script is to plot the peak velocity map in cartesian/polar
%coordinates and in the SC Map. Snippets of this code were taken from
%plot_data.m, written by Christian Quaia (Quaia et al, 1998) to plot the
%visual field in polar coordinates and logarithmic warping for the SC Map
%and a script for the SC grid


%% double check that the pre and the post data have the same targets

if sum(unique([rfindatapost.amplitude]) == unique([rfindatapre.amplitude])) ~= numel(unique([rfindatapre.amplitude]))
    display('-------------------ERROR: PRE AND POST DO NOT HAVE THE SAME AMPLITUDES--------------------');
end

if sum(unique([rfindatapost.angle]) == unique([rfindatapre.angle])) ~= numel(unique([rfindatapre.angle]))
    display('-------------------ERROR: PRE AND POST DO NOT HAVE THE SAME POLAR DEGREES--------------------');
end

%% Getting the site of injection (lesion) in angle and amplitude
lesr = sqrt((RFx/10)^2 + (RFy/10)^2);
lest = (atan2(RFy/10, RFx/10)) * 180/pi;
if lest < 0 %if angle is negative
    lest = lest + 360;
end 

%% Making another layer of data for the 360 degree which is the 0 degree to complete the circle
x = numel(rfindatapre) + 1;
y = numel(rfindatapost) + 1;
for r = 1:numel(rfindatapre)
    if rfindatapre(r).angle == 0   %Making redundant data points for 0 degree datapoints to have 360 degree datapoints in order to complete the interpolation cirlce
        rfindatapre(x) = rfindatapre(r);
        rfindatapre(x).angle = 360;
        x = x+1;
    end
    if rfindatapre(r).angle == 45    %Making redundant data points for 45 degree datapoints to have 405 degree datapoints in order to complete the interpolation cirlce
        rfindatapre(x) = rfindatapre(r);
        rfindatapre(x).angle = 405;
        x = x+1;
    end

    
    if rfindatapost(r).angle == 0
        rfindatapost(y) = rfindatapost(r);
        rfindatapost(y).angle = 360;
        y = y+1;
    end
    if rfindatapost(r).angle == 45
        rfindatapost(y) = rfindatapost(r);
        rfindatapost(y).angle = 405;
        y = y+1;
    end
end

%% Getting the difference in peak velocity between pre and post
DiffStruct = struct('DiffAvgPV', [], 'NormalizedDiffPV',[], 'angles', [],'amplitude', []);

for pre_r = 1:numel(rfindatapre)
    currentpolarpre = rfindatapre(pre_r).angle;
    currentamppre = rfindatapre(pre_r).amplitude;
    
    post_r = find([rfindatapost.angle] == currentpolarpre & [rfindatapost.amplitude] == currentamppre);
    
    Diff = rfindatapost(post_r).AveragePV - rfindatapre(pre_r).AveragePV;
    DiffStruct(pre_r).angles  = rfindatapre(pre_r).angle;
    DiffStruct(pre_r).amplitude = rfindatapre(pre_r).amplitude;
    DiffStruct(pre_r).DiffAvgPV = Diff;
    DiffStruct(pre_r).NormalizedDiffPV = (Diff/rfindatapre(pre_r).AveragePV)*100;
end
%% Interpolation -- do scatteredinterpolant linear method on the data
angles_q = 0:1:450;
amplitude_q = ceil(min([DiffStruct.amplitude])):1:ceil(max([DiffStruct.amplitude]));

[Xq, Yq] = meshgrid(angles_q,amplitude_q);  %datapoints for the interpolation; X for angle, Y for amplitude

FDiffPV_norm = scatteredInterpolant([DiffStruct.angles]',[DiffStruct.amplitude]', [DiffStruct.NormalizedDiffPV]', 'linear', 'none');
interp_DiffPVmatrix_norm = FDiffPV_norm(Xq,Yq);

%convert polar angles/degrees into x and y cartesian points for plotting
targ_h = Yq .* cos(pi/180*Xq);
targ_v = Yq .* sin(pi/180*Xq);
les_h = lesr * cos(pi/180*lest);
les_v = lesr * sin(pi/180*lest);

%% Plot polar map
figure
cmfixn = -60; % Color scale minimum
cmfixx = 60; % Color scale maximum

map = newmap(cmfixn,cmfixx,0,0);	% Color map % (max,min,bw,level)

s = pcolor(targ_h,targ_v,interp_DiffPVmatrix_norm); %plot by pcolor
% s.LineStyle = 'none';
caxis([-60 60]);
shading interp;
hold on;

h = [DiffStruct.amplitude] .* cos( (pi/180) * [DiffStruct.angles]);
v = [DiffStruct.amplitude] .* sin( (pi/180) * [DiffStruct.angles]);
scatter(h(:), v(:), 20, 'filled', 'k')
ylabel('Vertical Degrees');
xlabel('Horizontal Degrees');
title(['% Difference in peak velocity -- ' filename])
clear l1;
scatter(les_h, les_v, 200, 'x', 'w', 'linewidth', 2);

%plot the angle/amplitude grid 
amplitude_grid = [0,2.5, 5, 10, 20, ceil(max([DiffStruct.amplitude]))-1];
angle_grid = [0:30:360];
[gridangle, gridamplitude] = meshgrid(angle_grid, amplitude_grid);
grid_h = gridamplitude .* cos(pi/180*gridangle);
grid_v = gridamplitude .* sin(pi/180*gridangle);
plot(grid_h, grid_v, 'k', 'linewidth', 1)

for a = 1:numel(amplitude_grid)
[ang_cir, amp_cir] = meshgrid(0:1:360, amplitude_grid(a));
h_cir = amp_cir .* cos(pi/180*ang_cir);
v_cir = amp_cir .* sin(pi/180*ang_cir);
plot(h_cir, v_cir, 'k', 'linewidth', 1);
end 

%% Setting parameters for plotting logarithmic warping for SC

overlap = 1;  %set it to the interval at which the interpolation query points were set 
%some parameters for logarithmic warping
Bx = 1.4;
By = 1.8;
Amm = 3.0;

siz = size(Xq); %number of amplitudes (rows) and directions (col) from interpolated data

%creating empty matrices
left_t = [];
left_r = [];
left_d = [];

right_t = [];
right_r = [];
right_d = [];


for x=1:siz(2)  %for each direction/angle

    %if the datapoints are on the right side of the visual field, then
    %you're going to be drawing the left SC
    
    %theres a 1 degree overlap
    if (Xq(1,x) < (90+overlap)) | ((Xq(1,x) > (270-overlap)) & (Xq(1,x) < (450+overlap)))
        left_t = [left_t Xq(:,x)];  %all the angles on the right side of the visual field
        left_r = [left_r Yq(:,x)];  %all the amplitudes
        left_d = [left_d interp_DiffPVmatrix_norm(:,x)];  %all the data on the right side of the visual field
    end
    
    %same for the left visual field --> right SC data points
    if (Xq(1,x) > (90-overlap)) & (Xq(1,x) < (270+overlap))
        right_t = [right_t Xq(:,x)];
        right_r = [right_r Yq(:,x)];
        right_d = [right_d interp_DiffPVmatrix_norm(:,x)];
    end;
    
end;

%if lesion/RF is on the right side of the visual field, then it's on
%the left SC map
if (lest < (90+overlap)) | ((lest > (270-overlap)) & (lest < (450+overlap)))
    left_lesr = lesr;
    left_lest = lest;
end;

%if lesion/RF is on the left side of the visual field, then it's on
%the right SC map
if (lest > (90-overlap)) & (lest < (270+overlap))
    right_lesr = lesr;
    right_lest = lest;
end;

%convert these datapoints into cartesian for plotting
left_h = left_r .* cos(pi/180*left_t);
left_v = left_r .* sin(pi/180*left_t);
right_h = right_r .* cos(pi/180*right_t);
right_v = right_r .* sin(pi/180*right_t);

% Add 3 to these cartesian points (Amm = 3.0) for horizontal
left_h = left_h + Amm; %for the left SC datapoints (h is positive), add 3
right_h = (-1) * right_h + Amm; %for the right SC datapoints (h is negative), make it positive and add 3

%logarithmic warping
left_x = (Bx / 2) * log((left_h .^ 2 + left_v .^ 2) / (Amm * Amm)); % x is the output of the log function of horizontal and vertical plotting positions with Bx (1.4) and Amm (3) parameter
left_y = By * atan2(left_v, left_h);  %y is By parameter multiplied by the polar angle in radians of horiz and vert plotting
right_x = (Bx / 2) * log((right_h .^ 2 + right_v .^ 2) / (Amm * Amm));
right_y = By * atan2(right_v, right_h);

%make all left_x negative
left_x = (-1) * left_x;

[v,d]=max(max(left_x));  %max of the data point matrix for left side for x (least negative)
siz = size(left_x); %size of the data point matrix for left side for x

if (siz(1,2) > 1),
    if (left_t(1,d) < (90+overlap)) & ((left_t(1,d+1) > (270-overlap)) & (left_t(1,d+1) < (450+overlap)))
        d = d + 1; %increment d
    end;
end;
left_t = [left_t(:,d:siz(1,2)) left_t(:,1:(d-1))];
left_r = [left_r(:,d:siz(1,2)) left_r(:,1:(d-1))];
left_x = [left_x(:,d:siz(1,2)) left_x(:,1:(d-1))];
left_y = [left_y(:,d:siz(1,2)) left_y(:,1:(d-1))];
left_d = [left_d(:,d:siz(1,2)) left_d(:,1:(d-1))];

%getting the actual data point locations (non-interpolated) and putting
%them to their respective sides
for x=1:numel(DiffStruct),
    if (DiffStruct(x).angles < (90+overlap)) | (DiffStruct(x).angles> (270-overlap)) & (DiffStruct(x).angles < (450+overlap))
        left_angle(x) = DiffStruct(x).angles;
        left_amplitude(x) = DiffStruct(x).amplitude;
    end;
    
    if (DiffStruct(x).angles > (90-overlap)) & DiffStruct(x).angles < (270+overlap)
        right_angle(x) = DiffStruct(x).angles;
        right_amplitude(x) = DiffStruct(x).amplitude;
    end;
end;
%logarithmically warping the target positions to put it in the SC map format, like it
%was done with the interpolated points
left_horiz_point = left_amplitude .* cos(pi/180*left_angle);
left_vert_point = left_amplitude .* sin(pi/180*left_angle);
right_horiz_point = right_amplitude .* cos(pi/180*right_angle);
right_vert_point = right_amplitude .* sin(pi/180*right_angle);

left_horiz_point = left_horiz_point + Amm;
right_horiz_point = (-1) * right_horiz_point + Amm;

left_x_point = (Bx / 2) * log((left_horiz_point .^ 2 + left_vert_point .^ 2) / (Amm * Amm)); %(Bx / 2) * log((left_horiz_point .^ 2 + left_vert_point .^ 2) / (Amm * Amm));
left_y_point = By * atan2(left_vert_point, left_horiz_point);
right_x_point = (Bx / 2) * log((right_horiz_point .^ 2 + right_vert_point .^ 2) / (Amm * Amm));
right_y_point = By * atan2(right_vert_point, right_horiz_point);
left_x_point = (-1) * left_x_point;

if RFx > 0  %if RF is postive, then RF is on the right visual field, then it's the left SC
    l1h = left_lesr * cos(pi/180*left_lest);
    l1v = left_lesr * sin(pi/180*left_lest);
    l1h = l1h + Amm;
    left_lesx = (Bx / 2) * log((l1h ^ 2 + l1v ^ 2) / (Amm * Amm));
    left_lesy = By * atan2(l1v, l1h);
    left_lesx = (-1) * left_lesx;
    clear l1h l1v;
end;

if RFx < 0 %if RF is negative, then RF is on the left visual field, then it's the right SC
    r1h = right_lesr * cos(pi/180*right_lest);
    r1v = right_lesr * sin(pi/180*right_lest);
    r1h = (-1) * r1h + Amm;
    right_lesx = (Bx / 2) * log((r1h ^ 2 + r1v ^ 2) / (Amm * Amm));
    right_lesy = By * atan2(r1v, r1h);
    clear r1h r1v;
end;

%%
figure(7)
s = pcolor(left_x,left_y,left_d); %plotting left SC first
scgrid1(-max([DiffStruct.amplitude]), 'k'); %plotting the SC grid

%put in datapoints
hold on;
pnt = plot(left_x_point, left_y_point,'k.');
set(pnt,'MarkerSize',20);
hold off;
hold on;
pcolor(right_x,right_y,right_d);
scgrid1(max([DiffStruct.amplitude]), 'k');

hold on;
pnt = plot(right_x_point(:), right_y_point(:),'k.');
set(pnt,'MarkerSize',20);
hold off;
caxis([-60 60]);
axis normal;
xlim([-3.5 3.5]);

      xlabel('mm SC');
    ylabel('mm SC');
    ax = get(gcf,'CurrentAxes');
    set(ax,'XTick',[-4:4],'YTick',[-3:3]);
    shading interp;
    title(['% Difference in peak velocity -- ' filename])
    
    if RFx < 0 %if RFx is postive, which means RF is on right visual field which means left SC
    hold on;
    clear l1;
    scatter(right_lesx,right_lesy,200,'xw','linewidth', 2);

    hold off;
    
    if RFx > 0 %if RFx is postive, which means RF is on right visual field which means left SC
    hold on;
    clear l1;
    l1 = scatter(left_lesx,left_lesy,200,'xw','linewidth', 2);
    hold off;
end
end
end 
