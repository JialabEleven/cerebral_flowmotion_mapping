% LSPIV with no parallel processing enabled.
% 
% This file is to illustrate how the technique works, showing the gaussian 
% peak fitting in real time during calculation.  The speed of the technique
% is limited by this non-linear fitting step.  The speed of the technique can
% be sped up through parallelization (our approach -- see LSPIV_parallel.m) 
% or a more efficient algorithm.  
%
% For additional information, please see corresponding paper:
%
% 'Line-Scanning Particle Image Velocimetry: an Optical Approach for
% Quantifying a Wide Range of Blood Flow Speeds in Live Animals'
% by Tyson N. Kim, Patrick W. Goodwill, Yeni Chen, Steven M. Conolly, Chris
% B. Schaffer, Dorian Liepmann, Rong A. Wang
%
% PWG 3/28/2012
% 出来的数据是负数的可能原因是扫描频率不够导致的，这时候我们需要降低设置的shiftamt(如果为1还是算出负数，那只能再提高扫描频率重新选择数据了)
close all
clc
clear all
% Parameters to improve fits
maxGaussWidth = 100;  % maximum width of peak during peak fitting

% Judge correctness of fit
numstd        = 3;  %num of stdard deviation from the mean before flagging
windowsize    = 1024; %in # scans, this will be converted to velocity points
                      %if one scan is 1/2600 s, then windowsize=2600 means
                      %a 1 second moving window.  Choose the window size
                      %according to experiment.
                      %这里与扫描频率有关系，895Hz就选择895
distance_pixel=0.994;
total_time=124.999;
cur_p1 = mfilename('fullpath');%获得现在所打开的文件的路径
i=strfind(cur_p1,'\');%匹配 \
cur_p1 = mfilename('fullpath');
[cur_p,~,~] = fileparts(cur_p1);

buildingDir = fullfile(cur_p, 'model.tif');
imgData1 = imread(buildingDir);

savepath = cur_p1;
%%  settings
% Ask user for setting 
str = {'capillary', 'artery', 'user'};
[speedSetting,v] = listdlg('PromptString','Select a configuration',...
                'SelectionMode','single',...
                'ListString',str);
if v == 0; beep; disp('Cancelled'); %break; 
end

if speedSetting == 1   % CAPILLARY SETTING
    numavgs       = 100;  %up to 100 (or more) for noisy or slow data
    skipamt       = 10;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
    shiftamt      = 2; %shiftamt代表跳过的帧数，如果capillary(~1mm/s),当扫描频率为900Hz时，我认为使用1-2比较合适
elseif speedSetting == 2   % ARTERY SETTING
    numavgs       = 100;  %up to 100 (or more) for noisy or slow data
    skipamt       = 20;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
    shiftamt      = 1;
elseif speedSetting == 3   % USER SETTING
    disp('settings are hard coded in the script, see script!');
    numavgs       = 100;  %up to 100 (or more) for noisy or slow data
    skipamt       = 20;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
    shiftamt      = 2;
end


%% Import the data from a multi-frame tif and make into a single array
%  The goal is a file format that is one single array,
%  so modify this section to accomodate your raw data format.
%
%  This particular file format assumes

% disp('import raw data');
% [fname,pathname]=uigetfile('*.TIF','pick a linescan file');%loads file
% if fname == 0; beep; disp('Cancelled'); %break; 
% end
% imageLines = imimportTif([pathname fname])';
%filename = 'G:\20240624_new_papaer_stroke\Venule_arteriole_wwt\Laser960-Laser4%-25x-P9-line1-T10000(arteriole)(side).png.frames\9_right_side.png'; %path of picture
%filename = 'G:\20240624_new_papaer_stroke\Venule_arteriole_wwt\Laser960-Laser20%-25x-P3_line1-T10000(vein)(8X zoom in)(center).png.frames\center_1.png'; %path of picture
picture = imread(buildingDir);
img_double=double(picture(:,:,1));
imageLines = img_double / 255.0;
figure(),imshow(imageLines);

%% choose where in the image to process
imagesc(imageLines(1:size(imageLines,2),:,1))
colormap('gray')

title('Select the boundaries of the region of interest 1/2');
[X1,Y1] = ginput(1);
line([X1 X1],[1 size(imageLines,2)]);

title('Select the boundaries of the region of interest 2/2');
[X2,Y2] = ginput(1);
line([X2 X2],[1 size(imageLines,2)]);
refresh
pause(.01);

startColumn   = round(min(X1, X2));      % Defines what part of the image we perform LSPIV on.
endColumn     = round(max(X1, X2));

%% minus out background signal (PWG 6/4/2009)
tic
disp('DC correction')
DCoffset = mean(imageLines, 1);
imageLinesDC = bsxfun(@minus, imageLines, DCoffset);

%% do LSPIV correlation
disp('LSPIV begin');

scene_fft  = fft(imageLinesDC(1:end-shiftamt,:),[],2);
test_img   = zeros(size(scene_fft));
test_img(:,startColumn:endColumn)   = imageLinesDC(shiftamt+1:end, startColumn:endColumn);
test_fft   = fft(test_img,[],2);
W      = 1./sqrt(abs(scene_fft)) ./ sqrt(abs(test_fft)); % phase only

LSPIVresultFFT      = scene_fft .* conj(test_fft) .* W; 
LSPIVresult         = ifft(LSPIVresultFFT,[],2);
disp('LSPIV complete');
    
%% find shift amounts
disp('Find the peaks');
velocity = [];
maxpxlshift = round(size(imageLines,2)/2)-1;

index_vals = skipamt:skipamt:(size(LSPIVresult,1) - numavgs);
numpixels = size(LSPIVresult,2);
velocity  = nan(size(index_vals));
amps      = nan(size(index_vals));
sigmas    = nan(size(index_vals));
goodness  = nan(size(index_vals));

%% iterate through
close all
for index = 1:length(index_vals)
    
    if mod(index_vals(index),100) == 0
        fprintf('line: %d\n',index_vals(index))
    end
    
    LSPIVresultAVG   = fftshift(sum(LSPIVresult(index_vals(index):index_vals(index)+numavgs,:),1)) ...
                                      / max(sum(LSPIVresult(index_vals(index):index_vals(index)+numavgs,:),1));
    
    % find a good guess for the center
    c = zeros(1, numpixels);
    c(numpixels/2-maxpxlshift:numpixels/2+maxpxlshift) = ...
        LSPIVresultAVG(numpixels/2-maxpxlshift:numpixels/2+maxpxlshift);
    [maxval, maxindex] = max(c);
    
    % fit a guassian to the xcorrelation to get a subpixel shift
    options = fitoptions('gauss1');
    options.Lower      = [0    numpixels/2-maxpxlshift   0            0];
    options.Upper      = [1e9  numpixels/2+maxpxlshift  maxGaussWidth 1];
    options.StartPoint = [1 maxindex 10 .1];
        
    [q,good] = fit((1:length(LSPIVresultAVG))',LSPIVresultAVG','a1*exp(-((x-b1)/c1)^2) + d1',options);
    
    %save the data
    velocity(index)  = (q.b1 - size(LSPIVresult,2)/2 - 1)/shiftamt;
    amps(index)      = q.a1;
    sigmas(index)    = q.c1;
    goodness(index)  = good.rsquare;

    % plot the data out
    figure(1)
    subplot(1,2,1)
    hold off
    plot(LSPIVresultAVG);
    hold all
    plot(q);
    xlim([1 numpixels]);
    ylim([0 1]);
    
    
    subplot(1,2,2)
    plot(index_vals,velocity)
    title('displacement [pixels/scan]');
    ylabel('pixel shift');
        
    xlabel('scan number');        
    pause(.01);
end
%% find possible bad fits
toc

% Find bad velocity points using a moving window
% 这里如果设置的numstd为2，那么是velocity大于局部velocity平均值的2倍标准差为异常值并记录
pixel_windowsize = round(windowsize / skipamt);

badpixels = zeros(size(velocity));
for index = 1:1:length(velocity)-pixel_windowsize
    pmean = mean(velocity(index:index+pixel_windowsize-1)); %partial window mean
    pstd  = std(velocity(index:index+pixel_windowsize-1));  %partial std 
    
    pbadpts = find((velocity(index:index+pixel_windowsize-1) > pmean + pstd*numstd) | ...
                   (velocity(index:index+pixel_windowsize-1) < pmean - pstd*numstd));%2倍标准差被认为是异常点

    badpixels(index+pbadpts-1) = badpixels(index+pbadpts-1) + 1; %running sum of bad pts
end
badvals  = find(badpixels > 0); % turn pixels into indicies
goodvals = find(badpixels == 0);

meanvel  = mean(velocity(goodvals)); %overall mean
stdvel   = std(velocity(goodvals));  %overall std


velocity=abs(velocity*distance_pixel*(1/1000)/(total_time/size(img_double,1)));%单位为mm/s
%velocity=abs(velocity*(distance_pixel*(size(img_double,1)/1000)/total_time));
velocity_revise= velocity';
time_interval=index_vals'*(total_time/size(img_double,1));
result=[time_interval velocity_revise];

speed = result(:, 2);
mean_speed = mean(speed);
std_speed = std(speed);
velocity_threshold = 100;

% 初始化记录
outlier_count = 0;
outlier_positions = [];
outlier_values = [];


% 计算 Z-score
z_scores_velocity = (speed - mean_speed) / std_speed;

% 修正异常值(中位数作为修正值，也可以将平均值作为修正值)
for i = 1:length(speed)
    if abs(z_scores_velocity(i)) > velocity_threshold
        outlier_count=outlier_count+1;%记录异常值数量
        outlier_positions(end+1)=i;%记录异常值数量
        outlier_values(end+1)=speed(i);%记录异常值的值

        % 计算前后50帧的范围
        start_idx = max(1, i-50); % 起始索引
        end_idx = min(length(speed), i+50); % 结束索引
        
        % 获取前后帧的中位数
        speed(i) = median(speed(start_idx:end_idx));

    end
end

result(:, 2) = speed;

% % 标记异常值
% idx = abs(speed - mean_speed) <= threshold;
% 
% % 去除异常值
% data_cleaned = data(idx, :);

% total_time=13.299;
% distance_pixel=0.124;
% velocity=abs(velocity*(distance_pixel*10/total_time));
% velocity_revise= velocity';
% time_interval=index_vals'*(total_time/10000);
% result=[time_interval velocity_revise];


% show results
figure(2)
subplot(3,1,1)
imgtmp = zeros([size(imageLines(:,startColumn:endColumn),2) size(imageLines(:,startColumn:endColumn),1) 3]); % to enable BW and color simultaneously
imgtmp(:,:,1) = imageLines(:,startColumn:endColumn)'; imgtmp(:,:,2) = imageLines(:,startColumn:endColumn)'; imgtmp(:,:,3) = imageLines(:,startColumn:endColumn)';
imagesc(imgtmp/max(max(max(imgtmp))))
title('Raw Data');
ylabel('[pixels]');
%colormap('gray');

subplot(3,1,2)
imagesc(index_vals,-numpixels/2:numpixels/2,fftshift(LSPIVresult(:,:),2)');
title('LSPIV xcorr');
ylabel({'displacement'; '[pixels/scan]'});


subplot(3,1,3)
plot(index_vals, velocity,'.');
hold all
plot(index_vals(badvals), velocity(badvals), 'ro');
hold off
xlim([index_vals(1) index_vals(end)]);
ylim([meanvel-stdvel*4 meanvel+stdvel*4]);
title('Fitted Pixel Displacement');
ylabel({'displacement'; '[mm/s]'});
xlabel('index [pixel]');
saveas(gca,[savepath,'merge_1.tif']);


h = line([index_vals(1) index_vals(end)], [meanvel meanvel]);
set(h, 'LineStyle','--','Color','k');
h = line([index_vals(1) index_vals(end)], [meanvel+stdvel meanvel+stdvel]);
set(h, 'LineStyle','--','Color',[.5 .5 .5]);
h = line([index_vals(1) index_vals(end)], [meanvel-stdvel meanvel-stdvel]);
set(h, 'LineStyle','--','Color',[.5 .5 .5]);
fprintf('\nMean  Velocity %0.2f [pixels/scan]\n', meanvel);
fprintf('Stdev Velocity %0.2f [pixels/scan]\n', stdvel);

figure()
plot(result(:, 1),result(:, 2));
title('velocity');
ylabel({'velocity'; '[mm/s]'});
saveas(gca,[savepath,'Velocity_1.tif']);

% 输出异常值信息
fprintf('检测到的异常值数量: %d\n', outlier_count);
fprintf('异常值的位置和对应值:\n');
for i = 1:outlier_count
    fprintf('位置: %d, 时间: %f, 异常速度: %f\n', outlier_positions(i), result(outlier_positions(i), 1), outlier_values(i));
end


% 假设result是你的数据，其中第一列是时间，第二列是速度
time = result(:, 1);
speed_1 = result(:, 2);

% 初始化
segment_size = 10;  % 每个分段包含10行数据
num_segments = ceil(length(speed_1) / segment_size);  % 向上取整计算总的分段数量
SNR_values = zeros(num_segments, 1);  % 用于存储每段的SNR值
time_avg = zeros(num_segments, 1);  % 用于存储每段的平均时间

for j = 1:num_segments
    % 计算每个分段的索引范围
    idx_start = (j-1) * segment_size + 1;
    idx_end = min(j * segment_size, length(speed_1));  % 确保索引不超出数据长度
    
    % 计算该分段的平均速度和标准差
    segment_speed = speed_1(idx_start:idx_end);
    mean_speed = mean(segment_speed);
    std_speed = std(segment_speed);
    
    % 计算SNR
    SNR_values(j) = mean_speed / std_speed;
    
    % 计算该分段的平均时间
    time_avg(j) = mean(time(idx_start:idx_end));
end

z_scores = (SNR_values - mean(SNR_values)) / std(SNR_values);

% 设定阈值，通常可以选择Z-score大于3的值为离群值
outlier_threshold = 3;
outliers_idx = abs(z_scores) > outlier_threshold;

% 找出异常值的索引（坐标）
outliers_idx1 = find(abs(z_scores) > outlier_threshold);

% 保存异常值的坐标和对应的 SNR 值
outlier_values = SNR_values(outliers_idx1);

% 显示异常值及其坐标
disp('异常值的坐标（索引）及对应的 SNR 值:');
for p = 1:length(outliers_idx1)
    fprintf('坐标: %d, SNR值: %f\n', outliers_idx1(p), outlier_values(p));
end


% 使用相邻点的平均值来替换离群值
for n = 1:length(SNR_values)
    if outliers_idx(n)
        if n == 1
            SNR_values(n) = SNR_values(n+1);  % 如果第一个点是离群值，用下一个点的值替换
        elseif n == length(SNR_values)
            SNR_values(n) = SNR_values(n-1);  % 如果最后一个点是离群值，用前一个点的值替换
        else
            SNR_values(n) = mean(SNR_values([n-1, n+1]));  % 否则，用前后两个点的平均值替换
        end
    end
end

output_SNR_time=[time_avg SNR_values];
output_result=result;
% 绘制SNR变换曲线
figure;
plot(time_avg, SNR_values, '-o');
title('SNR Variation Over Time');
xlabel('Time');
ylabel('SNR');
grid on;
saveas(gca,[savepath,'SNR_1.tif']);