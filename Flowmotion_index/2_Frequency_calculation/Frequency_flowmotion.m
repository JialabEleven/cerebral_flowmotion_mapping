clc
clear all

% 获取当前脚本所在文件夹
script_path = mfilename('fullpath');
savepath = fileparts(script_path);

% Excel 文件路径
buildingDir = fullfile(savepath, 'Velocity.xlsx');

% 读取 Excel 数据
data = readmatrix(buildingDir);

% 第一列为 time
time = data(:,1);

% 初始化存储原始频谱数据的变量
all_original_spectrum = [];
all_frequencies = [];


for i = 1:(size(data,2)-1)

    amplitude = data(:,i+1);

    figure();
    plot(time, amplitude);

    response = amplitude;
    original_mean = mean(response);
    response = response - original_mean;

    % 添加零填充来减少边界效应
    zero_pad = 2000;
    response_padded = [zeros(zero_pad,1); response; zeros(zero_pad,1)];

    y = fft(response_padded);

    % define sampling characteristic
    sampleNumber = length(response);
    sampleFrequency = 1/time(1,1);

    % 计算频率轴
    delta_f = sampleFrequency / sampleNumber;
    f = (-sampleNumber/2 : sampleNumber/2-1) * delta_f;

    f_frequency = f(size(amplitude,1)/2+1:size(amplitude,1));

    % 计算未滤波的原始频谱
    X_original = fftshift(abs(fft(response))) / sampleNumber;
    original_spectrum = X_original(sampleNumber/2+2:end);
    frequencies = f(sampleNumber/2+2:end);

    % 保存频率和对应的原始频谱值
    if i == 1
        all_frequencies = frequencies;
    end

    all_original_spectrum = [all_original_spectrum, original_spectrum];

    frequency_add = [0 frequencies];
    original_spectrum_add = [0 original_spectrum'];

    % 原始频谱图绘制
    figure();
    X_original = fftshift(abs(y))/sampleNumber;
    plot(frequency_add, original_spectrum_add);
    xlim([0 20]);
    title(['原始频谱图 - ', num2str(i)]);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');

    saveas(gcf, fullfile(savepath, ...
        ['original_spectrum_', num2str(i), '.tif']));


    % ==========================
    % 低通滤波：0-0.2 Hz
    % ==========================
    wp_low = 0.2/(sampleFrequency/2);
    ws_low = 0.3/(sampleFrequency/2);

    [N_low, wc_low] = buttord(wp_low, ws_low, 3, 10);
    [b_low, a_low] = butter(N_low, wc_low, 'low');

    filtered_velocity_low = filter(b_low, a_low, response_padded);
    filtered_velocity_low = filtered_velocity_low(zero_pad+1:end-zero_pad);


    % ==========================
    % 带通滤波：0.6-2 Hz
    % ==========================
    wp_band = [0.6 2]/(sampleFrequency/2);
    ws_band = [0.4 2.5]/(sampleFrequency/2);

    alpha_p_band = 6;
    alpha_s_band = 9;

    [N_band, wc_band] = buttord( ...
        wp_band, ws_band, alpha_p_band, alpha_s_band);

    [b_band, a_band] = butter(N_band, wc_band);

    filtered_velocity_band = filter( ...
        b_band, a_band, response_padded);

    filtered_velocity_band = ...
        filtered_velocity_band(zero_pad+1:end-zero_pad);


    % ==========================
    % 带通滤波：5-9 Hz
    % ==========================
    wp_high = [5 9]/(sampleFrequency/2);
    ws_high = [4.5 9.5]/(sampleFrequency/2);

    [N_high, wc_high] = buttord( ...
        wp_high, ws_high, 6, 9);

    [b_high, a_high] = butter( ...
        N_high, wc_high, 'bandpass');

    filtered_velocity_high = filter( ...
        b_high, a_high, response_padded);

    filtered_velocity_high = ...
        filtered_velocity_high(zero_pad+1:end-zero_pad);


    filtered_velocity_low = ...
        filtered_velocity_low + original_mean;

    filtered_velocity_high = ...
        filtered_velocity_high + original_mean;

    filtered_velocity_band = ...
        filtered_velocity_band + original_mean;


    % ==========================
    % 绘制低通时域图
    % ==========================
    figure();

    t = (0:length(filtered_velocity_low)-1) ...
        * (1/sampleFrequency);

    plot(t, filtered_velocity_low);
    grid on;

    title(['低通滤波时域图形 - ', num2str(i)]);
    xlabel('Time (s)');
    ylabel('Velocity');

    saveas(gcf, fullfile(savepath, ...
        ['time_domain_low_', num2str(i), '.tif']));


    % ==========================
    % 绘制低通频域图
    % ==========================
    figure();

    X_low = fftshift(abs(fft(filtered_velocity_low))) ...
        / sampleNumber;

    filtered_data_low = ...
        X_low(sampleNumber/2+1:sampleNumber);

    f_add = [0 f(sampleNumber/2+2:sampleNumber)];

    filtered_data_low_add = ...
        [0 filtered_data_low(2:end,:)'];

    plot(f_add, filtered_data_low_add);

    xlim([0 0.3]);

    title('低通滤波频域幅度特性(0-0.3Hz)');

    saveas(gcf, fullfile(savepath, ...
        ['frequency_domain_low_', num2str(i), '.tif']));


    % ==========================
    % 绘制高通时域图
    % ==========================
    figure();

    plot(t, filtered_velocity_high);
    grid on;

    title(['高通滤波时域图形 - ', num2str(i)]);
    xlabel('Time (s)');
    ylabel('Velocity');

    saveas(gcf, fullfile(savepath, ...
        ['time_domain_high_', num2str(i), '.tif']));


    % ==========================
    % 绘制高通频域图
    % ==========================
    figure();

    X_high = fftshift(abs(fft(filtered_velocity_high))) ...
        / sampleNumber;

    filtered_data_high = ...
        X_high(sampleNumber/2+1:sampleNumber);

    plot( ...
        f(sampleNumber/2+2:sampleNumber), ...
        filtered_data_high(2:end,:) ...
        );

    xlim([5 14]);

    title('高通滤波频域幅度特性 (7-14Hz)');

    saveas(gcf, fullfile(savepath, ...
        ['frequency_domain_high_', num2str(i), '.tif']));


    % ==========================
    % 绘制带通时域图
    % ==========================
    t = (0:length(filtered_velocity_band)-1) ...
        * (1/sampleFrequency);

    figure();

    plot(t, filtered_velocity_band);
    grid on;

    title(['改进后的带通滤波时域图形 - ', num2str(i)]);
    xlabel('Time (s)');
    ylabel('Velocity');

    saveas(gcf, fullfile(savepath, ...
        ['time_domain_band_', num2str(i), '.tif']));


    % ==========================
    % 绘制带通频域图
    % ==========================
    figure();

    X_band = fftshift(abs(fft(filtered_velocity_band))) ...
        / sampleNumber;

    filtered_data_band = [ ...
        X_band(size(amplitude,1)/2+1: ...
        size(amplitude,1)/2+1); ...
        2*X_band(size(amplitude,1)/2+2: ...
        size(amplitude,1))];

    plot( ...
        f(size(amplitude,1)/2+2:size(amplitude,1)), ...
        filtered_data_band(2:end,:) ...
        );

    xlim([0.5 3]);

    title('带通滤波频域幅度特性 (0.5-2.5Hz)');

    saveas(gcf, fullfile(savepath, ...
        ['frequency_domain_band_', num2str(i), '.tif']));


    % ==========================
    % 汇总时域信号
    % ==========================
    band_velocity_low(:,i) = filtered_velocity_low;
    band_velocity_high(:,i) = filtered_velocity_high;
    band_velocity(:,i) = filtered_velocity_band;

    filtered_data_low_merge(:,i) = filtered_data_low;
    filtered_data_band_merge(:,i) = filtered_data_band;
    filtered_data_high_merge(:,i) = filtered_data_high;

end


% =============================================
% 保存到 Excel
% =============================================

n = size(band_velocity_high, 2);

column_names = {'Frequency'};

for j = 1:n
    column_names{end+1} = ['Low' num2str(j)];
end

for j = 1:n
    column_names{end+1} = ['Band' num2str(j)];
end

for j = 1:n
    column_names{end+1} = ['High' num2str(j)];
end


% ==========================
% 保存时域滤波结果
% ==========================
output_velocity_Data = [ ...
    time ...
    band_velocity_low ...
    band_velocity ...
    band_velocity_high ...
    ];

output_table_velocity = array2table( ...
    output_velocity_Data, ...
    'VariableNames', column_names);

filename = 'output_velocity_data.xlsx';

fullpath = fullfile(savepath, filename);

writetable( ...
    output_table_velocity, ...
    fullpath, ...
    'Sheet', 'FilteredData');


% ==========================
% 保存修正时间结果
% ==========================
N = height(output_table_velocity);

output_table_velocity_revise = ...
    output_table_velocity(1:N,:);

time_interval = time(1,1);

new_time = ...
    (1:height(output_table_velocity_revise)).' ...
    * time_interval;

output_table_velocity_revise.Frequency = new_time;

filename_revise = ...
    'output_velocity_revise_data.xlsx';

fullpath_revise = ...
    fullfile(savepath, filename_revise);

writetable( ...
    output_table_velocity_revise, ...
    fullpath_revise, ...
    'Sheet', 'FilteredData');


% ==========================
% 保存滤波后频谱
% ==========================
output_frequency_data = [ ...
    f_frequency(:,2:end)' ...
    filtered_data_low_merge(2:end,:) ...
    filtered_data_band_merge(2:end,:) ...
    filtered_data_high_merge(2:end,:) ...
    ];

output_table_frequency = array2table( ...
    output_frequency_data, ...
    'VariableNames', column_names);

filename1 = 'output_frequency_data.xlsx';

fullpath1 = fullfile(savepath, filename1);

writetable( ...
    output_table_frequency, ...
    fullpath1, ...
    'Sheet', 'FilteredData');


% ==========================
% 保存原始频谱
% ==========================
output_frequency_data_raw = [ ...
    all_frequencies', ...
    all_original_spectrum ...
    ];

column_names1 = {'Frequency'};

for jj = 1:n
    column_names1{end+1} = ['raw_' num2str(jj)];
end

output_table_frequency_raw = array2table( ...
    output_frequency_data_raw, ...
    'VariableNames', column_names1);

new_filename = 'original_spectrum_data.xlsx';

fullpath_new = ...
    fullfile(savepath, new_filename);

writetable( ...
    output_table_frequency_raw, ...
    fullpath_new, ...
    'Sheet', 'OriginalSpectrum');