%%% PART 1 %%%

[audio, Fs] = audioread("voice.wav");

% Plot the audio signal
t = (0:length(audio)-1) / Fs;
figure;
plot(t, audio);
xlabel('Time (s)');
ylabel('Amplitude');
title('Audio Signal');

% Double the sampling rate
Fs_2 = Fs * 2;

% Save the audio with the doubled sampling rate
audiowrite("./out_doubled.wav", audio, Fs_2);


%%% PART 2 %%%

% Generate decreasing exponential signal
n = length(audio);
exp_signal = exp(-0.00005 * (0:n-1))';

% Multiply and plot the result
modified_audio = audio .* exp_signal;

% Plot the result
figure;
plot(t, modified_audio);
xlabel('Time (s)');
ylabel('Amplitude');
title('Modified Audio Signal');

% Save the audio signal
audiowrite('modified_audio.wav', modified_audio, Fs);


%%% PART 3 %%%

first_echo_delay = 1.0;
second_echo_delay = 2.0;

first_echo_intensity = 0.8;
second_echo_intensity = 0.5;

% Create the echo filter
d1 = round(first_echo_delay * Fs);
d2 = round(second_echo_delay * Fs);
echo_filter = [1; zeros(d1-1, 1); first_echo_intensity; zeros(d2-d1-1, 1); second_echo_intensity];

% Convolve with the audio for each channel separately
echoed_audio = zeros(size(audio));
for ch = 1:size(audio, 2)
    echoed_audio(:, ch) = conv(audio(:, ch), echo_filter, 'same');
end

% Save the result
audiowrite('echoed_audio.wav', echoed_audio, Fs);


%%% PART 4 %%%

% Load concert_hall_IR.wav
[ir, Fs_ir] = audioread("concert_hall_IR.wav");

% Convert to mono if necessary
if size(ir, 2) > 1
    ir = mean(ir, 2);
end

% Ensure the sampling rates match
if Fs ~= Fs_ir
    error('Sampling rates do not match!');
end

% Plot the impulse response
t_ir = (0:length(ir)-1) / Fs_ir;
figure;
plot(t_ir, ir);
xlabel('Time (s)');
ylabel('Amplitude');
title('Impulse Response of Concert Hall');

% Convolve with the audio and save the result
synthesized_audio_Hall = zeros(size(audio));
for ch = 1:size(audio, 2)
    synthesized_audio_Hall(:, ch) = conv(audio(:, ch), ir, 'same');
end
% Normalize and save the result
synthesized_audio_Hall = synthesized_audio_Hall / max(abs(synthesized_audio_Hall(:)));
audiowrite('synthesized_audio_concert_hall.wav', synthesized_audio_Hall, Fs);

% Plot the synthesized_audio_bucket
t_syn = (0:length(synthesized_audio_Hall)-1) / Fs;
figure;
plot(t_syn, synthesized_audio_Hall);
xlabel('Time (s)');
ylabel('Amplitude');
title('Synthesized Audio with Concert Hall IR');

% Repeat for iron_bucket_IR.wav
[ir_bucket, Fs_ir_bucket] = audioread("iron_bucket_IR.wav");

% Convert to mono if necessary
if size(ir_bucket, 2) > 1
    ir_bucket = mean(ir_bucket, 2);
end

% Ensure the sampling rates match
if Fs ~= Fs_ir_bucket
    error('Sampling rates do not match!');
end

% Plot the impulse response
t_ir_bucket = (0:length(ir_bucket)-1) / Fs_ir_bucket;
figure;
plot(t_ir_bucket, ir_bucket);
xlabel('Time (s)');
ylabel('Amplitude');
title('Impulse Response of Iron Bucket');

% Convolve with the audio and save the result
synthesized_audio_bucket = zeros(size(audio));
for ch = 1:size(audio, 2)
    synthesized_audio_bucket(:, ch) = conv(audio(:, ch), ir_bucket, 'same');
end
% Normalize and save the result
synthesized_audio_bucket = synthesized_audio_bucket / max(abs(synthesized_audio_bucket(:)));
audiowrite('synthesized_audio_iron_bucket.wav', synthesized_audio_bucket, Fs);

% Plot the synthesized_audio_bucket
t_syn_bucket = (0:length(synthesized_audio_bucket)-1) / Fs;
figure;
plot(t_syn_bucket, synthesized_audio_bucket);
xlabel('Time (s)');
ylabel('Amplitude');
title('Synthesized Audio with Iron Bucket IR');



%%% PART 5 %%%

% Parameters
T = 6; % Period
k_values = -500:500; % Harmonic range
N = 1000; % Number of points for numerical integration
t = linspace(-3, 3, N); % Time vector

% Define the function x(t) numerically
x_t = -1 .* (-3 <= t & t < -1) + ...
       1 .* (-1 <= t & t < 1) + ...
      -1 .* (1 <= t & t <= 3);

% Numerical integration using trapezoidal rule to calculate Fourier coefficients
a_k = zeros(size(k_values));
for idx = 1:length(k_values)
    k = k_values(idx);
    integrand = x_t .* exp(-1j * 2 * pi * k * t / T);
    a_k(idx) = (1 / T) * trapezoidal_rule(@(x) interp1(t, integrand, x, 'linear'), -3, 3, N);
end

% Plot the Fourier series for different ranges of harmonics
ranges = [-20 20; -100 100; -500 500];
figure;
for i = 1:size(ranges, 1)
    subplot(3, 1, i);
    k_range = ranges(i, 1):ranges(i, 2);
    stem(k_range, real(a_k(k_range + 501)));
    xlabel('k');
    ylabel('a_k');
    title(['Fourier Series Coefficients for k = [', num2str(ranges(i, 1)), ' ', num2str(ranges(i, 2)), ']']);
end

% New x(t) for the second part
x_t_new =  1 .* (-3 <= t & t < -2) + ...
          -1 .* (-2 <= t & t < 0) + ...
           1 .* (0 <= t & t < 2) + ...
          -1 .* (2 <= t & t <= 3);

% Numerical integration to calculate new Fourier series coefficients
c_k = zeros(size(k_values));
for idx = 1:length(k_values)
    k = k_values(idx);
    integrand = x_t_new .* exp(-1j * 2 * pi * k * t / T);
    c_k(idx) = (1 / T) * trapezoidal_rule(@(x) interp1(t, integrand, x, 'linear'), -3, 3, N);
end

% Verify the relation between a_k and c_k for a_3 and c_3
a_3 = a_k(k_values == 3);
c_3 = c_k(k_values == 3);
disp(['a_3 = ', num2str(a_3)]);
disp(['c_3 = ', num2str(c_3)]);

% Plot the new Fourier series coefficients
figure;
stem(k_values, real(c_k));
xlabel('k');
ylabel('c_k');
title('Fourier Series Coefficients for New x(t)');


%%% PART 6 %%%

% Load the audio file
[original_audio, Fs] = audioread("original_audio.wav");

% Define parameters for 360° audio
duration = length(original_audio) / Fs;
t = linspace(0, duration, length(original_audio))';

% Create sinusoidal modulation for phase shift
frequency = 1; % Frequency of the 360° effect
left_channel = original_audio(:, 1) .* sin(2 * pi * frequency * t);
right_channel = original_audio(:, 2) .* cos(2 * pi * frequency * t);

% Combine the channels
audio_360 = [left_channel, right_channel];

% Normalize and save the 360° audio file
audio_360 = audio_360 / max(abs(audio_360(:)));
audiowrite('audio_360.wav', audio_360, Fs);