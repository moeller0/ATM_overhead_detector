function [ output_args ] = ATM_overhead_detector( sweep_fqn, up_Kbit, down_Kbit )
%ATM_OVERHEAD_DETECTOR Summary of this function goes here
%   try to read in the result from a ping sweep run
%	sweep_fqn (optional): the log file of the ping sweep against the first hop after
%		the DSL link
%	up_Kbit (optional): the uplink rate in Kilobits per second
%	down_Kbit (optional): the downlink rate in Kilobits per second
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License version 2 as
% published by the Free Software Foundation.
%
%       Copyright (C) 2015 Sebastian Moeller
%
% NOTES:
%	Under octave 3.8.2 under macosx the fltk backend crashes, and the
%	gnuplot backend only exports/saves black boxes... seems to work under
%	linux though
%
% TODO:
%	estimate the best MTU for the estimated protocol stack (how to test this?)
%		1) estimate the largest MTU that avoids fragmentation (default 1500 - 28 should be largest without fragmentation)
%		2) estimate the largest MTU that does not have padding in the last
%		   ATM cell, for this pick the MTU that no partial ATM cell remains
%	include the potential PACKET sizes for VLAN tagged packets as well?
%	try boxcox function to deskew the right-skewed RTT distribution
%
%
%Thoughts:
%	the sweep should be taken directly connected to the modem to reduce
%		non-ATM routing delays

if (ismac)
	octave_use_gnuplot = 1;
else
	octave_use_gnuplot = 1;	% the fltk backend seems to have issues with exporting data via ghostscript
end

if ~(isoctave)
	dbstop if error;
	timestamps.(mfilename).start = tic;
else
	tic();
	if (octave_use_gnuplot)
		graphics_toolkit('gnuplot');
		setenv GNUTERM wxt ; %sm: should be equivalent to setenv("GNUTERM","wxt"); but unlike the later will not confuse matlab
	else
		if (ismac)
			%graphics_toolkit fltk; __init_fltk__ quit
		end
	end
end
disp(['Starting: ', mfilename]);

output_args = [];

% control options
show_mean = 1;					% the means are noisier than the medians
show_robust_mean = 1;			% the mean gets a bit better after excluding the top and bottom 5%
show_median = 1;				% the median seems the way to go
show_min = 1;					% the min should be the best measure, but in the ATM test sweep it is too variable
show_max = 0;					% only useful for debugging
show_sem = 0;					% give some estimate of the variance
show_ci = 1;					% show the confidence interval of the mean, if the mean is shown
show_geomean = 1;				%ATTENTION: this requires the octave-statistics package, installed and loaded (pkg load statistics)
show_robust_geomean = 1;		%ATTENTION: this requires the octave-statistics package
show_delogged_logmean = 0;
ci_alpha = 0.05;				% alpha for confidence interval calculation
use_measure = 'median';			% median, or robust_mean
plot_output_format = 'png';		% what to save
use_processed_results = 1;		% do not parse the ASCII file containg the ping output again (as the parser is very slow)
max_samples_per_size = [];		% if not empty only use maximally that many samples per size
% max_samples_per_size = 1000;	% if not empty only use maximally that many samples per size

if (isoctave)
	if (show_geomean || show_robust_geomean)
		disp('Octave statistics package required.')
		require_octave_stats_pkg = 1;
	end
	
	if exist('require_octave_stats_pkg', 'var') && (require_octave_stats_pkg)
		
		[pkg_is_loadable, pkg_is_loaded] = check_octave_pkg_availability('statistics');
		if (pkg_is_loadable) && ~(pkg_is_loaded)
			disp('Attemptimg to load statistics package');
			pkg load statistics
		end
		% geomean lives in the statistics package so use it as 'canary'
		if ~exist('geomean')
			disp('Could not load statistics package.');
			% change all control parameters that would drag in the statistics package
			show_geomean = 0;
			show_robust_geomean = 0;
			if ismember(use_measure, {'geomean', 'robust_geomean'})
				disp(['Selected analaysis statistic (', use_measure, ') is not available, defaulting to median instead']);
				use_measure = 'median';
			end
		else
			disp('Octave statistics pkg loaded successfully.');
		end
	end
	fflush(stdout); % make octave display intermediate output...
end

% if not specified we try to estimate the per cell RTT from the data
default_up_Kbit = [];
default_down_KBit = [];

if (nargin == 0)
	sweep_fqn = '';
	%  	sweep_fqn = fullfile(pwd, 'ping_sweep_CABLE_20120801_001235.txt');
	if isempty(sweep_fqn)
		[sweep_name, sweep_dir] = uigetfile({'ping*.txt';'ping*.mat'});
		sweep_fqn = fullfile(sweep_dir, sweep_name);
	end
	up_Kbit = default_up_Kbit;
	down_Kbit = default_down_KBit;
end
if (nargin == 1)
	up_Kbit = default_up_Kbit;
	down_Kbit = default_down_KBit;
end
if (nargin == 2)
	down_Kbit = default_down_KBit;
end

%ATM
quantum.byte = 48;	% ATM packets are always 53 bytes, 48 thereof payload
quantum.bit = quantum.byte * 8;
ATM_cell.byte = 53;
ATM_cell.bit = ATM_cell.byte * 8;


% known packet size offsets in bytes
offsets.IPv4 = 20;		% assume no IPv4 options are used, IPv6 would be 40bytes?
offsets.IPv6 = 40;		% not used yet...
offsets.ICMP = 8;		% ICMP header
offsets.ethernet = 14;	% ethernet header
offset.ATM.max_encapsulation_bytes = 44; % see http://ace-host.stuart.id.au/russell/files/tc/tc-atm/, but note that due to VLAN tags we can reach 48 worst case...
MTU = 1500;	% the nominal MTU to the ping host should be 1500, but might be lower if using a VPN
max_MTU_for_overhead_determination = 1280;	% 1280 is true for IPv6, for IPv4 the minMTU is 576
% fragmentation will cause an addition relative large increase in RTT (not necessarily registered to the ATM cells)
% that will confuse the ATM quantisation offset detector, so exclude all
% ping sizes that are potentially affected by fragmentation
max_ping_size_without_fragmentation = MTU + offsets.ethernet - offsets.IPv4 - offset.ATM.max_encapsulation_bytes;
% unknown offsets is what we need to figure out to feed tc-stab...


[sweep_dir, sweep_name] = fileparts(sweep_fqn);
cur_parsed_data_mat = [sweep_fqn(1:end-4), '.mat'];

if (use_processed_results && ~isempty(dir(cur_parsed_data_mat)))
	disp(['Loading processed ping data from ', cur_parsed_data_mat]);
	load(cur_parsed_data_mat, 'ping');
else
	% read in the result from a ping sweep
	disp(['Processing ping data from ', sweep_fqn]);
	ping = parse_ping_output(sweep_fqn);
	if isempty(ping)
		disp('No useable ping data found, exiting...');
		return
	end
	if (isoctave)
		save('-v7', cur_parsed_data_mat, 'ping');
	else
		save(cur_parsed_data_mat, 'ping');
	end
end




% analyze the data
min_ping_size = min(ping.data(:, ping.cols.size)) - offsets.ICMP;
disp(['Minimum size of ping payload used: ', num2str(min_ping_size), ' bytes.']);
known_overhead = offsets.IPv4;	% ping reports the ICMP header already included in size
ping.data(:, ping.cols.size) = ping.data(:, ping.cols.size) + known_overhead;	% we know we used IPv4 so add the 20 bytes already, so that size are relative to the start of the IP header
size_list = unique(ping.data(:, ping.cols.size));	% this is the number of different sizes, but there might be holes/missing sizes
max_pingsize = max(size_list);

% packets larger than the pMTU will get fragmented, resulting in a extra-large step (roughly 2 to 3 times larger than usual) somewhere in the data
% which will confuse the simplistic stair finder, so limit the search space
% to <+ 1280 the min MTU for IPv6, hoping that this should work
% everywhere...
if (size_list(end) > max_MTU_for_overhead_determination)
	disp(['Restricting the ATM quantization search space to <= ', num2str(max_MTU_for_overhead_determination), ' bytes.']);
	tmp_idx = find(size_list <= max_MTU_for_overhead_determination);
	if (isempty(tmp_idx))
		disp(['No data with size <= ', num2str(max_MTU_for_overhead_determination), ' bytes found; ATM quantization can not be determined....']);
		return
	end
	measured_size_list = size_list;
	size_list = measured_size_list(tmp_idx);
	measured_max_pingsize = max_pingsize;
	max_pingsize = max(size_list);
end


per_size.header = {'size', 'mean', 'robust_mean', 'median', 'min', 'max', 'std', 'n', 'sem', 'ci', 'geomean', 'robust_geomean', 'delogged_logmean'};
per_size.cols = get_column_name_indices(per_size.header);
per_size.data = zeros([max_pingsize, length(per_size.header)]) / 0;	% NaNs
per_size.data(:, per_size.cols.size) = (1:1:max_pingsize);

if ~isempty(max_samples_per_size)
	disp(['Analysing only the first ', num2str(max_samples_per_size), ' samples.']);
end

for i_size = 1 : length(size_list)
	cur_size = size_list(i_size);
	
	% throw out negative numbers?
	cur_size_idx = find(ping.data(:, ping.cols.size) == cur_size);
	
	remove_impossible_times = 1;
	if (remove_impossible_times)
		cur_size_n_samples = length(cur_size_idx);
		cur_size_idx = find((ping.data(:, ping.cols.size) == cur_size) & (ping.data(:, ping.cols.time) >= 0));
		if (length(cur_size_idx) < cur_size_n_samples)
			disp(['Excluded ', num2str(cur_size_n_samples - length(cur_size_idx)), ' samples due to negative RTTs (invalid measurements)...']);
		end
	end
	
	if ~isempty(max_samples_per_size)
		n_selected_samples = min([length(cur_size_idx), max_samples_per_size]);
		cur_size_idx = cur_size_idx(1:n_selected_samples);
		%disp(['Analysing only the first ', num2str(max_samples_per_size), ' samples of ', num2str(length(cur_size_idx))]);
	end
	per_size.data(cur_size, per_size.cols.mean) = mean(ping.data(cur_size_idx, ping.cols.time));
	% robust mean, aka mean of 5 to 95 quantiles
	per_size.data(cur_size, per_size.cols.robust_mean) = robust_mean(ping.data(cur_size_idx, ping.cols.time), 0.1, 0.9);	% take the mean while excluding extreme values
	
	per_size.data(cur_size, per_size.cols.median) = median(ping.data(cur_size_idx, ping.cols.time));
	per_size.data(cur_size, per_size.cols.min) = min(ping.data(cur_size_idx, ping.cols.time));
	per_size.data(cur_size, per_size.cols.max) = max(ping.data(cur_size_idx, ping.cols.time));
	per_size.data(cur_size, per_size.cols.std) = std(ping.data(cur_size_idx, ping.cols.time), 0);
	per_size.data(cur_size, per_size.cols.n) = length(cur_size_idx);
	per_size.data(cur_size, per_size.cols.sem) = per_size.data(cur_size, per_size.cols.std) / sqrt(length(cur_size_idx));
	per_size.data(cur_size, per_size.cols.ci) = calc_cihw(per_size.data(cur_size, per_size.cols.std), per_size.data(cur_size, per_size.cols.n), ci_alpha);
	if (show_geomean)
		per_size.data(cur_size, per_size.cols.geomean) = geomean(ping.data(cur_size_idx, ping.cols.time));
	end
	if (show_robust_geomean)
		per_size.data(cur_size, per_size.cols.robust_geomean) = robust_geomean(ping.data(cur_size_idx, ping.cols.time), 0.1, 0.9);	% take the geomean while excluding extreme values
	end
	%per_size.data(cur_size, per_size.cols.delogged_logmean) = 10^(mean(log10(ping.data(cur_size_idx, ping.cols.time))));
	per_size.data(cur_size, per_size.cols.delogged_logmean) = exp(mean(log(ping.data(cur_size_idx, ping.cols.time))));
	
	
	
end

clear ping	% with large data sets 32bit matlab will run into memory issues...


data_fh = figure('Name', sweep_name);
hold on;
legend_str = {};
if (show_mean)
	% means
	legend_str{end + 1} = 'mean';
	plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.mean), 'Color', [0 1 0 ]);
	legend_str{end + 1} = 'robust mean';
	plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.robust_mean), 'Color', [0 0.75 0 ]);
	if  (show_sem)
		legend_str{end + 1} = '+sem';
		legend_str{end + 1} = '-sem';
		plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.mean) - per_size.data(:, per_size.cols.sem), 'Color', [0 0.66 0]);
		plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.mean) + per_size.data(:, per_size.cols.sem), 'Color', [0 0.66 0]);
	end
	if  (show_ci)
		legend_str{end + 1} = '+ci';
		legend_str{end + 1} = '-ci';
		plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.mean) - per_size.data(:, per_size.cols.ci), 'Color', [0 0.37 0]);
		plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.mean) + per_size.data(:, per_size.cols.ci), 'Color', [0 0.37 0]);
	end
	
end

if (show_geomean)
	legend_str{end + 1} = 'geomean';
	plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.geomean), 'Color', [0.5 0.5 0.5 ]);
end
if (show_robust_geomean)
	legend_str{end + 1} = 'robust geomean';
	plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.robust_geomean), 'Color', [0.2 0.2 0.2 ]);
end

if (show_delogged_logmean)
	legend_str{end + 1} = 'delogged logmean';
	plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.delogged_logmean), 'Color', [0 0 0.5]);
end


if(show_median)
	% median +- standard error of the mean, confidence interval would be
	% better
	legend_str{end + 1} = 'median';
	plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.median), 'Color', [1 0 0]);
	if (show_sem)
		legend_str{end + 1} = '+sem';
		legend_str{end + 1} = '-sem';
		plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.median) - per_size.data(:, per_size.cols.sem), 'Color', [0.66 0 0]);
		plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.median) + per_size.data(:, per_size.cols.sem), 'Color', [0.66 0 0]);
	end
	if(show_min)
		% minimum, should be cleanest, but for the test data set looks quite sad...
		legend_str{end + 1} = 'min';
		plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.min), 'Color', [0 0 1]);
	end
	if(show_max)
		% minimum, should be cleanest, but for the test data set looks quite sad...
		legend_str{end + 1} = 'max';
		plot(per_size.data(:, per_size.cols.size), per_size.data(:, per_size.cols.max), 'Color', [0 0 0.66]);
	end
end

title({['If this plot shows a (noisy) step function with a stepping of ', num2str(quantum.byte), ' bytes'], ['then the data carrier is quantised, make sure to use tc-stab']});
xlabel('Approximate packet size [bytes]');
ylabel('ICMP round trip times (ping RTT) [ms]');
legend(legend_str, 'Location', 'NorthWest', 'Interpreter', 'none');
hold off;

if ~isempty(plot_output_format)
	write_out_figure(data_fh, fullfile(sweep_dir, [sweep_name, '_data.', plot_output_format]));
end

% potentially clean up the data, by interpolating values with large sem
% from the neighbours or replacing those with NaNs?

% if the size of the ping packet exceeds the MTU the ping packets gets
% fragmented the step over this ping size will cause a RTT increaser >> one
% RTT_quantum, so exclude all sizes potentially affected by this from the
% search space, (for now assume that the route to the ping host actually can carry 1500 byte MTUs...)
measured_pingsize_idx = find(~isnan(per_size.data(:, per_size.cols.(use_measure))));
tmp_idx = find(measured_pingsize_idx <= max_ping_size_without_fragmentation);
last_non_fragmented_pingsize = measured_pingsize_idx(tmp_idx(end));
ping_sizes_for_linear_fit = measured_pingsize_idx(tmp_idx);

% fit a line to the data, to estimate the RTT per byte
[p, S] = polyfit(per_size.data(ping_sizes_for_linear_fit, per_size.cols.size), per_size.data(ping_sizes_for_linear_fit, per_size.cols.(use_measure)), 1);
RTT_per_byte = p(end - 1);
fitted_line = polyval(p, per_size.data(ping_sizes_for_linear_fit, per_size.cols.size), S);
input_data = per_size.data(ping_sizes_for_linear_fit, per_size.cols.(use_measure));
% estimate the goodness of the linear fit the same way as for the stair
% function
linear_cumulative_difference = sum(abs(input_data - fitted_line));

% figure
% hold on
% plot(per_size.data(ping_sizes_for_linear_fit, per_size.cols.size), per_size.data(ping_sizes_for_linear_fit, per_size.cols.(use_measure)), 'Color', [0 1 0]);
% plot(per_size.data(ping_sizes_for_linear_fit, per_size.cols.size), fitted_line, 'Color', [1 0 0]);
% hold off
% based on the linear fit we can estimate the average RTT per ATM cell
estimated_RTT_quantum_ms = RTT_per_byte * 48;


% just get an idea what range the RTTs per ATM quantum can be for different
% bandwidths
% "ATM" cell over full duplex gigabit ethernet
min_GE_RTT_quantum_ms = (ATM_cell.bit / (1000 * 1000 * 1000) + ATM_cell.bit / (1000 * 1000 * 1000) ) * 1000;	% this estimate is rather a lower bound for fastpath , so search for best fits
% "ATM" cell over theoretical G.fast.vectoring (best case?)
min_GfastV_RTT_quantum_ms = (ATM_cell.bit / (500 * 1000 * 1000) + ATM_cell.bit / (500 * 1000 * 1000) ) * 1000;	% this estimate is rather a lower bound for fastpath , so search for best fits
% the next three are  2014 extreme values fot Deutsche Telekom wired
% assume VDSL2.vectoring 100Mbit 40Mbit
min_VDSL2V_RTT_quantum_ms = (ATM_cell.bit / (100 * 1000 * 1000) + ATM_cell.bit / (40 * 1000 * 1000) ) * 1000;	% this estimate is rather a lower bound for fastpath , so search for best fits
% assume ADSL2+ annex J fallback profile 2J R
max_ADSL2aJ_RTT_quantum_ms = (ATM_cell.bit / (448 * 1000) + ATM_cell.bit / (288 * 1000) ) * 1000;	% this estimate is rather a lower bound for fastpath , so search for best fits
% assume ADSL2+ annex B fixed prifile dsl light 384
max_ADSL1aB_RTT_quantum_ms = (ATM_cell.bit / (384 * 1000) + ATM_cell.bit / (64 * 1000) ) * 1000;	% this estimate is rather a lower bound for fastpath , so search for best fits



% the RTT should equal the average RTT increase per ATM quantum
% estimate the RTT step size
% at ADSL down 3008kbit/sec up 512kbit/sec we expect, this does not include
% processing time
if ~isempty(down_Kbit) || ~isempty(up_Kbit)
	expected_RTT_quantum_ms = (ATM_cell.bit / (down_Kbit * 1000) + ATM_cell.bit / (up_Kbit * 1000) ) * 1000;	% this estimate is rather a lower bound for fastpath , so search for best fits
	%	sm network rates are base 10 nt base 2
	%	expected_RTT_quantum_ms = (ATM_cell.bit / (down_Kbit * 1024) + ATM_cell.bit / (up_Kbit * 1024) ) * 1000;	% this estimate is rather a lower bound for fastpath , so search for best fits
else
	expected_RTT_quantum_ms = estimated_RTT_quantum_ms;
end
disp(['lower bound estimate for one ATM cell RTT based of specified up and downlink is ', num2str(expected_RTT_quantum_ms), ' ms.']);
disp(['estimate for one ATM cell RTT based on linear fit of the ping sweep data is ', num2str(estimated_RTT_quantum_ms), ' ms.']);

% lets search from expected_RTT_quantum_ms to 1.5 * expected_RTT_quantum_ms
% in steps of expected_RTT_quantum_ms / 100
% to allow for interleaved ATM setups increase the search space up to 32
% times best fastpath RTT estimate, 64 interleave seems to add 25ms to the
% per packet latency, but not to the per quantum delta t, so revisit this
% TODO check with high interleave ATM data (if available)
min_search_RTT_ms = expected_RTT_quantum_ms / 2;	% in case the initial estimates are only in the ballpark
search_RTT_steps_ms = expected_RTT_quantum_ms / 100;
max_search_RTT_ms = min([(32 * expected_RTT_quantum_ms) (max_ADSL1aB_RTT_quantum_ms * 1.5)]);
RTT_quantum_list = (min_search_RTT_ms : search_RTT_steps_ms : max_search_RTT_ms);
quantum_list = (1 : 1 : quantum.byte);

% BRUTE FORCE search of best fitting stair...
differences = zeros([length(RTT_quantum_list) length(quantum_list)]);
cumulative_differences = differences;

disp('Starting brute-force search for optimal stair fit, might take a while...');
if (isoctave)
	fflush(stdout); % make octave display intermediate output...
end

all_stairs = zeros([length(RTT_quantum_list) length(quantum_list) length(per_size.data(1:last_non_fragmented_pingsize, per_size.cols.(use_measure)))]);
for i_RTT_quant = 1 : length(RTT_quantum_list)
	cur_RTT_quant = RTT_quantum_list(i_RTT_quant);
	for i_quant = 1 : quantum.byte
		[differences(i_RTT_quant, i_quant), cumulative_differences(i_RTT_quant, i_quant), all_stairs(i_RTT_quant, i_quant, :)] = ...
			get_difference_between_data_and_stair( per_size.data(1:last_non_fragmented_pingsize, per_size.cols.size), per_size.data(1:last_non_fragmented_pingsize, per_size.cols.(use_measure)), ...
			quantum_list(i_quant), quantum.byte, 0, cur_RTT_quant );
	end
end

% for the initial test DSL set the best x_offset was 21, corresponding to 32 bytes overhead before the IP header.
[min_cum_diff, min_cum_diff_idx] = min(cumulative_differences(:));
[min_cum_diff_row_idx, min_cum_diff_col_idx] = ind2sub(size(cumulative_differences),min_cum_diff_idx);
best_difference = differences(min_cum_diff_row_idx, min_cum_diff_col_idx);
disp(['Best staircase fit cumulative difference is: ', num2str(cumulative_differences(min_cum_diff_row_idx, min_cum_diff_col_idx))]);
disp(['Best linear fit cumulative difference is: ', num2str(linear_cumulative_difference)]);
% judge the quantization
if (cumulative_differences(min_cum_diff_row_idx, min_cum_diff_col_idx) < linear_cumulative_difference)
	% stair fits better than line
	quant_string = ['Quantized ATM carrier LIKELY (cummulative residual: stair fit ', num2str(cumulative_differences(min_cum_diff_row_idx, min_cum_diff_col_idx)), ' linear fit ', num2str(linear_cumulative_difference)];
else
	quant_string = ['Quantized ATM carrier UNLIKELY (cummulative residual: stair fit ', num2str(cumulative_differences(min_cum_diff_row_idx, min_cum_diff_col_idx)), ' linear fit ', num2str(linear_cumulative_difference)];
end
disp(quant_string);

disp(['remaining ATM cell length after ICMP header is ', num2str(quantum_list(min_cum_diff_col_idx)), ' bytes.']);
disp(['ICMP RTT of a single ATM cell is ', num2str(RTT_quantum_list(min_cum_diff_row_idx)), ' ms.']);


% as first approximation use the ATM cell offset and known offsets (ICMP
% IPv4 min_ping_size) to estimate the number of cells used for per packet
% overhead
% this assumes that no ATM related overhead is >= ATM cell size
% -1 to account for matlab 1 based indices
% what is the offset in the 2nd ATM cell
n_bytes_overhead_2nd_cell = quantum.byte - (quantum_list(min_cum_diff_col_idx) - 1);	% just assume we can not fit all overhead into one cell...
% what is the known overhead size for the first data point:
tmp_idx = find(~isnan(per_size.data(:, per_size.cols.mean)));
known_overhead_first_ping_size = tmp_idx(1);
%pre_IP_overhead = quantum.byte + (n_bytes_overhead_2nd_cell - known_overhead);	% ths is the one we are after in the end
pre_IP_overhead = quantum.byte + (n_bytes_overhead_2nd_cell - known_overhead_first_ping_size);	% ths is the one we are after in the end
disp(' ');
disp(['Estimated overhead preceding the IP header: ', num2str(pre_IP_overhead), ' bytes']);


res_fh = figure('Name', 'Comparing ping data with');
hold on
legend_str = {'ping data', 'fitted stair', 'fitted line'};
plot(per_size.data(1:last_non_fragmented_pingsize, per_size.cols.size), per_size.data(1:last_non_fragmented_pingsize, per_size.cols.(use_measure)), 'Color', [1 0 0]);
plot(per_size.data(1:last_non_fragmented_pingsize, per_size.cols.size), squeeze(all_stairs(min_cum_diff_row_idx, min_cum_diff_col_idx, :)) + best_difference, 'Color', [0 1 0]);

fitted_line = polyval(p, per_size.data(1:last_non_fragmented_pingsize, per_size.cols.size), S);
plot(per_size.data(1:last_non_fragmented_pingsize, per_size.cols.size), fitted_line, 'Color', [0 0 1]);

title({['Estimated RTT per quantum: ', num2str(RTT_quantum_list(min_cum_diff_row_idx)), ' ms; ICMP data offset in quantum ', num2str(quantum_list(min_cum_diff_col_idx)), ' bytes'];...
	['Estimated overhead preceding the IP header: ', num2str(pre_IP_overhead), ' bytes'];...
	quant_string});
xlabel('Approximate packet size [bytes]');
ylabel('ICMP round trip times (ping RTT) [ms]');
if (isoctave)
	legend(legend_str, 'Location', 'NorthWest');
else
	%annotation('textbox', [0.0 0.95 1.0 .05], 'String', ['Estimated overhead preceding the IP header: ', num2str(pre_IP_overhead), ' bytes'], 'FontSize', 9, 'Interpreter', 'none', 'Color', [1 0 0], 'LineStyle', 'none');
	legend(legend_str, 'Interpreter', 'none', 'Location', 'NorthWest');
end
hold off


%write_out_figure(res_fh, fullfile(sweep_dir, [sweep_name, '_results.pdf'));
if ~isempty(plot_output_format)
	write_out_figure(res_fh, fullfile(sweep_dir, [sweep_name, '_results.', plot_output_format]));
end


% if we have an ATM carrier pre_IP_overhead must be >= 8 byte, otherwise we
% probably are missing an ATM cell full of overhead
if (pre_IP_overhead < 8)
	pre_IP_overhead = pre_IP_overhead + 48;
	disp(['The ATM overhead can not really be smaller than 8 bytes,', sprintf('\n'),...
		'so it seems we have more than one ATM cell worth of overhead', sprintf('\n'),...
		'Adjusted estimated overhead preceding the IP header: ', num2str(pre_IP_overhead)]);
end

% use http://ace-host.stuart.id.au/russell/files/tc/tc-atm/ to present the
% most likely ATM encapsulation for a given overhead and present a recommendation
% for the tc stab invocation
display_protocol_stack_information(pre_IP_overhead);


% now turn this into tc-stab recommendations:
disp(['Add the following to both the egress root qdisc:']);
% disp(' ');
disp(['A) Assuming the router connects over ethernet to the DSL-modem:']);
disp(['stab mtu 2048 tsize 128 overhead ', num2str(pre_IP_overhead), ' linklayer atm']);	% currently tc stab does not account for the ethernet header
% disp(['stab mtu 2048 tsize 128 overhead ', num2str(pre_IP_overhead - offsets.ethernet), ' linklayer atm']);
% disp(' ');
% disp(['B) Assuming the router connects via PPP and non-ethernet to the modem:']);
% disp(['stab mtu 2048 tsize 128 overhead ', num2str(pre_IP_overhead), ' linklayer atm']);

disp(' ');
% on ingress do not exclude the the ethernet header?
disp(['Add the following to both the ingress root qdisc:']);
disp(' ');
disp(['A) Assuming the router connects over ethernet to the DSL-modem:']);
disp(['stab mtu 2048 tsize 128 overhead ', num2str(pre_IP_overhead), ' linklayer atm']);
disp(' ');
if ~(isoctave)
	timestamps.(mfilename).end = toc(timestamps.(mfilename).start);
	disp([mfilename, ' took: ', num2str(timestamps.(mfilename).end), ' seconds.']);
else
	toc
end

% and now the other end of the data, what is the max MTU for the link and
% what is the best ATM cell aligned MTU

disp('Done...');

return
end


function [ ping_data ] = parse_ping_output( ping_log_fqn )
%PARSE_PING_OUTPUT read the putput of a ping run/sweep
% for further processing
% TODO:
%	use a faster parser, using srtok is quite expensive
%
% This currently handles maxosx/linux ping, windows hrping and busybox ping
% windows hrping:
% C:\space\bin\hrping-v506>hrping -n 1 -l 16 www.heise.de
% This is hrPING v5.06.1143 by cFos Software GmbH -- http://www.cfos.de
%
% Source address is 134.2.91.182; using ICMP echo-request, ID=1883
% Pinging www.heise.de [193.99.144.85]
% with 16 bytes data (44 bytes IP):
%
% From 193.99.144.85: bytes=44 seq=0001 TTL=245 ID=b9e6 time=5.031ms
%
% Packets: sent=1, rcvd=1, error=0, lost=0 (0.0% loss) in 0.005031 sec
% RTTs in ms: min/avg/max/dev: 5.031 / 5.031 / 5.031 / 0.000
% Bandwidth in kbytes/sec: sent=8.745, rcvd=8.745
%
%macosx ping:
% hms-beagle:~ moeller$ ping -c 1 -s 16 www.heise.de
% PING www.heise.de (193.99.144.85): 16 data bytes
% 24 bytes from 193.99.144.85: icmp_seq=0 ttl=245 time=4.967 ms
%
% --- www.heise.de ping statistics ---
% 1 packets transmitted, 1 packets received, 0.0% packet loss
% round-trip min/avg/max/stddev = 4.967/4.967/4.967/0.000 ms


if ~(isoctave)
	timestamps.parse_ping_output.start = tic;
else
	tic();
end

verbose = 0;
n_rows_to_grow_table_by = 10000;	% grow table increment to avoid excessive memory copy ops


ping_data = [];
cur_sweep_fd = fopen(ping_log_fqn, 'r');
if (cur_sweep_fd == -1)
	disp(['Could not open ', ping_log_fqn, '.']);
	if isempty(dir(ping_log_fqn))
		disp('Reason: file does not seem to exist at the given directory...')
	end
	return
end
ping_data.header = {'size', 'icmp_seq', 'ttl', 'time'};
ping_data.field_names_list = {'bytes', 'size', 'icmp_seq', 'seq', 'TTL', 'ttl', 'time'};

ping_data.header = {'size', 'time'};	% save half the size...
ping_data.field_names_list = {'bytes', 'size', 'time'};

ping_data.cols = get_column_name_indices(ping_data.header);

ping_data.data = zeros([n_rows_to_grow_table_by, length(ping_data.header)]);
cur_data_lines = 0;
cur_lines = 0;

% skip the first line
% PING netblock-75-79-143-1.dslextreme.com (75.79.143.1): (16 ... 1000)
% data bytes
header_line = fgetl(cur_sweep_fd);

% detect hrping logs, as they data lines look slightly different from unix
% ping
% This is hrPING v5.06.1143 by cFos Software GmbH -- http://www.cfos.de
is_hrping = 0;
if strcmp('This is hrPING ', header_line(1:15))
	is_hrping = 1;
end


while ~feof(cur_sweep_fd)
	% grow the data table if need be
	if (size(ping_data.data, 1) == cur_data_lines)
		if (verbose)
			disp('Growing ping data table...');
		end
		ping_data.data = [ping_data.data; zeros([n_rows_to_grow_table_by, length(ping_data.header)])];
	end
	
	cur_line = fgetl(cur_sweep_fd);
	if ~(mod(cur_lines, 1000))
		disp([num2str(cur_lines +1), ' lines parsed...']);
		if (isoctave)
			fflush(stdout); % make octave display intermediate output...
		end
	end
	cur_lines = cur_lines + 1;
	
	[first_element, remainder] = strtok(cur_line);
	first_element_as_number = str2double(first_element);
	% skip empty & irrelevant lines early
	if isempty(first_element) || strcmp('Request', first_element) || strcmp('---', first_element) ...
			|| strcmp('Source', first_element) || strcmp('Pinging', first_element) || strcmp('with', first_element) || strcmp('Packets:', first_element) ...
			|| strcmp('RTTs', first_element) || strcmp('Bandwidth', first_element)
		% skip empty lines explicitly
		continue;
	end
	% the following will not work for merged ping
	%if strmatch('---', first_element)
	%	%we reached the end of sweeps
	%	break;
	%end
	% now read in the data
	%unix ping:		30 bytes from 75.79.143.1: icmp_seq=339 ttl=63 time=14.771 ms
	%hrping:		From 193.99.144.85: bytes=44 seq=0001 TTL=245 ID=b9e6 time=5.031ms
	if (~isempty(first_element_as_number) && ~isnan(first_element_as_number)) || (strcmp('From', first_element) && (is_hrping))
		% get the next element
		[tmp_next_item, tmp_remainder] = strtok(remainder);
		if strcmp(tmp_next_item, 'bytes') || is_hrping
			if ~(mod(cur_data_lines, 1000))
				disp(['Milestone ', num2str(cur_data_lines +1), ' ping packets reached...']);
				if (isoctave)
					fflush(stdout); % make octave display intermediate output...
				end
			end
			cur_data_lines = cur_data_lines + 1;
			
			% size of the ICMP package
			ping_data.data(cur_data_lines, ping_data.cols.size) = first_element_as_number;	% attention for hrping this is a NaN...
			% now process the remainder
			while ~isempty(remainder)
				[next_item, remainder] = strtok(remainder);
				equality_pos = strfind(next_item, '=');
				% data items are name+value pairs
				if ~isempty(equality_pos);
					cur_key = next_item(1: equality_pos - 1);
					cur_value = str2double(next_item(equality_pos + 1: end));
					%hr_ping reports time as time=5.031ms insted of unix's
					%time=14.771 ms, so handle the ms by hand
					if (is_hrping) && strcmp('ms', next_item(end-1:end))
						cur_value = str2double(next_item(equality_pos + 1: end-2));
					end
					if (ismember(cur_key, ping_data.field_names_list))
						switch cur_key
							% busybox ping and macosx ping return different key names
							case {'seq', 'icmp_seq'}
								ping_data.data(cur_data_lines, ping_data.cols.icmp_seq) = cur_value;
							case {'ttl', 'TTL'}
								ping_data.data(cur_data_lines, ping_data.cols.ttl) = cur_value;
							case 'time'
								ping_data.data(cur_data_lines, ping_data.cols.time) = cur_value;
							case 'bytes'
								ping_data.data(cur_data_lines, ping_data.cols.size) = cur_value;	%hrping reports the size as bytes=44
						end
					end
				end
			end
		else
			% skip this line
			if (verbose)
				disp(['Skipping: ', cur_line]);
			end
		end
	else
		if (verbose)
			disp(['Ping output: ', cur_line, ' not handled yet...']);
		end
	end
	
end

% remove empty lines
if (size(ping_data.data, 1) > cur_data_lines)
	ping_data.data = ping_data.data(1:cur_data_lines, :);
end

disp(['Found ', num2str(cur_data_lines), ' ping packets in ', ping_log_fqn]);
% clean up
fclose(cur_sweep_fd);

if ~(isoctave)
	timestamps.parse_ping_output.end = toc(timestamps.parse_ping_output.start);
	disp(['Parsing took: ', num2str(timestamps.parse_ping_output.end), ' seconds.']);
else
	toc
end

return
end


function [ difference , cumulative_difference, stair_y ] = get_difference_between_data_and_stair( data_x, data_y, x_size, stair_x_step_size, y_offset, stair_y_step_size )
% 130619sm: handle NaNs in data_y (marker for missing ping sizes)
% x_size is the flat part of the first stair, that is quantum minus the
% offset
% TODO: understand the offset issue and simplify this function
%		extrapolate the stair towards x = 0 again

debug = 0;
difference = [];

tmp_idx = find(~isnan(data_y));
x_start_val_idx = tmp_idx(1);
x_start_val = data_x(x_start_val_idx);
x_end_val = data_x(end);	% data_x is sorted...

% construct stair
stair_x = data_x;
proto_stair_y = zeros([x_end_val 1]);	% we need the final value in
% make sure the x_size values do not exceed the step size...
if (x_size > stair_x_step_size)
	if mod(x_size, stair_x_step_size) == 0
		x_size = stair_x_step_size;
	else
		x_size = mod(x_size, stair_x_step_size);
	end
end

%stair_y_step_idx = (x_start_val + x_size : stair_x_step_size : x_end_val);
%% we really want steps registered to x_start_val
%stair_y_step_idx = (mod(x_start_val, stair_x_step_size) + x_size : stair_x_step_size : x_end_val);
stair_y_step_idx = (mod(x_start_val + x_size, stair_x_step_size) : stair_x_step_size : x_end_val);
if stair_y_step_idx(1) == 0
	stair_y_step_idx(1) = [];
end

proto_stair_y(stair_y_step_idx) = stair_y_step_size;
stair_y = cumsum(proto_stair_y);
if (debug)
	figure
	hold on;
	title(['x offset used: ', num2str(x_size), ' with quantum ', num2str(stair_x_step_size)]);
	plot(data_x, data_y, 'Color', [0 1 0]);
	plot(stair_x, stair_y, 'Color', [1 0 0]);
	hold off;
end
% missing ping sizes are filled with NaNs, so skip those
notnan_idx = find(~isnan(data_y));
% estimate the best y_offset for the stair
difference = sum(abs(data_y(notnan_idx) - stair_y(notnan_idx))) / length(data_y(notnan_idx));
% calculate the cumulative difference between stair and data...
cumulative_difference = sum(abs(data_y(notnan_idx) - (stair_y(notnan_idx) + difference)));

return
end

% function [ stair ] = build_stair(x_vector, x_size, stair_x_step_size, y_offset, stair_y_step_size )
% stair = [];
%
% return
% end


function [columnnames_struct, n_fields] = get_column_name_indices(name_list)
% return a structure with each field for each member if the name_list cell
% array, giving the position in the name_list, then the columnnames_struct
% can serve as to address the columns, so the functions assitgning values
% to the columns do not have to care too much about the positions, and it
% becomes easy to add fields.
n_fields = length(name_list);
for i_col = 1 : length(name_list)
	cur_name = name_list{i_col};
	columnnames_struct.(cur_name) = i_col;
end
return
end




function [ci_halfwidth_vector] = calc_cihw(std_vector, n, alpha)
%calc_ci : calculate the half width of the confidence interval (for 1 - alpha)
%	the t_value lookup depends on alpha and the samplesize n; the relevant
%	calculation of the degree of freedom is performed inside calc_t_val.
%	ci_halfwidth = t_val(alpha, n-1) * std / sqrt(n)
%	Each groups CI ranges from mean - ci_halfwidth to mean - ci_halfwidth, so
%	the calling function has to perform this calculation...
%
% INPUTS:
%	std_vector: vector containing the standard deviations of all requested
%		groups
%	n: number of samples in each group, if the groups have different
%		samplesizes, specify each group's sample size in a vector
%	alpha: the desired maximal uncertainty/error in the range of [0, 1]
% OUTPUT:
%	ci_halfwidth_vector: vector containing the confidence intervals half width
%		for each group

% calc_t_val return one sided t-values, for the desired two sidedness one has
% to half the alpha for the table lookup
cur_alpha = alpha / 2;

% if n is scalar use same n for all elements of std_vec
if isscalar(n)
	t_ci = calc_t_val(cur_alpha, n);
	ci_halfwidth_vector = std_vector * t_ci / sqrt(n);
	% if n is a vector, prepare a matching vector of t_ci values
elseif isvector(n)
	t_ci_vector = n;
	% this is probably ugly, but calc_t_val only accepts scalars.
	for i_pos = 1 : length(n)
		t_ci_vector(i_pos) = calc_t_val(cur_alpha, n(i_pos));
	end
	ci_halfwidth_vector = std_vector .* t_ci_vector ./ sqrt(n);
end

return
end

%-----------------------------------------------------------------------------
function [t_val] = calc_t_val(alpha, n)
% the t value for the given alpha and n
% so call with the n of the sample, not with degres of freedom
% see http://mathworld.wolfram.com/Studentst-Distribution.html for formulas
% return values follow Bortz, Statistik fuer Sozialwissenschaftler, Springer
% 1999, table D page 775. That is it returns one sided t-values.
% primary author S. Moeller

% TODO:
%	sidedness of t-value???

% basic error checking
if nargin < 2
	error('alpha and n have to be specified...');
end

% probabilty of error
tmp_alpha = alpha ;%/ 2;
if (tmp_alpha < 0) || (tmp_alpha > 1)
	msgbox('alpha has to be taken from [0, 1]...');
	t_val = NaN;
	return
end
if tmp_alpha == 0
	t_val = -Inf;
	return
elseif tmp_alpha ==1
	t_val = Inf;
	return
end
% degree of freedom
df = n - 1;
if df < 1
	%msgbox('The n has to be >= 2 (=> df >= 1)...');
	% 	disp('The n has to be >= 2 (=> df >= 1)...');
	t_val = NaN;
	return
end


% only calculate each (alpha, df) combination once, store the results
persistent t_val_array;
% create the t_val_array
if ~iscell(t_val_array)
	t_val_array = {[NaN;NaN]};
end
% search for the (alpha, df) tupel, avoid calculation if already stored
if iscell(t_val_array)
	% cell array of 2d arrays containing alpha / t_val pairs
	if df <= length(t_val_array)
		% test whether the required alpha, t_val tupel exists
		if ~isempty(t_val_array{df})
			% search for alpha
			tmp_array = t_val_array{df};
			alpha_index = find(tmp_array(1,:) == tmp_alpha);
			if any(alpha_index)
				t_val = tmp_array(2, alpha_index);
				return
			end
		end
	else
		% grow t_val_array to length of n
		missing_cols = df - length(t_val_array);
		for i_missing_cols = 1: missing_cols
			t_val_array{end + 1} = [NaN;NaN];
		end
	end
end

% check the sign
cdf_sign = 1;
if (1 - tmp_alpha) == 0.5
	t_val = t_cdf;
elseif (1 - tmp_alpha) < 0.5 % the t-cdf is point symmetric around (0, 0.5)
	cdf_sign = -1;
	tmp_alpha = 1 - tmp_alpha; % this will be undone later
end

% init some variables
n_iterations = 0;
delta_t = 1;
last_alpha = 1;
higher_t = 50;
lower_t = 0;
% find a t-value pair around the desired alpha value
while norm_students_cdf(higher_t, df) < (1 - tmp_alpha);
	lower_t = higher_t;
	higher_t = higher_t * 2;
end

% search the t value for the given alpha...
while (n_iterations < 1000) && (abs(delta_t) >= 0.0001)
	n_iterations = n_iterations + 1;
	% get the test_t (TODO linear interpolation)
	% higher_alpha = norm_students_cdf(higher_t, df);
	% lower_alpha = norm_students_cdf(lower_t, df);
	test_t = lower_t + ((higher_t - lower_t) / 2);
	cur_alpha = norm_students_cdf(test_t, df);
	% just in case we hit the right t spot on...
	if cur_alpha == (1 - tmp_alpha)
		t_crit = test_t;
		break;
		% probably we have to search for the right t
	elseif cur_alpha < (1 - tmp_alpha)
		% test_t is the new lower_t
		lower_t = test_t;
		%higher_t = higher_t;	% this stays as is...
	elseif cur_alpha > (1 - tmp_alpha)
		%
		%lower_t = lower_t;	% this stays as is...
		higher_t = test_t;
	end
	delta_t = higher_t - lower_t;
	last_alpha = cur_alpha;
end
t_crit = test_t;

% set the return value, correct for negative t values
t_val = t_crit * cdf_sign;
if cdf_sign < 0
	tmp_alpha = 1 - tmp_alpha;
end

% store the alpha, n, t_val tupel in t_val_array
pos = size(t_val_array{df}, 2);
t_val_array{df}(1, (pos + 1)) = tmp_alpha;
t_val_array{df}(2, (pos + 1)) = t_val;

return
end

%-----------------------------------------------------------------------------
function [scaled_cdf] = norm_students_cdf(t, df)
% calculate the cdf of students distribution for a given degree of freedom df,
% and all given values of t, then normalize the result
% the extreme values depend on the values of df!!!

% get min and max by calculating values for extrem t-values (e.g. -10000000,
% 10000000)
extreme_cdf_vals = students_cdf([-10000000, 10000000], df);

tmp_cdf = students_cdf(t, df);

scaled_cdf =	(tmp_cdf - extreme_cdf_vals(1)) /...
	(extreme_cdf_vals(2) - extreme_cdf_vals(1));
return
end

%-----------------------------------------------------------------------------
function [cdf_value_array] = students_cdf(t_value_array, df)
%students_cdf: calc the cumulative density function for a t-distribution
% Calculate the CDF value for each value t of the input array
% see http://mathworld.wolfram.com/Studentst-Distribution.html for formulas
% INPUTS:	t_value_array:	array containing the t values for which to
%							calculate the cdf
%			df:	degree of freedom; equals n - 1 for the t-distribution

cdf_value_array = 0.5 +...
	((betainc(1, 0.5 * df, 0.5) / beta(0.5 * df, 0.5)) - ...
	(betainc((df ./ (df + t_value_array.^2)), 0.5 * df, 0.5) /...
	beta(0.5 * df, 0.5))) .*...
	sign(t_value_array);

return
end

%-----------------------------------------------------------------------------
function [t_prob_dist] = students_pf(df, t_arr)
%  calculate the probability function for students t-distribution

t_prob_dist =	(df ./ (df + t_arr.^2)).^((1 + df) / 2) /...
	(sqrt(df) * beta(0.5 * df, 0.5));

% % calculate and scale the cdf by hand...
% cdf = cumsum(t_prob_dist);
% discrete_t_cdf = (cdf - min(cdf)) / (max(cdf) - min(cdf));
% % numericaly get the t-value for the given alpha
% tmp_index = find(discrete_t_cdf > (1 - tmp_alpha));
% t_crit = t(tmp_index(1));

return
end

function in = isoctave ()
persistent inout;

if isempty(inout),
	inout = exist('OCTAVE_VERSION','builtin') ~= 0;
end;
in = inout;

return;
end


function [] = display_protocol_stack_information(pre_IP_overhead)
% use [1] http://ace-host.stuart.id.au/russell/files/tc/tc-atm/ to present the
% most likely ATM protocol stack setup for a given overhead so the user can
% compare with his prior knowledge

% how much data fits into ATM cells without padding? 32 cells would be 1519
% which is larger than the 1500 max MTU for ethernet
ATM_31_cells_proto_MTU = 31 * 48;	% according to [1] 31 cells are the optimum for all protocol stacks
ATM_32_cells_proto_MTU = 32 * 48;	% should be best for case 44

disp(' ');
disp('According to http://ace-host.stuart.id.au/russell/files/tc/tc-atm/');
disp(['', num2str(pre_IP_overhead), ' bytes overhead indicate']);

switch pre_IP_overhead
	case 8
		disp('Connection: IPoA, VC/Mux RFC-2684');
		disp('Protocol (bytes): ATM AAL5 SAR (8) : Total 8');
		overhead_bytes_around_MTU = 8;
		overhead_bytes_in_MTU = 0;
		
	case 16
		disp('Connection: IPoA, LLC/SNAP RFC-2684');
		disp('Protocol (bytes): ATM LLC (3), ATM SNAP (5), ATM AAL5 SAR (8) : Total 16');
		overhead_bytes_around_MTU = 16;
		overhead_bytes_in_MTU = 0;
		
	case 24
		disp('Connection: Bridged, VC/Mux RFC-1483/2684');
		disp('Protocol (bytes): Ethernet Header (14), ATM pad (2), ATM AAL5 SAR (8) : Total 24');
		overhead_bytes_around_MTU = 24;
		overhead_bytes_in_MTU = 0;
		
	case 28
		disp('Connection: Bridged, VC/Mux+FCS RFC-1483/2684');
		disp('Protocol (bytes): Ethernet Header (14), Ethernet PAD [8] (0), Ethernet Checksum (4), ATM pad (2), ATM AAL5 SAR (8) : Total 28');
		overhead_bytes_around_MTU = 28;
		overhead_bytes_in_MTU = 0;
		
	case 32
		disp('Connection: Bridged, LLC/SNAP RFC-1483/2684');
		disp('Protocol (bytes): Ethernet Header (14), ATM LLC (3), ATM SNAP (5), ATM pad (2), ATM AAL5 SAR (8) : Total 32');
		overhead_bytes_around_MTU = 32;
		overhead_bytes_in_MTU = 0;
		disp('OR');
		disp('Connection: PPPoE, VC/Mux RFC-2684');
		disp('Protocol (bytes): PPP (2), PPPoE (6), Ethernet Header (14), ATM pad (2), ATM AAL5 SAR (8) : Total 32');
		overhead_bytes_around_MTU = 24;
		overhead_bytes_in_MTU = 8;
		
	case 36
		disp('Connection: Bridged, LLC/SNAP+FCS RFC-1483/2684');
		disp('Protocol (bytes): Ethernet Header (14), Ethernet PAD [8] (0), Ethernet Checksum (4), ATM LLC (3), ATM SNAP (5), ATM pad (2), ATM AAL5 SAR (8) : Total 36');
		overhead_bytes_around_MTU = 36;
		overhead_bytes_in_MTU = 0;
		disp('OR');
		disp('Connection: PPPoE, VC/Mux+FCS RFC-2684');
		disp('Protocol (bytes): PPP (2), PPPoE (6), Ethernet Header (14), Ethernet PAD [8] (0), Ethernet Checksum (4), ATM pad (2), ATM AAL5 SAR (8) : Total 36');
		overhead_bytes_around_MTU = 28;
		overhead_bytes_in_MTU = 8;
		
	case 10
		disp('Connection: PPPoA, VC/Mux RFC-2364');
		disp('Protocol (bytes): PPP (2), ATM AAL5 SAR (8) : Total 10');
		overhead_bytes_around_MTU = 8;
		overhead_bytes_in_MTU = 2;
		
	case 14
		disp('Connection: PPPoA, LLC RFC-2364');
		disp('Protocol (bytes): PPP (2), ATM LLC (3), ATM LLC-NLPID (1), ATM AAL5 SAR (8) : Total 14');
		overhead_bytes_around_MTU = 12;
		overhead_bytes_in_MTU = 2;
		
	case 40
		disp('Connection: PPPoE, LLC/SNAP RFC-2684');
		disp('Protocol (bytes): PPP (2), PPPoE (6), Ethernet Header (14), ATM LLC (3), ATM SNAP (5), ATM pad (2), ATM AAL5 SAR (8) : Total 40');
		overhead_bytes_around_MTU = 32;
		overhead_bytes_in_MTU = 8;
		
	case 44
		disp('Connection: PPPoE, LLC/SNAP+FCS RFC-2684');
		disp('Protocol (bytes): PPP (2), PPPoE (6), Ethernet Header (14), Ethernet PAD [8] (0), Ethernet Checksum (4), ATM LLC (3), ATM SNAP (5), ATM pad (2), ATM AAL5 SAR (8) : Total 44');
		overhead_bytes_around_MTU = 36;
		overhead_bytes_in_MTU = 8;
		disp('OR');
		disp('Connection: PPPoE, LLC/SNAP RFC-2684 + VLAN tag terminated at modem');
		disp('Protocol (bytes): VLAN tag (4), PPP (2), PPPoE (6), Ethernet Header (14), ATM LLC (3), ATM SNAP (5), ATM pad (2), ATM AAL5 SAR (8) : Total 44');
		overhead_bytes_around_MTU = 36;
		overhead_bytes_in_MTU = 8;
		
	case {0, 48}
		disp('Overhead of 0 bytes is not possible so assume 1 full packet (48 bytes) overhead...');
		disp('Connection: Bridged, LLC/SNAP+FCS RFC-1483/2684 + VLAN tag terminated at modem');
		disp('Protocol (bytes): Ethernet Header (14), VLAN tag (4), Ethernet PAD [8] (0), Ethernet Checksum (4), ATM LLC (3), ATM SNAP (5), ATM pad (2), ATM AAL5 SAR (8) : Total 36');
		overhead_bytes_around_MTU = 36;
		overhead_bytes_in_MTU = 8;
		
		
	otherwise
		disp('a protocol stack this program does NOT know (yet)...');
end

disp(' ');
return;
end


function range_mean = robust_mean(value_list, lower_limit_ratio, upper_limit_ratio)

n_vals = length(value_list);
sorted_values = sort(value_list);

lowest_robust_idx = max([ceil(n_vals * lower_limit_ratio), 1]);
highest_robust_idx = max([floor(n_vals * upper_limit_ratio), 1]);

range_mean = mean(sorted_values(lowest_robust_idx:highest_robust_idx));

return
end

function range_geomean = robust_geomean(value_list, lower_limit_ratio, upper_limit_ratio)

n_vals = length(value_list);
sorted_values = sort(value_list);

lowest_robust_idx = max([ceil(n_vals * lower_limit_ratio), 1]);	% allow lower_limit_ratio = 0
highest_robust_idx = max([floor(n_vals * upper_limit_ratio), 1]);

range_geomean = geomean(sorted_values(lowest_robust_idx:highest_robust_idx));

return
end


function [ ret_val ] = write_out_figure(img_fh, outfile_fqn)
%WRITE_OUT_FIGURE save the figure referenced by img_fh to outfile_fqn,
% using .ext of outfile_fqn to decide which image type to save as.
%   Detailed explanation goes here
% write out the data

% check whether the path exists, create if not...
[pathstr, name, img_type] = fileparts(outfile_fqn);
if isempty(dir(pathstr)),
	mkdir(pathstr);
end

switch img_type(2:end)
	case 'pdf'
		% pdf in 7.3.0 is slightly buggy...
		print(img_fh, '-dpdf', outfile_fqn);
	case 'ps'
		print(img_fh, '-depsc2', outfile_fqn);
	case 'tiff'
		% tiff creates a figure
		print(img_fh, '-dtiff', outfile_fqn);
	case 'png'
		% tiff creates a figure
		print(img_fh, '-dpng', outfile_fqn);
		%       case 'tif'
		%               % tif only creates the image
		%               % write out the mosaic as image, sadly the compression does not work...
		%               imwrite(mos, outfile_fqn, img_type, 'Compression', 'none');
	case 'fig'
		%sm: allows to save figures for further refinements
		saveas(img_fh, outfile_fqn, 'fig');
	otherwise
		% default to uncompressed images
		disp(['Image type: ', img_type, ' not handled yet, passing on to print unchanged (might not work).']);
		print(img_fh, outfile_fqn);
		% write out the mosaic as image, sadly the compression does not work...
		%               imwrite(mos, [outfile_fqn, '.tif'], format, 'Compression', 'none');
end

if ~isnumeric(img_fh)
	disp(['Saved current figure to: ', outfile_fqn]);
else
	disp(['Saved figure (', num2str(img_fh), ') to: ', outfile_fqn]);
end

ret_val = 0;
return
end


function [pkg_is_loadable, pkg_is_loaded] = check_octave_pkg_availability(pkg_name)
% figure out whether an octave pkg can be loaded...
pkg_is_loadable = 0;
pkg_is_loaded = 0;

tmp_pkg_list = pkg('list');
n_avaliable_pkg = size(tmp_pkg_list, 2);

% search through all avalable pkgs
for i_pkg = 1 : n_avaliable_pkg
	cur_pkg_struct = tmp_pkg_list{1, i_pkg};
	if strcmp(pkg_name, cur_pkg_struct.name)
		pkg_is_loadable = 1;
		pkg_is_loaded = cur_pkg_struct.loaded;
	end
end

if ~(pkg_is_loadable)
	disp(['The requested octave pkg: ', pkg_name, ' is not available on this host']);
	tmp_forge_pkg_list = pkg('list', '-forge');
	n_avaliable_forge_pkg = size(tmp_forge_pkg_list, 2);
	% search through all avalable pkgs
	for i_forge_pkg = 1 : n_avaliable_forge_pkg
		cur_forge_pkg = tmp_forge_pkg_list{1, i_forge_pkg};
		if strcmp(pkg_name, cur_forge_pkg)
			disp(['The requested pkg: ',pkg_name , ' seems available on octave-forge; consider trying:']);
			disp(['pkg(''install'', ''-forge'', ''', pkg_name, ''', ''-verbose'')']);
		end
	end
end
fflush(stdout); % make octave display intermediate output...
return
end
