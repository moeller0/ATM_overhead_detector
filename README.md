# ATM_overhead_detector
Heuristically determine likelihood of ATM-AAL5 "quantization" and likely per packet overhead for ATM based links.

Basically, this is just following the ideas established by Jesper Dangaard Brouer in https://web.archive.org/web/20130811042603/http://adsl-optimizer.dk/thesis/main_final_hyper.pdf and by Russell Stuart (https://web.archive.org/web/20150606220856/http://ace-host.stuart.id.au/russell/files/tc/tc-atm/) so none if this is actually new or original.

In a nutshell, this project consists out of two stages, a data collection phase and an atm detection phase. 
	The first phase is performed by a small shell script, that, on unix machines, should collect a set ICMP (ping) probes of systematically differing sizes to a remote site.
	The second phase then parses the results and tries to detect whether one of the links on the network path was/is affected by ATM AAL5 quantization; in addition it also tries to estimate the amount of per packet overhead that is applied on the ATM-AAL5 link. Note that the 2nd phase will always generate an overhead estimate even on non-ATM links, where the overhead is going to be probably wrong. The user is advised to apply good judgment in seeing how well the estimated quantised "stair"-function fits the empirical data.
  
The wiki contains a bit more information: https://github.com/moeller0/ATM_overhead_detector/wiki  
  
  
  

Instructions:

0) Read ping_collector.sh (unix) or ping_collector.bat (windows)

1) optionally edit parameters in ping_collector.sh/ping_collector.bat

2) Exceute the following from a terminal "./ping_collector.sh suitable.remote.host.IP_or_address" for windows run "ping_collector.bat" inside a command window (CMD.EXE, should be started with the "Run as Administrator" option), make sure to first install hrping (https://www.cfos.de/en/ping/ping.htm)

3) wait until the script finishes (might take a few hours, basically PINGSPERSIZE * PINGPERIOD * n_SWEEPS in seconds)

4) run ATM_overhead_detector.m in either matlab or octave* and look at the output (load the resulting output file from the ping_collector run)

NOTE: The estimated overhead will only be correct for ATM links, so if you are certain your link does not use ATM/AAL5, don't bother to use use ATM_overhead_detector, it will not give you reliable information about your link's per-paket-overhead...



*) If the octave statistics package is selected the geometric mean becomes available as measure to use, if the pkg is not available the script simply disables geomean and selects median as default use_measure.

Tests with the fltk backend on both linux and macosx caused octave/ghostscript crashes, so this now defaults to gnuplot under octave

Especially the parser for the output file of the 1st phase is slow, especially under octave; case in point parsing a ~400 MB file contaning 1506625 ping packets took 6420.28 seconds (1 hour 47 minutes and 0.28 seconds) under octave but only 2524.7899 seconds (42 minutes and 4.7899 seconds) under matlab. Fortunetely each file only needs to be parsed once...


Now, it would be great if users of this code could post the "Quantized ATM carrier LIKELY (cummulative residual: stair fit 1.9738 linear fit 14.5684" information from the run as a new issue, together with the information whether the link truely is using an ATM carrier or not. The goal is to figure out a better heuristic to declare ATM carrier likelyhood than simply comparing the cummulative residuals...
