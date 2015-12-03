#! /bin/bash
# Note this does not work with dash, please either use "bash ./ping_collector.sh" or better yet
# "chmod u+x ./ping_collector.sh ; ./ping_collector.sh" to make sure this is interpreted by bash
#
# TODO use seq or bash to generate a list of the requested sizes (to allow for non-equdistantly spaced sizes)

# just an identifier for the ping log
TECH=""

# finding a proper target IP is somewhat of an art, just traceroute a remote site 
# and find the nearest host reliably responding to pings showing the smallet variation of pingtimes
# for this I typically run "traceroute 8.8.8.8", and then select the first host on the ISP side (typically after 
# the first large RTT increment) and test its response by "ping -c 10 -s 16 NNN.NNN.NNN.NNN", if this host does not repsond 
# I pick the next host along the route to 8.8.8.8. I assume the closer the host the less disturbed by other traffic the 
# response will be.


if [ ! $# == 1 ]; then
    echo "To run measurements supply the TARGET IP address as first agument to ${0} this script."
    echo "Use traceroute 8.8.8.8 to get a list of increasingly distant hosts, pick the first host out of your network (ideally the DSLAM)."
    echo "Test whether the selected host responds to ping: 'ping -s16 -c 1 target.IP.address.quad' : this needs to actually return non zero RTTs."
    echo "If the hosts does not reply to the pings take the next host from the traceroute (movin closer to 8.8.8.8), repeat until you find a replying host."
    echo "Once the main script is started have a quick look at the logfile, to see whether the RTTs stay close to the initial test RTT."
    echo "If the RTTs have increased a lot, the PINGPERIOD might be too short, and the host might have put us on a slow path; either increase PINGPERIOD or try the next host..."
    echo ""
    echo "Here is the traceroute (might take a while):"
    echo ""
    traceroute 8.8.8.8
    echo ""
    echo "Alternatively you might try to use googles infrastructure by running: ${0} gstatic.com "
    
    exit 0
else
    TARGET=${1}		# Replace by an appropriate host
fi


DATESTR=`date +%Y%m%d_%H%M%S`	# to allow multiple sequential records
LOG=ping_sweep_${TECH}_${DATESTR}.txt


# by default non-root ping will only end one packet per second, so work around that by calling ping independently for each package
# empirically figure out the shortest period still giving the standard ping time (to avoid being slow-pathed by our host)
# at 100 packets/s of 116 + 28 + 40 we would need 4 ATM cells = 192byte * 100/s = 150kbit/s
# at 100 packets/s of 16 + 28 + 40nwe would need 2 ATM cells = 96byte * 100/s = 75kbit/s
# on average we need 150 + 75 * 0.5 = 112.5 Kbit/s, increase the ping period if uplink < 112.5 Kbit/s
PINGPERIOD=0.01		# increase if uplink slower than roughly 200Kbit/s
PINGSPERSIZE=10000	# the higher the link rate the more samples we need to reliably detect the increasingly smaller ATM quantisation steps. Can be reduced for slower links

# Start, needed to find the per packet overhead dependent on the ATM encapsulation
# to reliably show ATM quantization one would like to see at least two steps, so cover a range > 2 ATM cells (so > 96 bytes)
# Note to be more robust use 3
#SWEEPMINSIZE=16		# 64bit systems seem to require 16 bytes of payload to include a timestamp... so use 16 as minimum
#SWEEPMAXSIZE=116
SWEEPMAXSIZE=166	# this contains 3 full cells so more transitions to pin the quantization offset to...
    

n_SWEEPS=`expr ${SWEEPMAXSIZE} - ${SWEEPMINSIZE}`


i_sweep=0
i_size=0

while [ ${i_sweep} -lt ${PINGSPERSIZE} ]
do
    (( i_sweep++ ))
    echo "Current iteration: ${i_sweep}"
    # now loop from sweepmin to sweepmax
    i_size=${SWEEPMINSIZE}
    while [ ${i_size} -le ${SWEEPMAXSIZE} ]
    do
	echo "${i_sweep}. repetition of ping size ${i_size}"
	ping -c 1 -s ${i_size} ${TARGET} >> ${LOG} &
	(( i_size++ ))
	# we need a sleep binary that allows non integer times (GNU sleep is fine as is sleep of macosx 10.8.4)
	sleep ${PINGPERIOD}
    done
done

#tail -f ${LOG}

echo "Done... ($0)
"