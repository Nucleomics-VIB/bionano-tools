#!/bin/sh

# A script to log the cpu watt-usage and temp of xeon-phi cards
# launch on Thinkmate when starting assembly job
# scriptname: logphicards.sh
#
# inspired by http://www.unix.com/shell-programming-and-scripting/
#   223177-shell-script-logging-cpu-memory-usage-linux-process.html
#
# Stephane Plaisance VIB-BITS march-14-2016 v1.0
# fixed small bugs (v1.1; 2016-03-15)
#
# visit our Git: https://github.com/Nucleomics-VIB

usage='# Usage: logphicards.sh
#    -t <log-frequency in sec (default 60sec)>'

while getopts "t:h" opt; do
  case $opt in
    t) timeint=${OPTARG} ;;
    h) echo "${usage}" >&2; exit 0 ;;
  esac
done

maxmic=5

# repeat loop every x sec
FREQ=${timeint:-60}

# test if time is a valid number
re='^[0-9]+$'
if ! [[ ${FREQ} =~ $re ]] ; then
echo "not an integer number of seconds"
echo "${usage}"
exit 1
fi

# log file name and header
LOG_FILE=Xeon_usage_"$(date +%s)".log
echo -e "logtime\tmic\tcpu_user\tcpuT\tmemT\ttotW" > $LOG_FILE

echo "# logging to $LOG_FILE every $FREQ sec"
echo "# press <Ctrl>-C to stop logging"

# infinite loop will run until ctrl-C is hit
while :; do

for mic in $(eval echo "mic{0..${maxmic}}"); do
t=$(date +%s)
(echo -ne "${t}\t${mic}\t"; micsmc -c ${mic} -t ${mic} -f ${mic} | \
egrep "Device Utilization:|Cpu Temp:|Memory Temp:|Total Power:" | \
awk '{if (FNR==1) {match($0, /([0-9\.]+)/, arr); if(arr[1] != "") print arr[1]} 
	else {print substr($0, 31, length($0)-29)}}' | \
cut -d "%" -f 1 | cut -d " " -f 1 | transpose -t) >> ${LOG_FILE}
done
sleep ${FREQ}

done

# when the app closes, no additional line will be added to the log file