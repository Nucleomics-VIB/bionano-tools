#!/bin/bash

# Compute stats for all cmap files after HS run using the BNG calc_cmap.pl
# return all results in a transposed tab-separated table
#
# script name: HS_stats.sh
# Requirements:
# run on a unix computer installed with working bionano code
# unix executables transpose and columns installed
# $BNG_HSSCRIPTS should point to the location of the BNG perl script

# Stephane Plaisance (VIB-NC+BITS) 2016/12/05; v1.0
#
# visit our Git: https://github.com/BITS-VIB

# run within the assembly folder (root of all HS runs)
infolder=$1
refcmap="$(ls *.cmap)"
infasta="$(basename $(ls ${infolder}/*.fa* | head -1))"
incmap="$(basename $(ls ${infolder}/output/fa2cmap/*.cmap))"
results="${infolder}/HS_stats_${infasta}.txt"

# edit if required
BNG_HSSCRIPTS=$BNG_SCRIPTS/HybridScaffold/scripts

echo "# computing HS stats for ${infasta} vs ${refmap}"

(echo -n "${infasta}"$'\t'; perl $BNG_HSSCRIPTS/calc_cmap_stats.pl ${refcmap} | tr "=" "\t" | transpose -t | head -1 | sed -e 's/ \+/_/g' | sed -e 's/_\b//g' | column -t) > ${results}
(echo -n "${refcmap}"$'\t'; perl $BNG_HSSCRIPTS/calc_cmap_stats.pl ${refcmap} | tr "=" "\t" | transpose -t | tail -1 | column -t) >>  ${results}
(echo -n "${incmap}"$'\t'; perl $BNG_HSSCRIPTS/calc_cmap_stats.pl ${infolder}/output/fa2cmap/${incmap} | tr "=" "\t" | transpose -t | tail -1 | column -t) >> ${results}
(echo -n "filtered NGS"$'\t'; perl $BNG_HSSCRIPTS/calc_cmap_stats.pl ${infolder}/output/align_final/filtered_NGS.cmap | tr "=" "\t" | transpose -t | tail -1 | column -t) >> ${results}
(echo -n "filtered-not-used NGS"$'\t'; perl $BNG_HSSCRIPTS/calc_cmap_stats.pl ${infolder}/output/align_final/filtered_not_used_NGS.cmap | tr "=" "\t" | transpose -t | tail -1 | column -t) >> ${results}
(echo -n "filtered BNG"$'\t'; perl $BNG_HSSCRIPTS/calc_cmap_stats.pl ${infolder}/output/align_final/filtered_BNG.cmap | tr "=" "\t" | transpose -t | tail -1 | column -t) >> ${results}
(echo -n "filtered-not-used BNG"$'\t'; perl $BNG_HSSCRIPTS/calc_cmap_stats.pl ${infolder}/output/align_final/filtered_not_used_BNG.cmap | tr "=" "\t" | transpose -t | tail -1 | column -t) >> ${results}
