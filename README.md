![ngstools](pictures/Irys_icon.png) - BioNano-Tools
bionano-tools
==========

# QC-tools

Tools to process and QC BioNanoGenomics data.

## **run_MQR.sh**

The bash script **[run_MQR.sh](qc-tools/run_MQR.sh)**

Perform molecule quality report (MQR) at CLI instead of running this under IrysView. One may prefer to perform the MQR directly on his/her Nix server. The main advantage is that one can launch this code in a bash loop and perform all MQR from a list of BNX files (single runs) without supervision.

The following loop will process all BNX files in the current folder and create a MQR (in its own folder) from each BNX and a common reference (more parameters are available)

<pre>
for b in *.bnx; do
  run_MQR.sh -i $b -r myreference.cmap;
done
</pre>

Type the script name followed by -h will list all available parameters

**run_MQR.sh -h**
<div>
# Usage: runMQR.sh -i <molecules.bnx> -r <reference.cmap>
#		[optional: -l <minlen|150> -p <pval|1e-9>]
#		[optional: -t <max-threads|32> -m <max-ram|64>]
#		[optional: -s <sample N molecules>]
</div>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

------------

![Creative Commons License](http://i.creativecommons.org/l/by-sa/3.0/88x31.png?raw=true)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/).
