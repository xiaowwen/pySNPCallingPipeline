                                        ABOUT
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Python implementation of Java SNP calling pipeline (https://github.com/DSGlab/SNPCallingPipeline).
Usage is descriped below. For consistency the configuration file is exactly as described in the Java code.
Descriptions of the various options in the Java code is given at https://github.com/DSGlab/SNPCallingPipeline.
An example configuration file is included in this directory.
                                       USAGE
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
usage: pySNPCallingPipeline.py [-h] -c CONF_FILE [--aln_only] [--no_aln]
                               [--local] [--submit] [--subset SUBSET]
                               [--def_run]

Run pySNPCallingPipeline

optional arguments:
  -h, --help            show this help message and exit
  -c CONF_FILE, --conf CONF_FILE
                        A configuration file is required to run
                        pySNPCallingPipeline.
  --aln_only            Alignment only, default is FALSE.
  --no_aln              Run all analysis using pre-run alignments, default is
                        FALSE.
  --local               Is this a LOCAL run or should SLURM files be created?
                        Default is True.
  --submit              If running a supercomputing cluster, should only slurm
                        files be created or should the jobs be submitted as
                        well.
  --subset SUBSET       Provide a comma seperated list of a subset of
                        "alignment, getHQSNPs, intraClonalSNPs, checkSNPs,
                        filterSNPs"
  --def_run             Default: run entire calculation locally.

