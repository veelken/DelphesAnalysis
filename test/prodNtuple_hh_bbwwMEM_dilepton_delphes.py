#!/usr/bin/env python
import os, logging, sys, getpass
from collections import OrderedDict as OD
from hhAnalysis.DelphesAnalysis.configs.prodNtupleConfig_hh_bbwwMEM_dilepton_delphes import prodNtupleConfig_hh_bbwwMEM_dilepton_delphes
from tthAnalysis.HiggsToTauTau.jobTools import query_yes_no
from tthAnalysis.HiggsToTauTau.runConfig import tthAnalyzeParser

# E.g.: ./prodNtuple_hh_bbwwMEM_dilepton_delphes.py -e 2016 -v 2021Nov22v2

parser = tthAnalyzeParser()
parser.add_use_home()
args = parser.parse_args()

# Common arguments
era                = args.era
version            = args.version
no_exec            = args.no_exec
auto_exec          = args.auto_exec
check_output_files = not args.not_check_input_files
debug              = args.debug
num_parallel_jobs  = args.num_parallel_jobs
running_method     = args.running_method

# Additional arguments
use_home           = args.use_home

samples = {
  'signal_lo'     : "/hdfs/local/karl/DelphesNtuples/2016/HH_DL_LO_PU40/0000/",
  'background_lo' : "/hdfs/local/karl/DelphesNtuples/2016/ttbar_DL_LO_PU40/0000/"
}

if __name__ == '__main__':
  logging.basicConfig(
    stream = sys.stdout,
    level  = logging.INFO,
    format = '%(asctime)s - %(levelname)s: %(message)s',
  )

  analysis = prodNtupleConfig_hh_bbwwMEM_dilepton_delphes(
    configDir           = os.path.join("/home",       getpass.getuser(), "hhAnalysis", era, version),
    outputDir           = os.path.join("/hdfs/local", getpass.getuser(), "hhAnalysis", era, version),
    executable_analyze  = "produceMEMNtuple_hh_bb2l_delphes",
    cfgFile_analyze     = "produceMEMNtuple_hh_bb2l_delphes_cfg.py",
    samples             = samples,
    ##max_jobs_per_sample = 1000, # CV: use for tests
    max_jobs_per_sample = 5000, # CV: use when making plots for paper
    max_events_per_job  = 500,
    era                 = era,
    check_output_files  = check_output_files,
    running_method      = running_method,
    num_parallel_jobs   = num_parallel_jobs,
    isDebug             = debug,
    use_home            = use_home,
  )

  job_statistics = analysis.create()
  for job_type, num_jobs in job_statistics.items():
    logging.info(" #jobs of type '%s' = %i" % (job_type, num_jobs))

  if auto_exec:
    run_analysis = True
  elif no_exec:
    run_analysis = False
  else:
    run_analysis = query_yes_no("Start jobs ?")
  if run_analysis:
    analysis.run()
  else:
    sys.exit(0)
