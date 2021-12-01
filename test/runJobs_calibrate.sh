#!/bin/bash

rm calibrate_jets_delphes_signal.log
calibrate_jets_delphes calibrate_jets_delphes_signal_cfg.py >& calibrate_jets_delphes_signal.log

rm calibrate_jets_delphes_background.log
calibrate_jets_delphes calibrate_jets_delphes_background_cfg.py >& calibrate_jets_delphes_background.log

rm calibrate_jets_delphes_all.root
hadd calibrate_jets_delphes_all.root calibrate_jets_delphes_signal.root calibrate_jets_delphes_background.root