#!/usr/bin/python
# need to add /opt/delft3dfm/lib to ld_path for successful run
import bmi.wrapper as bw


dflowfm = bw.BMIWrapper(engine='dflowfm', configfile='../tilted-v_empty-run.mdu')

dflowfm.initialize()

dflowfm.update()

dflowfm.finalize()