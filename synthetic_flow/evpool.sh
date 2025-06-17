#!/usr/bin/env bash
#JDL_OUTPUT=evtpool.root@disk=2,*stat@disk=1,*.json@disk=1
#JDL_PACKAGE=O2sim::v20241217-1

# the ini file is the only thing that should be given; it needs to be prepared from the PWGs
INIFILE=${O2DPG_ROOT}/MC/config/PWGLF/ini/GeneratorLFHyperNucleiPbPbGapWithFlow.ini
NEV=1
NWORKERS=8

# we prepare a generator config
${O2DPG_ROOT}/MC/bin/o2_hybrid_gen.py --iniFile ${INIFILE} --clone ${NWORKERS} --mode parallel --output gen_config.json

# nothing needs to be changed here; it's some standard configs; The seed is automatically picked up from ALIEN_PROC_ID
${O2DPG_ROOT}/MC/bin/o2dpg_sim_workflow.py -gen hybrid -eCM 1.0 -tf 1 -ns ${NEV} --make-evtpool -confKey "GeneratorHybrid.configFile=${PWD}/gen_config.json;GeneratorHybrid.num_workers=${NWORKERS}" -interactionRate 1000000 -seed ${ALIEN_PROC_ID:-11} -run 300100

# we launch the event pool production; It will automatically timeout gracefully when JOBTTL is found
${O2DPG_ROOT}/MC/bin/o2dpg_workflow_runner.py -f workflow.json -tt pool
