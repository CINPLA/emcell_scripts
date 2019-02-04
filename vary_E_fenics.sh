#!/bin/bash
for ez in -100000 -80000 -60000 -40000 -20000 -10000 -1000 -100 0 100 1000 10000 20000 30000; do
    python run_fenics_mfpt_sim.py funnel $ez
    python run_fenics_mfpt_sim.py sharp $ez
done
# COUNTER=-120000
# while [  $COUNTER -lt 50001 ]; do
#     python run_fenics_mfpt_sim.py funnel $COUNTER
#     python run_fenics_mfpt_sim.py sharp $COUNTER
#     let COUNTER=COUNTER+500
# done
