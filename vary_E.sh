for ez in -100000 -80000 -60000 -40000 -20000 -10000 -1000 -100 0 100 1000 10000 20000 30000; do
    # python run_emcell_mfpt_sim.py funnel $ez 1e-7
    python run_emcell_mfpt_sim.py sharp $ez 1e-8
done
