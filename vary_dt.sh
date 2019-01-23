for dt in 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8; do
    python run_emcell_mfpt_sim.py funnel 0 $dt
    python run_emcell_mfpt_sim.py sharp 0 $dt
    python run_emcell_mfpt_sim.py funnel -30000 $dt
    python run_emcell_mfpt_sim.py sharp -30000 $dt
done
