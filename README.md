# semiempirical_transparency

Semiempirical_transparency is a repository for running radiography simulations and comparing them to different analytic and semiempirical models. This analysis reveals that the free streaming model shows a clear bias compared to the simulated results. The semiempirical model shows excellent agreement with experimental data.

### To run:
```console
cd run
sh run.sh calib 20
sbatch calcEnergyDeposited.sh calib
```
