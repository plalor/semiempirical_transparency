Workflow:

### Beam simulations
# first, modify run_phi.sh to submit an appropriate number of jobs
sbatch run_phi.sh phi_4MeV
sbatch run_phi.sh phi_6MeV
sbatch run_phi.sh phi_10MeV
sbatch run_phi.sh phi_9.7MeV
sbatch run_phi.sh phi_6.3MeV

# Then, analyze the output to produce beam spectra
sbatch calc_phi.sh phi_4MeV
sbatch calc_phi.sh phi_6MeV
sbatch calc_phi.sh phi_10MeV
sbatch calc_phi.sh phi_9.7MeV
sbatch calc_phi.sh phi_6.3MeV

# Then, run Bremmstrahlung.ipynb to calculate input_spectrum.txt

### Detector response simulations
# first, modify run_D.sh to submit an appropriate number of jobs
sbatch run_D.sh

# Then, analyze the output to produce beam spectra
sbatch calc_D.sh

### Open beam runs
# modify run_open_beam.sh to submit an appropriate number of jobs
sbatch run_open_beam.sh

### Perform calibration runs
# modify run_calib.sh to submit an appropriate number of jobs
sbatch run_calib.sh

# Next, merge the output
sbatch merge_files.sh

### Perform test runs
# modify run_test_materials.sh to submit an appropriate number of jobs
sbatch run_test_materials.sh

### Perform low lambda runs
# modify run_low_lmbda.sh to submit an appropriate number of jobs
sbatch run_low_lmbda.sh

### Perform high Z phantom runs
# modify run_high_Z.sh to submit an appropriate number of jobs
sbatch run_high_Z.sh
