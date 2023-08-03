# ScaffoldImplementation

### Installation on Polaris
```bash
qsub -I -l select=1 -l walltime=1:00:00 -A datascience -q debug -l filesystems=home:eagle
module load conda/2023-01-10-unstable
conda activate
conda create -n scaffold_imp --clone base
conda activate scaffold_imp
make install
