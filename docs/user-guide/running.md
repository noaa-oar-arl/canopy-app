# Running the Model

This guide covers how to run the Canopy-App model in different environments and configurations.

## Basic Execution

### Single Point Simulation

For a basic single-point simulation:

```bash
# Navigate to the model directory
cd /path/to/canopy-app

# Run the model
./canopy_app.exe
```

The model will read the `namelist.canopy` file in the current directory and process the specified input files.

### Multi-Point Simulation

For simulations over multiple grid points, modify the namelist parameters:

```fortran
&CANOPY_OPTIONS
 nlat = 10
 nlon = 15
 ntime = 24
 ! ... other parameters
/
```

## Running with Different Input Formats

### NetCDF Input

```fortran
&CANOPY_OPTIONS
 infmt_opt = 0
 file_vars = 'namelist.canopy'
 ! NetCDF files will be read automatically
/
```

### Text Input

```fortran
&CANOPY_OPTIONS
 infmt_opt = 1
 file_vars = 'namelist.canopy'
 ! Text files will be read automatically
/
```

## Parallel Execution

### OpenMP

The model supports OpenMP threading for shared-memory parallelism:

```bash
# Set number of threads
export OMP_NUM_THREADS=4

# Run the model
./canopy_app.exe
```

## SLURM Job Submission

For running on HPC systems with SLURM:

```bash
#!/bin/bash
#SBATCH --job-name=canopy_app
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --output=canopy_%j.out
#SBATCH --error=canopy_%j.err

# Load required modules
module load intel/2021.2
module load netcdf/4.8.1

# Set OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the model
cd $SLURM_SUBMIT_DIR
./canopy_app.exe
```

Submit the job:

```bash
sbatch run_canopy_slurm.sh
```

## Python Interface

The model can also be run through the Python interface:

```python
import canopy_app

# Initialize the model
model = canopy_app.CanopyModel()

# Load configuration
model.load_namelist('namelist.canopy')

# Run simulation
results = model.run()

# Process results
print(f"Simulation completed with {len(results)} time steps")
```

## Performance Optimization

### Compiler Optimization

Compile with optimization flags for better performance:

```bash
make clean
make FFLAGS="-O3 -march=native -fopenmp"
```

### Grid Size Optimization

Choose grid sizes that are multiples of common factors for better cache performance:

- Good: nlat=64, nlon=128 (powers of 2)
- Good: nlat=60, nlon=120 (multiples of 12)
- Avoid: nlat=61, nlon=127 (prime numbers)

### I/O Optimization

- Use netCDF input format for large simulations
- Enable compression in output files
- Use parallel I/O for multi-process runs

## Monitoring Progress

### Real-time Monitoring

Monitor simulation progress:

```bash
# Watch output file
tail -f canopy_output.log

# Monitor memory usage
watch -n 1 'ps aux | grep canopy_app'
```

### Progress Indicators

The model provides progress information:

```
Canopy-App Model v1.0
Starting simulation...
Time step    1/24 (04.2%) - 2022-06-30 12:00:00
Time step    2/24 (08.3%) - 2022-06-30 13:00:00
...
Simulation completed successfully!
Total runtime: 00:02:35
```

## Common Runtime Issues

### Memory Issues

```
Error: Out of memory during allocation
```

Solutions:
- Reduce grid size or number of time steps
- Increase system memory or swap
- Use `ulimit -v` to check memory limits

### Input File Issues

```
Error: Cannot open input file 'input.nc'
```

Solutions:
- Check file path and permissions
- Verify file format matches `infmt_opt` setting
- Ensure all required variables are present

### Numerical Issues

```
Warning: Solution did not converge
```

Solutions:
- Reduce time step size
- Check input data for unrealistic values
- Adjust solver tolerance parameters

## Batch Processing

### Multiple Simulations

Process multiple cases:

```bash
#!/bin/bash
for case in case1 case2 case3; do
    echo "Running $case"
    cd $case
    ./canopy_app.exe > output_$case.log 2>&1
    cd ..
done
```

### Parameter Sweeps

Automate parameter sensitivity studies:

```python
import os
import subprocess

# Parameter ranges
dx_values = [50, 100, 200, 500]
lai_values = [2.0, 4.0, 6.0, 8.0]

for dx in dx_values:
    for lai in lai_values:
        # Create run directory
        run_dir = f"run_dx{dx}_lai{lai}"
        os.makedirs(run_dir, exist_ok=True)

        # Modify namelist
        # ... (modify parameters)

        # Run simulation
        subprocess.run(["./canopy_app.exe"], cwd=run_dir)
```

For more information on troubleshooting runtime issues, see the [Troubleshooting Guide](troubleshooting.md).
