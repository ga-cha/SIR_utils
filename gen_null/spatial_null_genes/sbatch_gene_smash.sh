#!/bin/bash
# A job manager for generating spatial null genes using brainSMASH
# This script submits jobs to the SLURM scheduler for each null index

function queue_jobs {
    while true; do
        n_pending=$(squeue -u "$USER" | grep ' PD ' | wc -l)
        if [ "$n_pending" -lt 2 ]; then  # Leave room for new array
            break
        fi
        sleep 1m
    done
}

out=sbatch_out.log
  
    sbatch --job-name="gene_smash" --output="$out" --qos=shortq --array=3-1000 <<'EOL'
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gabriella.chan@monash.edu

# define python executable and input file inside the compute script so they are available at runtime
BASE_DIR="${SLURM_SUBMIT_DIR}"
py_sir=/fs03/kg98/gchan/conda_envs/sir/bin/python
in_file="${BASE_DIR}/data/gene_expr_s132_indiv.csv"

IDX=${SLURM_ARRAY_TASK_ID}
echo "Processing gene null index: ${IDX}"

results_dir="${BASE_DIR}/results/gene_null"
results_file="${results_dir}/null_genes_${IDX}.mat"
mkdir -p "${results_dir}"

# run the Python generator
"$py_sir" "${BASE_DIR}/gene_brainsmash.py" "${IDX}" "${in_file}"

# convert the saved matrix into a MATLAB table with variable names from the CSV header
MATLAB_MFILE="/tmp/convert_null_${IDX}.m"
cat > "$MATLAB_MFILE" <<MML
try
    tbl = readtable('${in_file}', 'VariableNamingRule', 'preserve');
    gene_names = tbl.Properties.VariableNames;
    data = load('${results_file}');
    null_slice = data.null_genes;
    null_genes = array2table(null_slice, 'VariableNames', gene_names);
    save('${results_file}', 'null_genes');
catch ME
    disp(getReport(ME));
    exit(1);
end
exit(0);
MML

matlab -nodisplay -nosplash -r "run('$MATLAB_MFILE');"

EOL
