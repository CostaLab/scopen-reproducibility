import os
import subprocess

input_dir = "../../../Signac/data/Fibroblast/filtered_peak_bc_matrix/"
output_dir = "../"

job_name = "scOpen"
output_prefix = "scOpen"
subprocess.run(["sbatch", "-J", job_name,
                "-o", f"./cluster_out/{job_name}.txt",
                "-e", f"./cluster_err/{job_name}.txt",
                "--time", "120:00:00",
                "--mem", "180G",
                "-c", "6",
                "run.zsh", input_dir, output_dir, output_prefix])
