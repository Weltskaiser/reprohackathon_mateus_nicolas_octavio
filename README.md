Run ReproHackathon Project

To run the ReproHackathon project, start by cloning the repository on a machine that has at least 16 CPUs, 32 GB RAM, and 200 GB of storage:

git clone <repository-url>
cd <repository-folder>


Next, activate your Conda environment configured for Nextflow:

conda activate nextflow


Finally, run the main Nextflow pipeline:

nextflow run main.nf
