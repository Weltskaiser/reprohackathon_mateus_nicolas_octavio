<h1 align="center">🚀 ReproHackathon — How to Run the Project</h1>

In light of the challenges of reproducibility, the objective of this project was to reproduce the differential expression analysis presented in the article by Peyrusson and collaborators (2020) under the same conditions described in the body of the work, respecting—so far as possible—the data, experimental design, criteria, and steps reported in the document. The team’s work aimed to show how implementation decisions—such as pipeline structure, organization of inputs and outputs, and versions of software and libraries—impact the reproducibility and traceability of the final result, specifically in the reproduced figures.
In operational terms, an executable workflow (Nextflow) and a containerized setup (Doker) were used. The executed code was well documented and made available on Github.

Take your seat, and grab a coffee!
Execution: 4hr

Before starting, ensure you are using a machine with:

 - 16+ CPUs;
 - 64 GB RAM;
 - At least 200 GB available storage (if you have less, it might work too)

🔧 Step 1 — Clone the Repository

  `git clone <repository-url>`\
  `cd <repository-folder>`

▶️ Step 2 — Run the Pipeline

  `. run.sh`

## You can find the results inside "Results"
