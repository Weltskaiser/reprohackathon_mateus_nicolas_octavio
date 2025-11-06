Run the following bash codes to build the required images:
  docker build -f Dockerfile.sratoolkit sratoolkit
  docker build -f Dockerfile.cutadapt cutadapt .
  docker build -f Dockerfile.bowtie bowtie .
  docker build -f Dockerfile.featurecounts featurecounts .
  docker build -f Dockerfile.deseq2 deseq2 .

