-> Run the following bash codes from the `containers` folder to build the required images:
```sh
docker build -f Dockerfile.fasterq_dump -t fasterq_dump .
docker build -f Dockerfile.trim_galore_cutadapt -t trim_galore_cutadapt:1.11 .
docker build -f Dockerfile.bowtie -t bowtie:0.12.7 .
docker build -f Dockerfile.feature_counts -t feature_counts:1.4.6-p3 .
docker build -f Dockerfile.deseq2_r -t deseq2_r .
```
