# calc_heritability

This pipeline calculates both broad sense and narrow sense heritability and is used as the heritability calculator for [CeNDR](elegansvariation.org).

## Usage on QUEST

```
# set up environment
module load python/anaconda3.6
source activate /projects/b1059/software/conda_envs/nf20_env
module add singularity

# run pipeline
nextflow run andersenlab/calc_heritability --vcf 20210121 --traitfile <path_to_file>
```

## Usage on local machine (not QUEST)

*Note: this pipeline requires previous installation of [docker](https://docs.docker.com/get-docker/) and [nextflow](https://www.nextflow.io/docs/latest/getstarted.html)*

* By failing to specify a vcf, the pipeline will automatically download the current *C. elegans* release from CeNDR

```
# run pipeline
nextflow run andersenlab/calc_heritability --traitfile <path_to_file>
```


## Usage on GCP/CeNDR

Provision a virtual machine from the existing template and connect with SSH, then switch to the root user and install the required packages:

```
sudo su
apt-get install docker.io git nano
```

Use 'nano' to create a *.env file (see: example.env) with variables pointing to your Google Storage locations:

```
nano test.env
```

test.env:
```
TRAIT_FILE="gs://elegansvariation.org/reports/heritability/kse_test_hert/ExampleTraitData.csv"
OUTPUT_DIR="gs://elegansvariation.org/reports/heritability/kse_test_hert2/results"
WORK_DIR="gs://nf-pipelines/workdir/kse_test_hert"
VCF_VERSION="20210121"
```

```
docker run -i -t \
  --env-file test.env \
  andersenlab/calc_heritability:v0.1 \
  heritability-nxf.sh
```

You can also pass them in as part of the command:
```
docker run -i -t \
  -e TRAIT_FILE="gs://elegansvariation.org/reports/heritability/kse_test_hert/ExampleTraitData.csv" \
  -e OUTPUT_DIR="gs://elegansvariation.org/reports/heritability/kse_test_hert2/results" \
  -e WORK_DIR="gs://nf-pipelines/workdir/kse_test_hert" \
  -e VCF_VERSION="20210121" \
  andersenlab/calc_heritability:v0.1 \
  heritability-nxf.sh
```