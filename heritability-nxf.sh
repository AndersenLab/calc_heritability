#!/bin/bash
#
# This script acts as a wrapper around the execution of Nextflow passing environment variables as arguments
#
###################################################################################################################

DEFAULT_VCF_VERSION="20220216"
DEFAULT_SPECIES="c_elegans"
DEFAULT_GOOGLE_PROJECT="andersen-lab"
DEFAULT_QUEUE_REGION="us-east4"
DEFAULT_DATA_BUCKET="caendr-site-private"
DEFAULT_GOOGLE_SERVICE_ACCOUNT_EMAIL="nscalc-201573431837@andersen-lab.iam.gserviceaccount.com"
# Environment variables with default values:
export GOOGLE_ZONE=us-east4-a

if [[ -z "${VCF_VERSION}" ]]; then
  VCF_VERSION=${DEFAULT_VCF_VERSION}
  echo "VCF_VERSION environment variable is not set - defaulting to ${VCF_VERSION}"
fi

if [[ -z "${SPECIES}" ]]; then
  SPECIES=${DEFAULT_SPECIES}
  echo "SPECIES environment variable is not set - defaulting to ${DEFAULT_SPECIES}"
fi

if [[ -z "${GOOGLE_PROJECT}" ]]; then
  GOOGLE_PROJECT=${DEFAULT_GOOGLE_PROJECT}
  echo "GOOGLE_PROJECT environment variable is not set - defaulting to ${GOOGLE_PROJECT}"
fi

if [[ -z "${QUEUE_REGION}" ]]; then
  QUEUE_REGION=${DEFAULT_QUEUE_REGION}
  echo "QUEUE_REGION environment variable is not set - defaulting to ${QUEUE_REGION}"
fi

if [[ -z "${DATA_BUCKET}" ]]; then
  DATA_BUCKET=${DEFAULT_DATA_BUCKET}
  echo "DATA_BUCKET environment variable is not set - defaulting to ${DEFAULT_DATA_BUCKET}"
fi

if [[ -z "${GOOGLE_SERVICE_ACCOUNT_EMAIL}" ]]; then
  GOOGLE_SERVICE_ACCOUNT_EMAIL=${DEFAULT_GOOGLE_SERVICE_ACCOUNT_EMAIL}
  echo "GOOGLE_SERVICE_ACCOUNT_EMAIL environment variable is not set - defaulting to ${GOOGLE_SERVICE_ACCOUNT_EMAIL}"
fi

# Environment variables that MUST be set

if [[ -z "${TRAIT_FILE}" ]]; then
  echo "TRAIT_FILE environment variable must be set to the Google Storage path of the data"
  exit 1
fi

if [[ -z "${OUTPUT_DIR}" ]]; then
  echo "OUTPUT_DIR environment variable must be set to the Google Storage path of the output directory"
  exit 1
fi

if [[ -z "${WORK_DIR}" ]]; then
  echo "WORK_DIR environment variable must be set to the Google Storage path of the working directory"
  exit 1
fi

  
echo "profile:                      gcp"
echo "google_project:               ${GOOGLE_PROJECT}"
echo "google_zone:                  ${QUEUE_REGION}"
echo "google_service_account_email: ${GOOGLE_SERVICE_ACCOUNT_EMAIL}"
echo "data_bucket:                  ${DATA_BUCKET}"
echo "traitfile:                    ${TRAIT_FILE}"
echo "vcf:                          ${VCF_VERSION}"
echo "species:                      ${SPECIES}"
echo "work_dir:                     ${WORK_DIR}"
echo "out:                          ${OUTPUT_DIR}"

nextflow run main.nf \
  -profile gcp \
  --google_project "${GOOGLE_PROJECT}" \
  --google_zone "${QUEUE_REGION}" \
  --google_service_account_email "${GOOGLE_SERVICE_ACCOUNT_EMAIL}" \
  --traitfile "${TRAIT_FILE}" \
  --data_bucket "${DATA_BUCKET}" \
  --vcf "${VCF_VERSION}" \
  --species "${SPECIES}" \
  --work_dir "${WORK_DIR}" \
  --out "${OUTPUT_DIR}"
