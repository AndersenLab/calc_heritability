#! usr/bin/env nextflow

if( !nextflow.version.matches('20+') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )
// reps = 10000
params.reps = 500
params.binDir = "${workflow.projectDir}"
params.species = "c_elegans"
params.fix = "fix"
params.maf = 0.05

/*
~ ~ ~ > * Parameters
*/

download_vcf = null
params.R_libpath = ""

// VCF param
if(params.debug) {
    println """
        *** Using debug mode ***
    """
    // debug for now with small vcf
    params.vcf = "WI.20220216.small.hard-filter.isotype.vcf.gz"

    vcf_file = Channel.fromPath("/projects/b1059/data/${params.species}/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz")
    vcf_index = Channel.fromPath("/projects/b1059/data/${params.species}/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz.tbi")
    params.traitfile = "${params.binDir}/input_data/testdata/ExampleTraitData.tsv"

    // lower number of reps for debug
    // reps = 10
        
} else if(params.gcp) { 
    // use the data directly from google on gcp
    // vcf_file = Channel.fromPath("gs://elegansvariation.org/releases/20220216/variation/WI.20220216.small.hard-filter.isotype.vcf.gz")
    // vcf_index = Channel.fromPath("gs://elegansvariation.org/releases/20220216/variation/WI.20220216.small.hard-filter.isotype.vcf.gz.tbi")

    vcf_file = Channel.fromPath("gs://caendr-site-public-bucket/dataset_release/c_elegans/20220216/variation//WI.20220216.small.hard-filter.isotype.vcf.gz")
    vcf_index = Channel.fromPath("gs://caendr-site-public-bucket/dataset_release/c_elegans/20220216/variation//WI.20220216.small.hard-filter.isotype.vcf.gz.tbi")

} else if(!params.vcf) {
    // if there is no VCF date provided, pull the latest vcf from cendr.
    params.vcf = "20220216"
    download_vcf = true
    
} else {
    // use the vcf data from QUEST when a cendr date is provided

    // Check that params.vcf is valid
    if("${params.vcf}" == "20220216" || "${params.vcf}" == "20210121" || "${params.vcf}" == "20200815" || "${params.vcf}" == "20180527" || "${params.vcf}" == "20170531" || "${params.vcf}" == "20210901" || "${params.vcf}" == "20210803") {
        // if("${params.vcf}" in ["20210121", "20200815", "20180527", "20170531", "20210901"]) {
        // check to make sure 20210901 is tropicalis
        if("${params.vcf}" == "20210901") {
            if("${params.species}" == "c_elegans" || "${params.species}" == "c_briggsae") {
                println """
                Error: VCF file (${params.vcf}) does not match species ${params.species} (should be c_tropicalis). Please enter a new vcf date or a new species to continue.
                """
                System.exit(1)
            }
        }
        // check to make sure vcf matches species for briggsae
        if("${params.vcf}" == "20210803") {
            if("${params.species}" == "c_elegans" || "${params.species}" == "c_tropicalis") {
                println """
                Error: VCF file (${params.vcf}) does not match species ${params.species} (should be c_briggsae). Please enter a new vcf date or a new species to continue.
                """
                System.exit(1)
            }
        }
        // check to make sure vcf matches species for elegans
        if("${params.vcf}" == "20220216" || "${params.vcf}" == "20210121" || "${params.vcf}" == "20200815" || "${params.vcf}" == "20180527" || "${params.vcf}" == "20170531") {
            if("${params.species}" == "c_briggsae" || "${params.species}" == "c_tropicalis") {
                println """
                Error: VCF file (${params.vcf}) does not match species ${params.species} (should be c_elegans). Please enter a new vcf date or a new species to continue.
                """
                System.exit(1)
            }
        }
        // use the vcf data from QUEST when a cendr date is provided
        vcf_file = Channel.fromPath("/projects/b1059/data/${params.species}/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz")
        vcf_index = Channel.fromPath("/projects/b1059/data/${params.species}/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz.tbi")
    } else {
        // check that vcf file exists, if it does, use it. If not, throw error
        if (!file("${params.vcf}").exists()) {
            println """
            Error: VCF file (${params.vcf}) does not exist. Please provide a valid filepath or a valid CeNDR release date (i.e. 20210121)
            """
            System.exit(1)
        } else {
            // if it DOES exist
            println """
            WARNING: Using a non-CeNDR VCF for analysis. 
            """
            vcf_file = Channel.fromPath("${params.vcf}")
            vcf_index = Channel.fromPath("${params.vcf}.tbi")

        }
    }
}



/*
~ ~ ~ > * WORKFLOW
*/

workflow {

	// if no VCF is provided, download the latest version from CeNDR
	if(download_vcf) {
	    pull_vcf()

	    vcf_file = pull_vcf.out.hard_vcf
	    vcf_index = pull_vcf.out.hard_vcf_index
	}

	// Fix strain names
    // allow for non c_elegans, briggsae, tropicalis species
    if(!("${params.species}" == "c_elegans" || "${params.species}" == "c_tropicalis" || "${params.species}" == "c_briggsae")) {
        Channel.fromPath("${params.traitfile}")
            .combine(Channel.fromPath("${params.binDir}/input_data/${params.species}/strain_isotype_lookup.tsv"))
            .combine(Channel.fromPath("${params.binDir}/bin/Fix_Isotype_names_bulk_h2.R"))
            .combine(Channel.from("false")) | fix_strain_names_bulk 
    } else {
        Channel.fromPath("${params.traitfile}")
            .combine(Channel.fromPath("${params.binDir}/input_data/${params.species}/strain_isotype_lookup.tsv"))
            .combine(Channel.fromPath("${params.binDir}/bin/Fix_Isotype_names_bulk_h2.R"))
            .combine(Channel.from("${params.fix}")) | fix_strain_names_bulk 
    }
           
    traits_to_map = fix_strain_names_bulk.out.fixed_strain_phenotypes
            .flatten()
            .map { file -> tuple(file.baseName.replaceAll(/pr_/,""), file) }

    // Genotype matrix
    pheno_strains = fix_strain_names_bulk.out.phenotyped_strains_to_analyze

    vcf_file.combine(vcf_index)
            .combine(pheno_strains) | vcf_to_geno_matrix

    // calclate heritability and generate report output
    traits_to_map 
    	.combine(vcf_to_geno_matrix.out)
        .combine(Channel.from("${params.reps}"))
        .combine(Channel.fromPath("${params.binDir}/bin/20210716_H2_script.R")) | heritability

    //generate html report
    traits_to_map
    	.join(heritability.out)
    	.combine(fix_strain_names_bulk.out.strain_issues)
        .combine(Channel.fromPath("${params.binDir}/bin/20210716_hert_report.Rmd")) | html_report

}

/*
==============================================
~ > *                                    * < ~
~ ~ > *                                * < ~ ~
~ ~ ~ > *  DOWNLOAD VCF FROM CENDR   * < ~ ~ ~
~ ~ > *                                * < ~ ~
~ > *                                    * < ~
==============================================
*/

process pull_vcf {

    tag {"PULLING VCF FROM CeNDR"}

    output:
        path "*hard-filter.isotype.vcf.gz", emit: hard_vcf 
        path "*hard-filter.isotype.vcf.gz.tbi", emit: hard_vcf_index 
        path "*impute.isotype.vcf.gz", emit: impute_vcf 
        path "*impute.isotype.vcf.gz.tbi", emit: impute_vcf_index 
        path "*.strain-annotation*.tsv", emit: ann_vcf

    """
        wget https://storage.googleapis.com/elegansvariation.org/releases/${params.vcf}/variation/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz
        tabix -p vcf WI.${params.vcf}.small.hard-filter.isotype.vcf.gz

        wget https://storage.googleapis.com/elegansvariation.org/releases/${params.vcf}/variation/WI.${params.vcf}.impute.isotype.vcf.gz
        tabix -p vcf WI.${params.vcf}.impute.isotype.vcf.gz

        wget https://storage.googleapis.com/elegansvariation.org/releases/${params.vcf}/variation/WI.${params.vcf}.strain-annotation.bcsq.tsv

    """
}

/*
==============================================================
~ > *                                                    * < ~
~ ~ > *                                                * < ~ ~
~ ~ ~ > *  FIX STRAIN NAMES TO MATCH THOSE ON CENDR  * < ~ ~ ~
~ ~ > *                                                * < ~ ~
~ > *                                                    * < ~
==============================================================
*/

/*
THIS WILL NEED TO BE UPDATED TO HANDLE OTHER SPECIES
*/


process fix_strain_names_bulk {

    tag {"BULK TRAIT"}

    publishDir "${params.out}/Phenotypes", mode: 'copy', pattern: "*pr_*.tsv"
    publishDir "${params.out}/Phenotypes", mode: 'copy', pattern: "strain_issues.txt"

    input:
        tuple file(phenotypes), file(isotype_lookup), file(fix_isotype_script), val(fix)

    output:
        path "pr_*.tsv", emit: fixed_strain_phenotypes 
        path "Phenotyped_Strains.txt", emit: phenotyped_strains_to_analyze 
        path "strain_issues.txt", emit: strain_issues

    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${fix_isotype_script} > Fix_Isotype_names_bulk.R 
        Rscript --vanilla Fix_Isotype_names_bulk.R ${phenotypes} $fix $isotype_lookup
    """

}

/*
===========================================================
~ > *                                                 * < ~
~ ~ > *                                             * < ~ ~
~ ~ ~ > *  CONVERT THE VCF TO A GENOTYPE MATRIX   * < ~ ~ ~
~ ~ > *                                             * < ~ ~
~ > *                                                 * < ~
===========================================================
*/

process vcf_to_geno_matrix {

    publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

    label "large"

    input:
        tuple file(vcf), file(index), file(strains)

    output:
        file("Genotype_Matrix.tsv") 

    """
        bcftools view -S ${strains} -Ou ${vcf} |\\
        bcftools filter -i N_MISSING=0 -Oz --threads 5 -o Phenotyped_Strain_VCF.vcf.gz
        tabix -p vcf Phenotyped_Strain_VCF.vcf.gz
        plink --vcf Phenotyped_Strain_VCF.vcf.gz \\
              --threads 5 \\
              --snps-only \\
              --biallelic-only \\
              --maf ${params.maf} \\
              --set-missing-var-ids @:# \\
              --indep-pairwise 50 10 0.8 \\
              --geno \\
              --allow-extra-chr
        awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
        sort -k1,1d -k2,2n > markers.txt
        bcftools query -l Phenotyped_Strain_VCF.vcf.gz |\\
        sort > sorted_samples.txt 
        bcftools view -v snps \\
        -S sorted_samples.txt \\
        -R markers.txt -Ou \\
        Phenotyped_Strain_VCF.vcf.gz |\\
        bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' |\\
        sed 's/[[# 0-9]*]//g' |\\
        sed 's/:GT//g' |\\
        sed 's/0|0/-1/g' |\\
        sed 's/1|1/1/g' |\\
        sed 's/0|1/NA/g' |\\
        sed 's/1|0/NA/g' |\\
        sed 's/.|./NA/g'  |\\
        sed 's/0\\/0/-1/g' |\\
        sed 's/1\\/1/1/g'  |\\
        sed 's/0\\/1/NA/g' |\\
        sed 's/1\\/0/NA/g' |\\
        sed 's/.\\/./NA/g' > Genotype_Matrix.tsv
    """

}


/*
===========================================================
~ > *                                                 * < ~
~ ~ > *                                             * < ~ ~
~ ~ ~ > * CALCULATE BROAD AND NARROW HERITABILITY * < ~ ~ ~
~ ~ > *                                             * < ~ ~
~ > *                                                 * < ~
===========================================================
*/

process heritability {

    label "medium"

	publishDir "${params.out}/", mode: 'copy'

	input:
		tuple val(TRAIT), file(phenotype), file(geno_matrix), val(reps), file(h2_script)


	output:
		tuple val(TRAIT), file("heritability_result.tsv")


	"""
		# add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${h2_script} > H2_script.R 
        Rscript --vanilla H2_script.R ${phenotype} ${geno_matrix} ${reps}

	"""


}


/*
===========================================================
~ > *                                                 * < ~
~ ~ > *                                             * < ~ ~
~ ~ ~ > *           GENERATE HTML REPORT          * < ~ ~ ~
~ ~ > *                                             * < ~ ~
~ > *                                                 * < ~
===========================================================
*/

process html_report {

    label "small"

	publishDir "${params.out}/", mode: 'copy'

	input:
		tuple val(TRAIT), file(phenotype), file(hert), file(strain_issues), file(html_report)

	output:
		tuple file("*.html"), file("h2_plot.png")


	"""
		cat "${html_report}" | \\
		sed 's+Phenotypes/pr_TRAITNAME.tsv+${phenotype}+g' | \\
		sed "s+TRAITNAME+${TRAIT}+g" | \\
		sed 's+heritability_result.tsv+${hert}+g' | \\
		sed 's+Phenotypes/strain_issues.txt+${strain_issues}+g' > hert_report_${TRAIT}.Rmd 
	    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" > .Rprofile
	    Rscript -e "rmarkdown::render('hert_report_${TRAIT}.Rmd')"

	"""


}






