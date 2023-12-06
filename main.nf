#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import processes from src:
include {
    mixcr_align_bulk; mixcr_assemblePartial_1; mixcr_assemblePartial_2;
    mixcr_assemble; mixcr_assembleContigs; mixcr_exportClones; 
    mixcr_align_10X_5p; mixcr_refineTags; platypus
} from './src/processes.nf'

// Print pipeline parameters to stdout:
log.info """\
RNASeq - BCR Pipeline

===================================

This pipeline was developed by Brent Schlegel, Analyst at the Bioinformatics Core @ the UPMC Children's Hospital of Pittsburgh (2023).
It was inspired by the work of Rubio et. al (A Nextflow pipeline for T-cell receptor repertoire reconstruction and analysis from RNA sequencing data,
ImmunoInformatics,2022) and their subsequent MiXCR-based pipeline for extracting TCR repertoires from bulk RNA-seq data (https://github.com/ConesaLab/TCR_nextflow).

The use of MiXCR and IGMT in this pipeline are for strictly non-commerical, academic purposes.
In order to obtain a non-commercial license to use MiXCR, please see: https://licensing.milaboratories.com/

===================================

Project parameters:
- Project Name          : ${params.project_name}
- Sample List           : ${params.readsfile}
- IGMT Library          : ${params.igmt}
- Species               : ${params.species}
- Chain                 : ${params.chain}
- Input Directory       : ${params.indir}
- Output Directory      : ${params.outdir}
- Single Cell?          : ${params.is_sc}

"""
// Bulk RNAseq Workflow:
workflow mixcr_bulk {

    take: samples_channel
    main:
        mixcr_align_bulk(samples_channel) | 
        mixcr_assemblePartial_1 | 
        mixcr_assemblePartial_2 | 
        mixcr_assemble | 
        mixcr_assembleContigs | 
        mixcr_exportClones
}
// 10X 5' GEX Workflow:
workflow mixcr_10X_GEX {

    take: samples_channel
    main:
        mixcr_align_10X_5p(samples_channel) | 
        mixcr_assemblePartial_1 | 
        mixcr_assemblePartial_2 |
        mixcr_assemble |
        mixcr_assembleContigs |
        mixcr_exportClones
}

workflow {

    // Input validation
    def valid_species = ['hsa','mmu'] // pipeline only supports Human and Mouse samples
    is_valid_specie = params.species in valid_species
    if (!is_valid_specie) {
        log.error "`params.species` must be one of ${valid_species}"
        exit 1,   "`params.species` must be one of ${valid_species}"
    }

    Channel.fromPath("${params.readsfile}") // check that the master CSV file exists, pass to Input channel
        .ifEmpty { exit 1, "File not foud: ${params.readsfile}" }
        .set { sampleInfoChannel }

   /*
   * Create a channel that emits tuples containing three elements:
   * the SampleID, the first read-pair file and the second read-pair file
   */
    samples_channel = sampleInfoChannel
        .splitCsv(header:true)                  // Read in the CSV
        .map { row -> tuple (row.SampleID,      // Map to tuple, Sample ID
        file(row.R1.trim()),                    // R1 FASTQ
        file(row.R2.trim()),                    // R2 FASTQ
        file(params.igmt))                      // IGMT database file (MiXCR allele reference)
    }

    /*
    * Create a channel that emits a tuple containing four elements:
    * the location of the MiXCR output, the sample info file,
    * the downstream Platypus script, and a script containing utility functions 
    * for processing the MiXCR results
    */
    Channel.fromPath("${params.outdir}/MiXCR")                  // Input directory
      .combine(Channel.fromPath("${params.readsfile}"))         // Reads file (.csv)
      .combine(Channel.fromPath("${baseDir}/src/platypus.R"))   // Path to platypus.R
      .combine(Channel.fromPath("${baseDir}/src/utils.R"))      // Path to utils.R
      .map { outdir, readsFile, platypusScript, utilsScript ->  // Map to tuple
          tuple(outdir, readsFile, platypusScript, utilsScript)
      }                                                         
    .set { platypus_channel }
    
    // If the sample is 10X 5' GEX, use single-cell MiXCR preset, else use the bulk preset:
    (params.is_sc ? mixcr_10X_GEX(samples_channel) : mixcr_bulk(samples_channel))

    // Downstream processing w/ Platypus:
    platypus(platypus_channel)
}
