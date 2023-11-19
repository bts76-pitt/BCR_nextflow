process mixcr_align_bulk{
    cpus "$params.cpus"
    time "$params.time"
    memory "$params.memory"

    tag "$SampleID"
    label "mixcr"

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID), file(R1), file(R2), file(igmt)

    output:
    tuple val(SampleID),file(igmt), path("${SampleID}.vdjca"), emit: vdjca

    script:
    """
    mixcr align \
    --library imgt.202214-2.sv8 \
    --preset rna-seq \
    --species $params.species \
    -OvParameters.geneFeatureToAlign="VTranscriptWithout5UTRWithP" \
    -OvParameters.parameters.floatingLeftBound=false \
    -OjParameters.parameters.floatingRightBound=false \
    -OcParameters.parameters.floatingRightBound=false \
    -OallowPartialAlignments=true \
    --threads $params.threads \
    --force-overwrite \
    --report ${SampleID}.report.txt \
    --json-report ${SampleID}.report.json \
    ${R1} ${R2} ${SampleID}.vdjca
    """
}

process mixcr_align_sc{
    cpus "$params.cpus"
    memory "$params.memory"
    time "$params.time"

    tag "$SampleID"
    label "mixcr"

    publishDir "$params.outdir/MiXCR",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID), file(R1), file(R2), file(igmt)

    output:
    tuple val(SampleID),file(igmt), path("${SampleID}.vdjca"), emit: vdjca

    script:
    """
    mixcr align \
    --preset 10x-sc-5gex \
    --library imgt.202214-2.sv8 \
    --species $params.species \
    -OvParameters.geneFeatureToAlign="VTranscriptWithout5UTRWithP" \
    -OvParameters.parameters.floatingLeftBound=false \
    -OjParameters.parameters.floatingRightBound=false \
    -OcParameters.parameters.floatingRightBound=false \
    -OallowPartialAlignments=true \
    --threads $params.threads \
    --force-overwrite \
    --report ${SampleID}.report.txt \
    --json-report ${SampleID}.report.json \
    ${R1} ${R2} \
    ${SampleID}.vdjca
    """
}

process mixcr_assemblePartial_1{
    cpus "$params.cpus"
    memory "$params.memory"
    time "$params.time"

    tag "$SampleID"
    label "mixcr"

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID), file(igmt), path(vdjca)

    output:
    tuple val(SampleID), file(igmt), path("${SampleID}.passembled.1.vdjca"), emit: vdjca_1

    script:
    """
    #assemble overlapping fragmented sequencing reads:
    ## First round:
    mixcr assemblePartial \
      --report ${SampleID}.report.txt \
      --json-report ${SampleID}.report.json \
      --force-overwrite \
      ${vdjca} \
      ${SampleID}.passembled.1.vdjca
    """
}

process mixcr_assemblePartial_2{
    cpus "$params.cpus"
    memory "$params.memory"
    time "$params.time"

    tag "$SampleID"
    label "mixcr"

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID), file(igmt), path(vdjca_1)

    output:
    tuple val(SampleID), file(igmt), path("${SampleID}.passembled.2.vdjca"), emit: vdjca_2

    script:
    """
    #assemble overlapping fragmented sequencing reads:
    ## Second round
    mixcr assemblePartial \
      --report ${SampleID}.report.txt \
      --json-report ${SampleID}.report.json \
      --force-overwrite \
      ${vdjca_1} \
      ${SampleID}.passembled.2.vdjca
    """
}

process mixcr_assemble{
       cpus "$params.cpus"
    memory "$params.memory"
    time "$params.time"

    tag "$SampleID"
    label "mixcr"

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID), file(igmt), path(vdjca_2)

    output:
    tuple val(SampleID), file(igmt), path("${SampleID}.clna"), emit: clna

    script:
    """
    #Full assemble:
    mixcr assemble \
        -OassemblingFeatures="CDR3" \
        -OseparateByV=true \
        -OseparateByJ=true \
        -OseparateByC=true \
        --write-alignments \
        --force-overwrite \
        --report ${SampleID}.assembled.report.txt \
        --json-report ${SampleID}.assembled.json \
        ${vdjca_2} \
        ${SampleID}.clna
    """
}

process mixcr_assembleContigs {
    cpus "$params.cpus"
    memory "$params.memory"
    time "$params.time"

    tag "$SampleID"
    label "mixcr"

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID),file(igmt), path(clna)

    output:
    tuple val(SampleID),file(igmt), path("${SampleID}.clns"), emit: clns

    script:
    """
    #Assemble contigs based on VDJ region:
    mixcr assembleContigs \
        -OsubCloningRegions=VDJRegion \
        --report ${SampleID}.report.txt \
        --json-report ${SampleID}.report.json \
        --force-overwrite \
        ${clna} \
        ${SampleID}.clns
    """
}

process mixcr_exportClones {
    cpus "$params.cpus"
    memory "$params.memory"
    time "$params.time"

    tag "$SampleID"
    label "mixcr"

    publishDir "$params.outdir/MiXCR/$SampleID",
        mode: 'copy',
        overwrite: true
    input:
    tuple val(SampleID),file(igmt), path(clns)

    output:
    path("*")

    script:
    """
    #Export clones:
    mixcr exportClones \
            --impute-germline-on-export \
            -allAAFeaturesImputed CDR1Begin FR4End \
            -allNFeaturesImputed CDR1Begin FR4End \
            -nFeatureImputed FR1 \
            -aaFeatureImputed FR1 \
            --force-overwrite \
            -nMutations VRegion \
            -nMutations JRegion \
            -nMutations DRegion \
        ${clns} \
        ${SampleID}.tsv

    #Export MiXCR Reports:
    #json-formatted:
    mixcr exportReports \
    --json \
    ${clns} \
    ${SampleID}_MiXCR.report.json

    #txt-formatted:
    mixcr exportReports \
    ${clns} \
    ${SampleID}_MiXCR.report.txt
    """
}

process platypus {
    cpus "$params.cpus"
    time "$params.time"
    memory "$params.memory"

    label 'platypus'
    tag "Platypus"

    publishDir "$params.outdir/Platypus",
        pattern: '*',
        mode: 'copy',
        overwrite: true
    input:
    tuple path(inputDir), path(sampleInfo), file(platypus), file(utils)

    output:
    path("*")

    script:
    """
    Rscript ${platypus} -s ${params.species} -m ${sampleInfo}
    """
}