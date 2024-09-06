#!/usr/bin/env nextflow

process coverageIntersect {
  tag "${caseBed}"
  publishDir "${params.outdir}", mode: 'copy'
  container 'stithi/dnanexus-cocorv-nextflow-python:v1'

  memory { 4.GB * task.attempt }

  input:
    path caseBed
    path controlBed

  output:
    path "intersect.coverage10x.bed.gz"

  script:
  """
  bedtools intersect -sorted -a ${caseBed} -b ${controlBed} > \
        "intersect.coverage10x.bed.gz"
  """
}

process normalizeQC {
  //tag "${caseVCFPrefix}_${chr}"
  publishDir "${params.outdir}/vcf_vqsr_normalizedQC", mode: 'copy'
  container 'stithi/dnanexus-cocorv-nextflow-python:v1'

  memory { 4.GB * task.attempt }

  input:
    path(case_vcf_file)
    path(ref_fasta)

  output:
    path("*.biallelic.leftnorm.ABCheck.vcf.gz")
    path("*.biallelic.leftnorm.ABCheck.vcf.gz.tbi")

  script:
  """
  bash /opt/utilities/vcfQCAndNormalize.sh ${case_vcf_file} ${ref_fasta}
  """
}

process annotate {
  tag "${chr}"
  publishDir "${params.outdir}/annotation", mode: 'copy'
  container 'stithi/dnanexus-cocorv-nextflow-python:v1'

  memory { 16.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 5

  input:
    path(normalizedQCedVCFFile)
    path(indexFile)
    val reference
    path(annovarFolder)

  output:
    tuple val("${chr}"),
          path("${chr}.annotated.vcf.gz")

  script:
  chr = normalizedQCedVCFFile.simpleName
  refbuild = null
  if (reference == "GRCh37") {
    refbuild="hg19"  
  }
  else {
    if (reference == "GRCh38") {
     refbuild="hg38"   
    }
  }
  """
  if [[ ${params.addVEP} != "T" ]]; then
    outputPrefix="${chr}.annotated"
  else 
    outputPrefix="${chr}.annotated.annovar"
  fi
  bash /opt/utilities/annotate.sh ${normalizedQCedVCFFile} ${annovarFolder} ${refbuild} \${outputPrefix} ${params.VCFAnno} ${params.toml} ${params.protocol} ${params.operation}

  if [[ ${params.addVEP} == "T" ]]; then
    module load perl/5.28.1
    module load htslib/1.10.2
    module load bcftools/1.15.1
    module load samtools/1.10
    bash /opt/utilities/annotateVEPWithOptions.sh ${chr}.annotated.annovar.vcf.gz ${reference} ${chr}.annotated ${params.reference} ${params.vepFolder} ${params.cache} ${params.lofteeFolder} ${params.lofteeDataFolder} ${params.caddSNV} ${params.caddIndel} ${params.spliceAISNV} ${params.spliceAIIndel} ${params.perlThread} ${params.AM} ${params.REVEL} 1 ${params.VEPAnnotations}
  fi
  """
}

process skipAnnotation {
  tag "${chr}"
  container 'stithi/dnanexus-cocorv-nextflow-python:v1'

  memory { 4.GB * task.attempt }

  input:
    tuple val(chr), path(normalizedQCedVCFFile), path(indexFile)

  output:
    tuple val("${chr}"),
          path("${chr}.annotated.vcf.gz")

  script:
  if (chr == "NA") {
    annotated=params.caseAnnotatedVCFPrefix+params.caseAnnotatedVCFSuffix
  } else {
    annotated=params.caseAnnotatedVCFPrefix+chr+params.caseAnnotatedVCFSuffix
  }
  """
  ln -s ${annotated} ${chr}.annotated.vcf.gz
  """
}

process caseGenotypeGDS {
  tag "${chr}"
  publishDir "${params.outdir}/vcf_vqsr_normalizedQC", mode: 'copy'
  container 'stithi/dnanexus-cocorv-nextflow-r:v1'

  cpus params.cpus
  memory { 16.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 5

  input: 
    path(normalizedQCedVCFFile)
    path(indexFile)

  output: 
    tuple val("${chr}"),
          path("${chr}.biallelic.leftnorm.ABCheck.vcf.gz.gds")

  script:
  chr = normalizedQCedVCFFile.simpleName
  """
  Rscript /opt/utilities/vcf2gds.R ${normalizedQCedVCFFile} ${chr}.biallelic.leftnorm.ABCheck.vcf.gz.gds ${params.cpus}
  """
}

process caseAnnotationGDS {
  tag "${chr}"
  publishDir "${params.outdir}/annotation", mode: 'copy'
  container 'stithi/dnanexus-cocorv-nextflow-r:v1'

  memory { 16.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 1

  input: 
    tuple val(chr), path(annotatedFile)

  output: 
    tuple val("${chr}"),
          path("${chr}.annotated.vcf.gz.gds")

  script:
  """
  Rscript /opt/utilities/vcf2gds.R ${annotatedFile} ${chr}.annotated.vcf.gz.gds 1
  """
}

process extractGnomADPositions {
  tag "${chr}"
  publishDir "${params.outdir}/gnomADPosition", mode: 'copy'
  container 'stithi/dnanexus-cocorv-nextflow-python:v1'

  memory { 4.GB * task.attempt }

  input: 
    path(normalizedQCedVCFFile)
    path(indexFile)
    path(gnomADPCPosition)

  output: 
    path "${chr}.extracted.vcf.gz"

  script:
  chr = normalizedQCedVCFFile.simpleName
  """
  bcftools view -R ${gnomADPCPosition} -Oz -o ${chr}.extracted.vcf.gz ${normalizedQCedVCFFile}
  """
}

process mergeExtractedPositions {
  publishDir "${params.outdir}/gnomADPosition", mode: 'copy'
  container 'stithi/dnanexus-cocorv-nextflow-python:v1'

  memory { 4.GB * task.attempt }

  input: 
    path(extractedVCFFile)

  output: 
    path("all.extracted.vcf.bgz")

  script:
  """
  bcftools concat -Oz -o "all.extracted.vcf.bgz" ${extractedVCFFile}
  """
}

process RFPrediction {
  publishDir "${params.outdir}/gnomADPosition", mode: 'copy'
  container 'stithi/dnanexus-cocorv-nextflow-python:v1'

  memory { 4.GB * task.attempt }

  input: 
    path(VCFForPrediction)
    val reference

  output: 
    path "PC.population.output.gz"
    path "casePopulation.txt"

  script:
  """
  if [[ ${reference} == "GRCh37" ]]; then
    bash /opt/utilities/gnomADPCAndAncestry.sh /opt/utilities/hail_hg19/gnomad.r2.1.pca_loadings.ht /opt/utilities/hail_hg19/gnomad.r2.1.RF_fit.onnx ${VCFForPrediction} ${reference} "PC.population.output.gz" 0.9 "casePopulation.txt"
  else
    bash /opt/utilities/gnomADPCAndAncestry.sh /opt/utilities/hail_hg38/gnomad.v4.0.pca_loadings.ht /opt/utilities/hail_hg38/gnomad.v4.0.RF_fit.onnx ${VCFForPrediction} ${reference} "PC.population.output.gz" 0.75 "casePopulation.txt"
  fi
  """
}

process CoCoRV {
  tag "${chr}"
  container 'stithi/dnanexus-cocorv-nextflow-r:v1'

  memory { 64.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 1

  input: 
    tuple val(chr), path(caseGenotypeGDS), path(caseAnnoGDS) 
    path intersectBed
    path ancestryFile
    path controlGDS
    path caseSample
    val ACANConfig
    val variantExclude
    val highLDVariantFile
    val reference


  output: 
    path("${chr}.association.tsv") 
    path("${chr}.case.group")
    path("${chr}.control.group")

  script:
  if (params.gnomADVersion == "v4exome") {
    controlAnnoGDS = controlGDS + "/chr" + chr + ".annovar.vep.vcf.gz.gds"
    controlCountGDS = controlGDS + "/gnomad.exomes.v4.1.sites.chr" + chr + ".vcf.bgz.gds"
  } 
  else if (params.gnomADVersion == "v4genome") {
    controlAnnoGDS = controlGDS + "/chr" + chr + ".annovar.vep.vcf.gz.gds"
    controlCountGDS = controlGDS + "/gnomad.genomes.v4.1.sites.chr" + chr + ".vcf.bgz.gds"
  }
  else if (params.gnomADVersion == "v2exome") {
    controlAnnoGDS = controlGDS + "/chr" + chr + ".vep.fixed.vcf.gz.gds"
    controlCountGDS = controlGDS + "/gnomad.exomes.r2.1.1.sites." + chr + ".vcf.bgz.gds"
  }
  """
  otherOptions=""
  if [[ "${params.CoCoRVOptions}" != "NA" ]]; then
    otherOptions="${params.CoCoRVOptions}"
  fi
  if [[ "${params.groupFunctionConfig}" != "NA" ]]; then
    otherOptions="\${otherOptions} --variantGroupCustom ${params.groupFunctionConfig}"
  fi
  if [[ "${params.extraParamJason}" != "NA" ]]; then
    otherOptions="\${otherOptions} --extraParamJason ${params.extraParamJason}"
  fi
  if [[ "${params.annotationUsed}" != "NA" ]]; then
    otherOptions="\${otherOptions} --annotationUsed ${params.annotationUsed}"
  fi
  if [[ "${highLDVariantFile}" != "NA" ]]; then
    otherOptions="\${otherOptions} --highLDVariantFile ${highLDVariantFile}"
  fi
  if [[ "${variantExclude}" != "NA" ]]; then
    otherOptions="\${otherOptions} --variantExclude ${variantExclude}"
  fi
  if [[ "${params.variantInclude}" != "NA" ]]; then
    otherOptions="\${otherOptions} --variantInclude ${params.variantInclude}"
  fi
 
  Rscript /opt/utilities/CoCoRV_wrapper.R \
      --sampleList ${caseSample} \
      --outputPrefix ${chr} \
      --AFMax ${params.AFMax} \
      --bed ${intersectBed} \
      --variantMissing ${params.variantMissing} \
      --variantGroup ${params.variantGroup} \
      --removeStar \
      --ACANConfig ${ACANConfig} \
      --caseGroup ${ancestryFile} \
      --minREVEL ${params.REVELThreshold} \
      --checkHighLDInControl \
      --pLDControl ${params.pLDControl} \
      --fullCaseGenotype \
      --reference ${reference} \
      --gnomADVersion ${params.gnomADVersion} \
      --controlAnnoGDSFile ${controlAnnoGDS} \
      --caseAnnoGDSFile ${caseAnnoGDS} \
      --batchSize ${params.batchSize} \
      \${otherOptions} \
      ${controlCountGDS} \
      ${caseGenotypeGDS}
  """
}

process mergeCoCoRVResults {
  publishDir "${params.outdir}/CoCoRV", mode: 'copy'
  container 'stithi/dnanexus-cocorv-nextflow-r:v1'

  memory { 4.GB * task.attempt }

  input:
    path associationResult
    path caseVariants
    path controlVariants

  output:
    tuple path("association.tsv"), 
     path("kept.variants.case.txt"), path("kept.variants.control.txt")

  script:
  """
  body() {
    IFS= read -r header
    printf '%s\n' "\$header"
    "\$@"
  }

  i=0
  cat ${caseVariants} > "kept.variants.case.txt"
  cat ${controlVariants} > "kept.variants.control.txt"
  for file in ${associationResult}; do
    echo \${file}
    if [[ \${i} == 0 ]]; then
      cat \${file} > "association.tsv.tmp"
      i=1
    else
      tail -n+2 \${file} >> "association.tsv.tmp"
    fi
  done
  cat "association.tsv.tmp" | sort -gk2,2 > "association.tsv"
  """
}

process QQPlotAndFDR {
  publishDir "${params.outdir}/CoCoRV", mode: 'copy'
  container 'stithi/dnanexus-cocorv-nextflow-r:v1'

  memory { 16.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
  maxRetries 1  

  input: 
    tuple path("association.tsv"), 
     path("kept.variants.case.txt"), path("kept.variants.control.txt")

  output:
    path "association.tsv.dominant.nRep1000*"

  script:
  """ 
  Rscript /opt/utilities/QQPlotAndFDR.R "association.tsv" \
           "association.tsv.dominant.nRep1000" --setID gene \
      --outColumns gene --n 1000 \
      --pattern "case.*Mutation.*_DOM\$|control.*Mutation.*_DOM\$" \
      --FDR
  """
}

// main pipeline
workflow {
  if (params.gnomADVersion == "v2exome") {
    reference = "GRCh37"
    refFasta = params.resourceFiles + "/GRCh37-lite.fa"
    annovarFolder = params.resourceFiles + "/annovarFolder"
    controlBed = params.resourceFiles + "/gnomadv2_exome/coverage10x.bed.gz"
    controlGDS = params.resourceFiles + "/gnomadv2_exome"
    gnomADPCPosition = params.resourceFiles + "/gnomadv2_exome/hail_positions.GRCh37.chr.pos.tsv"
    ACANConfig = params.Build_hg19.ACANConfig
    variantExclude = params.Build_hg19.variantExclude
    highLDVariantFile = params.Build_hg19.highLDVariantFile
  }
  else if (params.gnomADVersion == "v4exome") {
    reference = "GRCh38"
    refFasta = params.resourceFiles + "/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
    annovarFolder = params.resourceFiles + "/annovarFolder"
    controlBed = params.resourceFiles + "/gnomADv4.1_exome/coverage10x.bed.gz"
    controlGDS = params.resourceFiles + "/gnomADv4.1_exome"
    gnomADPCPosition = params.resourceFiles + "/gnomADv4.1_genome/hail_positions.GRCh38.v4.chr.pos.tsv"
    ACANConfig = params.Build_hg38.ACANConfig
    variantExclude = params.Build_hg38.variantExclude
    highLDVariantFile = params.Build_hg38.highLDVariantFile
  }
  else if (params.gnomADVersion == "v4genome") {
    reference = "GRCh38"
    refFasta = params.resourceFiles + "/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
    annovarFolder = params.resourceFiles + "/annovarFolder"
    controlBed = params.resourceFiles + "/gnomADv4.1_genome/coverage10x.bed.gz"
    controlGDS = params.resourceFiles + "/gnomADv4.1_genome"
    gnomADPCPosition = params.resourceFiles + "/gnomADv4.1_genome/hail_positions.GRCh38.v4.chr.pos.tsv"
    ACANConfig = params.Build_hg38.ACANConfig
    variantExclude = "NA"
    highLDVariantFile = params.Build_hg38.highLDVariantFile
  }
  
  // coverage  
  if (params.caseBed != null) {
    coverageIntersect(params.caseBed, controlBed)
    intersectChannel = coverageIntersect.out
  }
  else {
    intersectChannel = Channel.value(controlBed)
  }

  // normalize and QC
  case_vcf_ch = Channel
        .fromPath(params.caseVCF + "/*.vcf.gz")
  normalizeQC(case_vcf_ch, refFasta)

  // extract VQSR failed variants
  //extractVQSRFailedVariants(normalizeQC.out[0], normalizeQC.out[1])

  // merge VQSR failed variants with gnomad VQSR failed variants
  //mergeVQSRFailedVariants(extractVQSRFailedVariants.out.collect())
  
  // annotate
  if (params.caseAnnotatedVCFPrefix == "NA") {
    annotate(normalizeQC.out[0], normalizeQC.out[1], reference, annovarFolder)
    annotateChannel = annotate.out
  } else {
    skipAnnotation(normalizeQC.out)
    annotateChannel = skipAnnotation.out
  }

  // case genoypte vcf to gds
  caseGenotypeGDS(normalizeQC.out[0], normalizeQC.out[1])

  // case annotation to gds
  caseAnnotationGDS(annotateChannel)
  
  // run gnomAD based population prediction
  if (params.casePopulation == null) {
    // extract gnomAD positions
    extractGnomADPositions(normalizeQC.out[0], normalizeQC.out[1], gnomADPCPosition)

    // merge extracted gnomAD positions
    mergeExtractedPositions(extractGnomADPositions.out.collect())

    // run gnomAD based population prediction
    RFPrediction(mergeExtractedPositions.out, reference)
    populationChannel = RFPrediction.out[1]
  } else {
    populationChannel = Channel.fromPath(params.casePopulation)
  }

  // run CoCoRV
  // RFPrediction.out.view()
  CoCoRV(caseGenotypeGDS.out.join(caseAnnotationGDS.out), 
    intersectChannel,
    populationChannel,
    controlGDS,
    params.caseSample,
    ACANConfig,
    variantExclude,
    highLDVariantFile,
    reference)

  // merge CoCoRV results
  mergeCoCoRVResults(CoCoRV.out[0].collect(), CoCoRV.out[1].collect(), 
    CoCoRV.out[2].collect())

  // QQ plot and FDR
  QQPlotAndFDR(mergeCoCoRVResults.out)
}

