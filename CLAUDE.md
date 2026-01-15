# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics pipeline for prioritizing disease predisposition genes using large public genetic data. It's part of the St. Jude KIDS24 BioHackathon project, implementing the CoCoRV (Case-Control Rare Variant) association testing pipeline in a cloud-based environment.

## Architecture

### Core Technology Stack
- **Orchestration**: Nextflow (DSL2) workflow engine
- **Containerization**: Docker with two specialized containers
- **Cloud Platform**: DNAnexus for deployment and execution
- **Languages**: Groovy (Nextflow), Bash, R (statistical analysis), Python (Hail)
- **Key Tools**: bcftools, bedtools, ANNOVAR, VEP, Hail, SeqArray

### Workflow Structure
The pipeline follows a chromosome-parallel processing pattern:
1. **Preprocessing**: VCF normalization and quality control
2. **Annotation**: ANNOVAR/VEP variant annotation (optional skip if pre-annotated)
3. **Format Conversion**: VCF to GDS conversion for R analysis
4. **Population Genetics**: gnomAD-based ancestry prediction using Hail
5. **Association Testing**: CoCoRV statistical analysis (64GB RAM required)
6. **Postprocessing**: Result merging, QQ plots, FDR correction

### Key Architectural Patterns
- **Chromosome-level parallelization**: Input VCFs must be split by chromosome (`*.chr{N}.vcf.gz`)
- **Optional workflow branches**: Skip expensive steps if intermediate files provided
- **Memory-aware scaling**: Dynamic memory allocation with retry strategies
- **Container-first design**: All processes run in Docker containers

## Development Commands

### Container Management
```bash
# Build Docker images
docker build -t stithi/dnanexus-cocorv-nextflow-python:v1 docker-files/dockerfile-for-python/
docker build -t stithi/dnanexus-cocorv-nextflow-r:v1 docker-files/dockerfile-for-r/

# Pull pre-built images
docker image pull stithi/dnanexus-cocorv-nextflow-python:v1
docker image pull stithi/dnanexus-cocorv-nextflow-r:v1
```

### Local Testing
```bash
# Run pipeline locally
nextflow run DNAnexus-implementation/cocorv-nextflow/main.nf -params-file params.json

# Run with specific parameters
nextflow run DNAnexus-implementation/cocorv-nextflow/main.nf \
  --caseVCF /path/to/vcf/dir \
  --caseSample /path/to/sample/list \
  --resourceFiles /path/to/resources \
  --gnomADVersion v4exome
```

### DNAnexus Deployment
```bash
# Deploy to DNAnexus (requires dx-toolkit)
dx build -a --nextflow cocorv-nextflow --destination CoCoRV_dev:/applets/nextflow-workflows/cocorv-nextflow_v2
```

## Configuration

### Primary Config Files
- `nextflow.config`: Main pipeline configuration with Docker settings, parallelization limits (queue size: 22), and parameter defaults
- `nextflow_schema.json`: JSON Schema for parameter validation and documentation
- `dx-commands.txt`: DNAnexus deployment commands

### Key Parameters
- **Required**: `caseVCF`, `caseSample`, `resourceFiles`
- **Genome builds**: Supports GRCh37/hg19 and GRCh38/hg38
- **gnomAD versions**: v2exome, v4exome, v4genome
- **Memory requirements**: CoCoRV process needs 64GB, GDS conversion needs 16GB

## Data Formats and Inputs

### Expected Input Structure
- **Case VCFs**: Split by chromosome (`case.chr1.vcf.gz`, `case.chr2.vcf.gz`, etc.)
- **Sample list**: Tab-delimited file with sample IDs and case/control status
- **Resource files**: Pre-staged reference data in container at `/opt/utilities/`
- **Coverage BED**: Optional genomic regions for analysis restriction

### Output Structure
- **Association results**: Per-chromosome statistics merged with headers
- **Variant groups**: Gene-level variant aggregations
- **QQ plots and FDR**: Statistical visualizations and multiple testing correction

## Memory and Performance

### Resource Requirements
- **CoCoRV process**: 64GB RAM minimum
- **GDS conversion**: 16GB RAM
- **Parallelization**: Up to 22 concurrent processes
- **Retry strategy**: Memory scales with attempt number (4GB × attempt)

### Performance Considerations
- Chromosome splitting enables parallel processing
- Optional workflow branches reduce redundant computation
- Container overhead requires adequate system resources

## Common Development Patterns

### Error Handling
```groovy
errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'terminate' }
maxRetries 5
memory { 4.GB * task.attempt }
```

### Optional Process Execution
The pipeline includes skip logic for expensive processes:
- Skip annotation if `--annotation` files provided
- Skip population prediction if `--ancestry` file provided

### Reference Genome Handling
Configuration automatically resolves paths based on `gnomADVersion`:
- v2exome → GRCh37/hg19 references
- v4exome/v4genome → GRCh38/hg38 references

## Testing and Validation

This is a scientific pipeline without traditional unit tests. Validation occurs through:
- Parameter schema validation (`nextflow_schema.json`)
- Process-level input/output checks
- Manual testing with real genomic datasets
- Retry mechanisms for transient failures

## Important Notes

- This is **not a traditional software application** but a scientific workflow
- Heavy reliance on external bioinformatics tools and reference datasets
- Container environments must include all required tools and reference files
- Pipeline expects specific file naming conventions for chromosome-split VCFs
- Memory requirements are substantial due to genomic data processing