# ENCODE ChIP-seq Methylation Data

This documents the dataset selection for analysis with Laura's SNV datasets.

## Open Chromatin Markers

Data were downloaded for both lung and liver.

**NB**: Not all data were available for mm10; so mm9 data were selected.
**NB**: As with the histone markers, the data quality were very low.

Liver: ENCSR000CNI
Lung: ENCSR000CNM

### Mouse Liver DNase-seq

#### ENCODE Filters:

| Filter               | Selection
|:---------------------|:---------
| Assay Type           | `DNA accessibility`
| Assay Title          | `DNase-seq`
| Status               | `released`
| Organism             | `Mus musculus`
| Biosample term name  | `liver`
| Organ                | `liver`
| Life stage           | `adult`
| Platform             | `Illumina Genome Analyzer I`

# [Metadata URL](https://www.encodeproject.org/metadata/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&biosample_ontology.term_name=liver&biosample_ontology.organ_slims=liver&replicates.library.biosample.life_stage=adult&assay_slims=DNA+accessibility&assay_title=DNase-seq&files.platform.term_name=Illumina+Genome+Analyzer+I)

#### Resulting Data Sets

Experiment ID                                                         | Line    | Target       | BED ID              | S3 URL              |
:---------------------------------------------------------------------|:--------|:-------------|:--------------------|---------------------|
[ENCSR000CNI](https://www.encodeproject.org/experiments/ENCSR000CNI/) | C57BL/6 | `DNase-seq`  | (NO PROCESSED DATA) | (NO PROCESSED DATA) |


### Mouse Lung DNase-seq

#### ENCODE Filters:

| Filter               | Selection
|:---------------------|:---------
| Assay Type           | `DNA accessibility`
| Assay Title          | `DNase-seq`
| Status               | `released`
| Organism             | `Mus musculus`
| Biosample term name  | `lung`
| Organ                | `lung`
| Life stage           | `adult`
| Platform             | `Illumina Genome Analyzer I`

# [Metadata URL](https://www.encodeproject.org/metadata/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&replicates.library.biosample.life_stage=adult&assay_slims=DNA+accessibility&assay_title=DNase-seq&biosample_ontology.term_name=lung&biosample_ontology.organ_slims=lung&files.platform.term_name=Illumina+Genome+Analyzer+I
)

#### Resulting Data Sets

| Experiment ID                                                         | Line    | Target       | Replicate | Hotspot broadPeak BED ID                                        | S3 URL                                                                                |
|:----------------------------------------------------------------------|:--------|:-------------|----------:|:----------------------------------------------------------------|---------------------------------------------------------------------------------------|
| [ENCSR000CNM](https://www.encodeproject.org/experiments/ENCSR000CNM/) | C57BL/6 | `DNase-seq`  |         1 | [ENCFF769WJE](https://www.encodeproject.org/files/ENCFF769WJE/) | s3://encode-public/2018/07/25/6c726abc-5a15-40b5-8989-c39aabbeb1d2/ENCFF769WJE.bed.gz |
| [ENCSR000CNM](https://www.encodeproject.org/experiments/ENCSR000CNM/) | C57BL/6 | `DNase-seq`  |         2 | [ENCFF036UVG](https://www.encodeproject.org/files/ENCFF036UVG/) | s3://encode-public/2018/07/25/b0217ab5-1d16-4e95-84c7-96ce98c2a966/ENCFF036UVG.bed.gz |
| [ENCSR000CNM](https://www.encodeproject.org/experiments/ENCSR000CNM/) | C57BL/6 | `DNase-seq`  |         3 | [ENCFF749VUS](https://www.encodeproject.org/files/ENCFF749VUS/) | s3://encode-public/2018/07/25/98a47fdb-e05d-42a8-baf8-3469c91237d9/ENCFF749VUS.bed.gz |
