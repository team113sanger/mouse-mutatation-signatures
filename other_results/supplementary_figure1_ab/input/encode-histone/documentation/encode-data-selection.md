# ENCODE ChIP-seq Methylation Data

This documents the dataset selection for analysis with Laura's SNV datasets.

## Specific Histone Markers

Data were downloaded for both lung and liver.

**NB**: Very different data qualities and sources were present for the diffferent organs, so it will be difficult to compare between them.

### Mouse Liver Histone ChIP-seq

#### ENCODE Filters:

| Filter               | Selection
|:---------------------|:---------
| Status               | `released`
| Genome assembly      | `mm10`
| Target category      | `histone`
| Organism             | `Mus musculus`
| Biosample term name  | `liver`
| Organ                | `liver`
| Life stage           | `adult`
| Available file types | `bed narrowPeak`
| Platform             | `Illumina HiSeq 2000`

[Metadata URL](https://www.encodeproject.org/metadata/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&biosample_ontology.term_name=liver&target.investigated_as=histone&replicates.library.biosample.life_stage=adult&files.platform.term_name=Illumina+HiSeq+2000&assembly=mm10)

#### Resulting Data Sets

Experiment ID                                                         | Line    | Target     | Replicated Peak BED ID                                          | S3 URL                                                                                  |
:---------------------------------------------------------------------|:--------|:-----------|:----------------------------------------------------------------|-----------------------------------------------------------------------------------------|
[ENCSR000CDH](https://www.encodeproject.org/experiments/ENCSR000CDH/) | C57BL/6 | `H3K27ac`  | [ENCFF545NYE](https://www.encodeproject.org/files/ENCFF545NYE/) | `s3://encode-public/2016/05/03/7b6b3773-64be-41df-963b-ff9d6c165a41/ENCFF545NYE.bed.gz` |
[ENCSR000CEN](https://www.encodeproject.org/experiments/ENCSR000CEN/) | C57BL/6 | `H3K27me3` | [ENCFF506DLC](https://www.encodeproject.org/files/ENCFF506DLC/) | `s3://encode-public/2016/11/15/2ca228ed-600b-4b40-bddd-9b0ee7f24176/ENCFF506DLC.bed.gz` |
[ENCSR000CEO](https://www.encodeproject.org/experiments/ENCSR000CEO/) | C57BL/6 | `H3K36me3` | [ENCFF718JST](https://www.encodeproject.org/files/ENCFF718JST/) | `s3://encode-public/2016/11/15/0aae6c4b-656f-459d-97eb-3efc83895dbb/ENCFF718JST.bed.gz` |
[ENCSR000CEP](https://www.encodeproject.org/experiments/ENCSR000CEP/) | C57BL/6 | `H3K79me2` | [ENCFF306QLL](https://www.encodeproject.org/files/ENCFF306QLL/) | `s3://encode-public/2016/11/15/33ca85ae-cd83-403f-a129-458563431186/ENCFF306QLL.bed.gz` |
[ENCSR000CEQ](https://www.encodeproject.org/experiments/ENCSR000CEQ/) | C57BL/6 | `H3K9ac`   | [ENCFF787EAZ](https://www.encodeproject.org/files/ENCFF787EAZ/) | `s3://encode-public/2016/11/15/f401d83a-14c1-487a-a704-49cebc3a056b/ENCFF787EAZ.bed.gz` |

### Mouse Lung Histone ChIP-seq

#### ENCODE Filters:

| Filter               | Selection
|:---------------------|:---------
| Status               | `released`
| Genome assembly      | `mm10`
| Target category      | `histone`
| Organism             | `Mus musculus`
| Biosample term name  | `lung`
| Organ                | `lung`
| Life stage           | `adult`
| Available file types | `bed narrowPeak`
| Platform             | `Illumina Genome Analyzer II`

**NB**: No Illumina HiSeq 2000 data available; so we use Illumina Genome Analyser II.

**NB**: These samples have low read count warnings.

[Metadata URL](https://www.encodeproject.org/metadata/?type=Experiment&biosample_ontology.organ_slims=lung&biosample_ontology.term_name=lung&replicates.library.biosample.life_stage=adult&assembly=mm10&target.investigated_as=histone&status=released&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&files.file_type=bed+narrowPeak&files.platform.term_name=Illumina+Genome+Analyzer+II)

#### Resulting Data Sets

Experiment ID                                                         | Line    | Target    | Replicated Peak BED ID                                          | S3 URL                                                                                  |
:---------------------------------------------------------------------|:--------|:----------|:----------------------------------------------------------------|-----------------------------------------------------------------------------------------|
[ENCSR000CAR](https://www.encodeproject.org/experiments/ENCSR000CAR/) | C57BL/6 | `H3K4me3` | [ENCFF508WEP](https://www.encodeproject.org/files/ENCFF508WEP/) | `s3://encode-public/2016/05/03/bbb83850-f8f8-4dd8-a16a-00cf08444398/ENCFF508WEP.bed.gz` |
[ENCSR000CAQ](https://www.encodeproject.org/experiments/ENCSR000CAQ/) | C57BL/6 | `H3K4me1` | [ENCFF865RFC](https://www.encodeproject.org/files/ENCFF865RFC/) | `s3://encode-public/2016/11/15/5dee72be-740b-4a98-8877-3efbfbfcc6f9/ENCFF865RFC.bed.gz` |
