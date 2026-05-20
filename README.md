## AMALGAM (v1.0.0)

AMALGAM (Mending Long-read Gene Annotations Method) is a pipeline designed to perform transcript isoform discovery and quantification from long-read RNA-seq data generated on **human** samples, regardless of how fragmented the data may be.

**Note:** AMALGAM is still under active development. Please reach out to Robert Wang (wangr5@email.chop.edu) if you have any feedback or run into any issues!

### Dependencies:

* [GffCompare](https://github.com/gpertea/gffcompare) (v0.12.6)
  + v0.12.10 should not be used due to a formatting change to the `.tracking` file
* [HTSlib](https://www.htslib.org/) (v1.18 or above)
  * Only need to install `bgzip` and `tabix` (make these available in `$PATH`)
* [Python](https://www.python.org/downloads) (v3.8.18 or above)
  * [biopython](https://biopython.org/) (v1.78 or above)
  * [cvxpy](https://www.cvxpy.org/) (v1.5.2 or above)
  * [networkx](https://networkx.org/en/) (v2.6.3 or above)
  * [numpy](https://numpy.org/) (v1.24.4 or above)
  * [pandas](https://pandas.pydata.org/) (v2.0.3 or above)
  * [pysam](https://pysam.readthedocs.io/en/latest/api.html) (v0.23.0 or above)
  * [pytabix](https://pypi.org/project/pytabix/) (v0.1 or above)
* [StringTie](https://github.com/gpertea/stringtie) (v3.0.0 or above)

A conda environment can be created with the necessary dependencies:
```
conda create --prefix ./conda_env
conda activate ./conda_env
conda install -c conda-forge -c bioconda --file ./conda_requirements.txt
```

### Usage:

AMALGAM assumes that you have already aligned long-read RNA-seq reads to the GRCh38/hg38 human reference genome (indexed BAM file) using your favorite aligner of choice. AMALGAM also extensively works with gene annotations from [GENCODE](https://www.gencodegenes.org/).

**Step 1:** Perform de-novo transcriptome assembly from long-read RNA-seq alignments using StringTie (short-read mode, default settings). If you have multiple samples, run StringTie individually on each sample.

```
[/path/to/stringtie] -o [/path/to/stringtie/output/GTF] -p [number of threads] [/path/to/input/BAM/file]
```

The `stringtie` command with one thread is expected to use about 0.5GB of memory and take about 15 minutes (based on a 47GB BAM file with 55 million reads).

**Step 2:** Use GffCompare (v0.12.6) to combine the output GTF file(s) from StringTie with existing GENCODE comprehensive annotations

```
# Construct an input list of GTF files for GffCompare
# NOTE: Make sure that the first file in this list is the GENCODE comprehensive annotations file!
echo "/path/to/gencode/GTF" > "/path/to/input/GTF/list"
echo "/path/to/stringtie/output/GTF" >> "/path/to/input/GTF/list"

# If you have multiple samples, append paths to StringTie output GTF files to the end of your input list
#   echo "/path/to/stringtie/output/GTF" >> "/path/to/input/GTF/list"

[/path/to/gffcompare] -i [/path/to/input/GTF/list] -T -o [/path/to/gffcompare/output/prefix]
```

The `gffcompare` command is expected to use about 2.5GB of memory and take about 5 minutes.

**Step 3:** Use `Build_Transcriptome.py` to identify high-confidence, full-length transcripts from the merged GTF file produced by GffCompare. This script can be found in the `scripts` folder and requires the following two files, which can be found in the `assets` folder:
* `human.refTSS_v4.1.hg38.bed.gz`: A tabix-indexed and bgzip-compressed BED file containing genomic coordinates of transcription start sites in the human reference genome (GRCh38/hg38) based on annotations from the [RefTSS database](https://reftss.riken.jp/) (v4.1)
* `atlas.clusters.2.0.GRCh38.bed.gz`: A tabix-indexed and bgzip-compressed BED file containing genomic coordinates of polyadenylation sites in the human reference genome (GRCh38/hg38) based on annotations from the [PolyASite database](https://polyasite.unibas.ch/) (v2.0)

```
python [/path/to/Build_Transcriptome.py] -i [/path/to/gffcompare/output/prefix] \
    -g [/path/to/gencode/GTF] \
    -f [/path/to/reference/genome/FASTA] \
    -x [/path/to/human.refTSS_v4.1.hg38.bed.gz] \
    -y [/path/to/atlas.clusters.2.0.GRCh38.bed.gz] \
    -o [/path/to/output.transcripts.gtf]
```

The `Build_Transcriptome.py` command is expected to use about 4.5GB of memory and take about 20 minutes.

**Step 4:** Use `Annotate_ORF.py` to annotate open-reading frames for each transcript in the GTF file produced by `Build_Transcriptome.py`. This script can also be found in the `scripts` folder.

```
python [/path/to/Annotate_ORF.py] -i [/path/to/output.transcripts.gtf] \
    -a [/path/to/gencode/GTF] \
    -f [/path/to/reference/genome/FASTA] \
    -o [/path/to/output.transcripts.updated.gtf]

# Compress and index the updated GTF file produced by Annotate_ORF.py
sort -k1,1V -k4,4g -k5,5g [/path/to/output.transcripts.updated.gtf] | bgzip > [/path/to/output.transcripts.updated.gtf.gz]
tabix -p gff [/path/to/output.transcripts.updated.gtf.gz]
```

The `Annotate_ORF.py` command is expected to use about 8GB of memory and take about 5 hours.

**Note:** If you want to skip this step and go straight to transcript quantification, you can compress/index the GTF file produced by `Build_Transcriptome.py` as follows:

```
sort -k1,1V -k4,4g -k5,5g [/path/to/output.transcripts.gtf] | bgzip > [/path/to/output.transcripts.updated.gtf.gz]
tabix -p gff [/path/to/output.transcripts.updated.gtf.gz]
```

**Step 5:** Use `Quantify_Transcripts.py` to quantify transcripts from `/path/to/output.transcripts.updated.gtf.gz` in an individual long-read RNA-seq BAM file. This script can also be found in the `scripts` folder.

```
python [/path/to/Quantify_Transcripts.py] -i [/path/to/input/BAM/file] \
    -g [/path/to/output.transcripts.updated.gtf.gz] \
    -o [/path/to/transcript/counts/file]
```

The `Quantify_Transcripts.py` command is expected to use about 3.5GB of memory and take about 16 hours.
