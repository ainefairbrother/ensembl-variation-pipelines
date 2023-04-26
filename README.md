# ensembl-variation-pipelines

## Usage

Example command:

```
module add nextflow-22.10.1-gcc-11.2.0-ju5saqw

nextflow \
-C ${ENSEMBL_ROOT}/ensembl-variation-pipelines/nextflow/nf_config/nextflow.config \
  run ${ENSEMBL_ROOT}/ensembl-variation-pipelines/nextflow/workflows/vcf_prepper.nf \
  -profile lsf \
  --input_config <path/to/input_config.json> \
  --output_dir <path/to/output_dir> \
  --bin_size 250000 \
  --remove_patch 1 \
  --skip_vep 0 \
  --skip_create_config 0 \
  --ini_file <path/to/ini_file> \
  --rank_file <path/to/variation_consequnce_rank.json> \
  --version 108 \
  -resume
```

## Options

- `input_config` : (optional) Give the full path of the input configuration file, default: `ensembl-variation-pipelines/nextflow/nf_config/input_sources.json`.

This is a json configuration file containing the source file information. You can specify the location of VCF files and also their priority in the config. The priority is used for generating focus track files. Variants from VCF file that have higher priority value will be stacked on top of variants from VCF file with lower priority value. e.g. - 

```
{
  "homo_sapiens_grch38" : [
    {
        "source_name": "dbSNP",
        "file_location": "/path/to/GCF_000001405.39_VEP.vcf.gz",
        "priority": 1
    },
    {
        "source_name": "GWAS",
        "file_location": "/path/to/input_gwas_catalog_v1.0.2-associations_e107_r2022-09-14_test.vcf.gz"
    }
  ]
}
```

In the above case, GWAS variants will be stacked on top of dbSNP variant. What it means is that, if GWAS and dbSNP both have a variant with name `rs100` only the variant from dbSNP will be kept. Priority is optional and by default it's value is 100. 

- `output_dir` : (optional) Give the full path of the dir where the outputs will be generated, default: `/nfs/production/flicek/ensembl/variation/new_website`

The generated VEP VCF and track files will be stored there. The directory structure would be - 

```
<output_dir>
  |
  -- <genome 1>
      |
      -- <source 1>
          |
          -- vcfs
              |
              -- *.vcf.gz
          -- tracks
              |
              -- *.bb
              -- *.bw
      -- <source 2>
      .
      .
      -- focus
          |
          -- vcfs
              |
              -- *.vcf.gz
          -- tracks
              |
              -- *.bb
              -- *.bw
  -- <genome 2>
  .
  .
```

- `bin_size` : (optional) This is a nextflow-vep pipeline parameter and determines the number of variants to split the VCF file by, default: `250000`.

- `remove_patch` : (optional) If value is 1, the pipeline will ignore sequence region with `PATCH`, `TEST`, `CTG` in name, default: `1`.

- `skip_vep` : (optional) If value is 1, the pipeline will skip running `nextflow-vep`, default: `0`. 

If `nextflow-vep` is skipped, the VCF file provided in the `input_config.json` should be VCF files that have already been run through VEP (with appropriate configuration required by vcf_prepper pipeline)

- `skip_create_config` : (optional) If value is 1, the pipeline will skip creating the config files, default: `0`.

The config files that falls under this condition are -

`synonyms files`: File containing sequence region synonyms and their original names in tab limited file created from Ensembl core database. Needed for `renameChr` step.

`chrom sizes files`: File containing sequence regions and their lengths in tab limited file. Used by UCSC tools to create bigWig and bogBeds.

- `ini_file` : (optional) A INI file that is used by the `createConfigs` step to generate the config files, default: `ensembl-variation-pipelines/nextflow/nf_config/DEFAULT.ini`.

The file format is -
```
[database]
host = <hostname>
port = <port>
user = <user>

[grch37_database]
host = <hostname>
port = <port>
user = <user>
```

Add the appropriate server information so that the core database related to the species and Ensembl version can be found. If `skip_create_config` is set to `1` this option is ignored.

- `rank_file` : (optional) Give the full path of the rank file to be generated and used, default: `ensembl-variation-pipelines/nextflow/nf_config/variation_consequnce_rank.json`

This rank file contains the rank of variant consequence and used to determine the most severe consequence of a variant. The pipeline generate this file automatically using Ensembl Variation api and later use it in the `vcfToBed` step.

Make sure you are checked out to appropriate `ensembl-variation` repository to get the appropriate ranks. ideally it would be the same version as provided by the `--version` parameter.

- `version` : (optional) Give the Ensembl version to be used, default: `108`.

Currently, used to get the appropriate Ensembl core database when creating the configs in `createConfigs` step.

## Pipeline Flow

```mermaid
flowchart TD
    p0((Channel.fromList))
    p1((Channel.fromList))
    p2((Channel.fromList))
    p3((Channel.fromList))
    p4((Channel.fromList))
    p5(( ))
    p6(( ))
    p7[createRankFile]
    p8[createConfigs]
    p9([count])
    p10([count])
    p11([combine])
    p12([subscribe])
    p13[vep:processInput]
    p14[vep:checkVCF]
    p15(( ))
    p16[vep:generateSplits]
    p17([transpose])
    p18[vep:splitVCF]
    p19([transpose])
    p20[vep:runVEP]
    p21(( ))
    p22([groupTuple])
    p23[vep:mergeVCF]
    p24([map])
    p25([map])
    p26(( ))
    p27[renameChr]
    p28[removeDupIDs]
    p29[indexVCF]
    p30[readChrVCF]
    p31([transpose])
    p32[splitChrVCF]
    p33([transpose])
    p34[vcfToBed]
    p35([groupTuple])
    p36[concatBed]
    p37[bedToBigBed]
    p38[bedToBigWig]
    p39([groupTuple])
    p40[createFocusTrack]
    p0 -->|vcf| p9
    p1 -->|output_dir| p13
    p2 -->|vep_config| p10
    p3 -->|genomes| p8
    p4 -->|sources| p5
    p6 -->|rank_file| p7
    p7 --> p8
    p8 --> p34
    p9 --> p11
    p10 --> p11
    p11 --> p12
    p0 -->|vcf| p13
    p2 -->|vep_config| p13
    p13 --> p14
    p14 --> p16
    p15 -->|bin_size| p16
    p16 --> p17
    p17 --> p18
    p18 --> p19
    p19 --> p20
    p20 --> p22
    p20 --> p21
    p22 --> p23
    p23 --> p24
    p24 --> p27
    p23 --> p25
    p25 --> p27
    p23 -->|input_file| p27
    p26 -->|priorities| p27
    p27 --> p28
    p28 --> p29
    p29 --> p30
    p30 --> p31
    p31 --> p32
    p32 --> p33
    p33 --> p34
    p34 --> p35
    p35 --> p36
    p36 --> p37
    p36 --> p38
    p36 --> p39
    p39 --> p40
```
