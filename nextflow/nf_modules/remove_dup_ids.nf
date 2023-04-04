#!/usr/bin/env nextflow

/*
* This script rename all the duplicated rsID to '.' in a VCF file 
*/

process removeDupIDs {
  input:
  tuple val(original), path(vcfFile)
  
  output:
  tuple val(original), path("processed-${original}-*.vcf.gz"), path("processed-${original}-*.vcf.gz.tbi"), emit: file
  
  shell:
  '''
  # format output file name
  output_file=processed-!{original}-!{vcfFile}
  
  # get the duplicated rsID list and load it into an array
  bcftools view --no-version !{vcfFile} | awk '!/^#/{if (uniq[$3]++ && $3 != ".") print($3)}' | sort -u > duplicated_ids.txt 
  export duplicated_id=`cat duplicated_ids.txt | xargs`
  
  # only keep the regions that we want to keep
  awk '
    BEGIN {
      split(ENVIRON["duplicated_id"], t_ids, " ")
      for (i in t_ids) ids[t_ids[i]] = ""
    }
    /^##/ {
      print
    }
    /^#/ {
      OFS="\t"
      print
    }
    !/^#/ {
      OFS="\t"
      if ($3 in ids) {
        $3 = ".";
      }
      print;
    }
  ' !{vcfFile} > ${output_file}
  
  # bgzip and index for next step
  bgzip ${output_file}
  bcftools index -t ${output_file}.gz
  '''
}