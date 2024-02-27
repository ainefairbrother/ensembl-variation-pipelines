#!/usr/bin/env python3

import sys
from cyvcf2 import VCF, Writer
from Bio import bgzf
import argparse

HEADERS = [
    {'ID': 'NTCSQ', 'Description': 'Number of regulatory consequences', 'Type':'Integer', 'Number': 'A'},
    {'ID': 'NRCSQ', 'Description': 'Number of transcripts consequences', 'Type':'Integer', 'Number': 'A'},
    {'ID': 'NGENE', 'Description': 'Number of overlapped gene', 'Type':'Integer', 'Number': 'A'},
    {'ID': 'NVPHN', 'Description': 'Number of associated variant-linked phenotypes', 'Type':'Integer', 'Number': 'A'},
    {'ID': 'NGPHN', 'Description': 'Number of associated gene-linked phenotypes', 'Type':'Integer', 'Number': 'A'},
    {'ID': 'NCITE', 'Description': 'Number of citations', 'Type':'Integer', 'Number': 'A'}
]

FIELDS = {
    "transcipt_consequence": "NTCSQ",
    "regulatory_consequence": "NRCSQ",
    "gene": "NGENE",
    "variant_phenotype": "NVPHN",
    "gene_phenotype": "NGPHN",
    "citation": "NCITE"
}

SKIP_CONSEQUENCE = [
    "downstream_gene_variant",
    "upstream_gene_variant",
    "intergenic_variant",
    "TF_binding_site_variant",
    "TFBS_ablation",
    "TFBS_amplification"
]

def parse_args(args = None, description: bool = None):
    parser = argparse.ArgumentParser(description = description)
    
    parser.add_argument(dest="input_file", type=str, help="input VCF file")
    parser.add_argument('-O', '--output_file', dest="output_file", type=str)
    
    return parser.parse_args(args)

def main(args = None):
    args = parse_args(args)

    input_file = args.input_file
    output_file = args.output_file or input_file.replace(".vcf.gz", "_renamed.vcf.gz")

    with bgzf.open(output_file, "wt") as o_file:
        input_vcf = VCF(input_file)

        # add to header and write header to output vcf
        for header in HEADERS:
            input_vcf.add_info_to_header(header)
        o_file.write(input_vcf.raw_header)

        # parse csq header and get index of each field
        csq_list = input_vcf.get_header_type("CSQ")['Description'].split("Format: ")[1].split("|")
        csq_header_idx = {}
        for index, value in enumerate(csq_list):
            csq_header_idx[value] = index
        
        # iterate through the file
        for variant in input_vcf:
            csqs = variant.INFO["CSQ"]
            items_per_allele = {}
            for csq in csqs.split(","):
                csq_values = csq.split("|")
                
                allele = csq_values[csq_header_idx["Allele"]]
                if allele not in items_per_allele:
                    items_per_allele[allele] = {item: set() for item in FIELDS}

                consequences = csq_values[csq_header_idx["Consequence"]]
                feature_stable_id = csq_values[csq_header_idx["Feature"]]
                for csq in consequences.split("&"):
                    if csq in SKIP_CONSEQUENCE:
                        continue

                    # genes
                    gene = csq_values[csq_header_idx["Gene"]]               
                    items_per_allele[allele]["gene"].add(gene)

                    # regualtory and transcript consequences
                    if csq.startswith("regulatory"):
                        items_per_allele[allele]["regulatory_consequence"].add(f"{feature_stable_id}:{csq}")
                    else:
                        items_per_allele[allele]["transcipt_consequence"].add(f"{feature_stable_id}:{csq}")

                # phenotype
                phenotypes = csq_values[csq_header_idx["PHENOTYPES"]]
                for phenotype in phenotypes.split("&"):
                    pheno_fields = phenotype.split("+")
                    if len(pheno_fields) != 3:
                        continue
                    
                    (name, source, feature) = pheno_fields

                    if feature.startswith("ENS"):
                        items_per_allele[allele]["gene_phenotype"].add(f"{name}:{source}:{feature}")
                    else:
                        items_per_allele[allele]["variant_phenotype"].add(f"{name}:{source}:{feature}")

                # citations
                citations = csq_values[csq_header_idx["PUBMED"]]
                for citation in citations.split("&"):
                    items_per_allele[allele]["citation"].add(citation)

            # create summary info fields
            for field in FIELDS:
                field_nums = []
                for allele in items_per_allele:
                    field_len = len(items_per_allele[allele][field])
                    if field_len > 0:
                        field_nums.append(str(field_len)) 

                if field_nums:
                    variant.INFO[FIELDS[field]] = ",".join(field_nums)

            o_file.write(str(variant))
            
        input_vcf.close()
    
if __name__ == "__main__":
    sys.exit(main())
