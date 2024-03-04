#!/usr/bin/env python3

import sys
import os
from cyvcf2 import VCF, Writer
from Bio import bgzf
import argparse

HEADERS = [
    {'ID': 'RAF', 'Description': 'Allele frequencies from representative population', 'Type':'Float', 'Number': 'A'},
    {'ID': 'NTCSQ', 'Description': 'Number of regulatory consequences', 'Type':'Integer', 'Number': '1'},
    {'ID': 'NRCSQ', 'Description': 'Number of transcripts consequences', 'Type':'Integer', 'Number': '1'},
    {'ID': 'NGENE', 'Description': 'Number of overlapped gene', 'Type':'Integer', 'Number': '1'},
    {'ID': 'NVPHN', 'Description': 'Number of associated variant-linked phenotypes', 'Type':'Integer', 'Number': 'A'},
    {'ID': 'NGPHN', 'Description': 'Number of associated gene-linked phenotypes', 'Type':'Integer', 'Number': 'A'},
    {'ID': 'NCITE', 'Description': 'Number of citations', 'Type':'Integer', 'Number': '1'}
]

PER_ALLELE_FIELDS = {
    "variant_phenotype": "NVPHN",
    "gene_phenotype": "NGPHN",
}

PER_VARIANT_FIELDS = {
    "transcipt_consequence": "NTCSQ",
    "regulatory_consequence": "NRCSQ",
    "gene": "NGENE",
    "citation": "NCITE"
}

FREQUENCY_FIELD = "RAF"
# [csq_field, diplay_name]
FREQUENCY_META = {
    "homo_sapiens": {
        "GRCh38": ["gnomAD_genomes_AF", "gnomAD genomes v3.1.2"],
        "GRCh37": ["gnomAD_exomes_AF", "gnomAD exomes v2.1.1"]
    }
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
    
    parser.add_argument(dest="species", type=str, help="species production name")
    parser.add_argument(dest="assembly", type=str, help="assembly default")
    parser.add_argument(dest="input_file", type=str, help="input VCF file")
    parser.add_argument('-O', '--output_file', dest="output_file", type=str)
    
    return parser.parse_args(args)

def minimise_allele(ref: str, alt: str) -> str:
    minimised_allele_string = alt
    if ref[0] == alt[0]:
        minimised_allele_string = alt[1:] if len(alt) > 1 else "-" 
    return minimised_allele_string

def main(args = None):
    args = parse_args(args)

    species = args.species
    assembly = args.assembly
    input_file = os.path.realpath(args.input_file)
    output_file = args.output_file or os.path.join(os.path.dirname(input_file), "UPDATED_SS_" + os.path.basename(input_file))

    # frequency meta
    (freq_csq_field, freq_info_display) = (None, "")
    if species in FREQUENCY_META and assembly in FREQUENCY_META[species] and len(FREQUENCY_META[species][assembly]) == 2:
        (freq_csq_field, freq_info_display) = FREQUENCY_META[species][assembly]

    input_vcf = VCF(input_file)

    # add to header and write header to output vcf
    if freq_info_display != "":
        HEADERS[0]['Description'] = HEADERS[0]['Description'] + f" ({freq_info_display})"
    for header in HEADERS:
        input_vcf.add_info_to_header(header)

    output_vcf = Writer(output_file, input_vcf, mode="w")

    # parse csq header and get index of each field
    csq_list = input_vcf.get_header_type("CSQ")['Description'].split("Format: ")[1].split("|")
    csq_header_idx = {}
    for index, value in enumerate(csq_list):
        csq_header_idx[value] = index
    
    # iterate through the file
    for variant in input_vcf:
        # create minimalized allele order
        allele_order = []
        ref = variant.REF
        for alt in variant.ALT:
            allele_order.append(minimise_allele(ref, alt))

        items_per_variant = {item: set() for item in PER_VARIANT_FIELDS}
        items_per_allele = {}

        # travers through each csq entry
        csqs = variant.INFO["CSQ"]
        for csq in csqs.split(","):
            csq_values = csq.split("|")
            
            allele = csq_values[csq_header_idx["Allele"]]
            if allele not in items_per_allele:
                items_per_allele[allele] = {item: set() for item in PER_ALLELE_FIELDS}
                
            consequences = csq_values[csq_header_idx["Consequence"]]
            feature_stable_id = csq_values[csq_header_idx["Feature"]]
            for csq in consequences.split("&"):
                if csq in SKIP_CONSEQUENCE:
                    continue

                # genes
                gene = csq_values[csq_header_idx["Gene"]]               
                items_per_variant["gene"].add(gene)

                # regualtory and transcript consequences
                if csq.startswith("regulatory"):
                    items_per_variant["regulatory_consequence"].add(f"{feature_stable_id}:{csq}")
                else:
                    items_per_variant["transcipt_consequence"].add(f"{feature_stable_id}:{csq}")
            
            csq_values_len = len(csq_values)

            # phenotype
            phenotype_csq_idx = csq_header_idx["PHENOTYPES"]
            if phenotype_csq_idx < csq_values_len:
                phenotypes = csq_values[phenotype_csq_idx]
                for phenotype in phenotypes.split("&"):
                    pheno_PER_ALLELE_FIELDS = phenotype.split("+")
                    if len(pheno_PER_ALLELE_FIELDS) != 3:
                        continue
                    
                    (name, source, feature) = pheno_PER_ALLELE_FIELDS
                    if feature.startswith("ENS"):
                        items_per_allele[allele]["gene_phenotype"].add(f"{name}:{source}:{feature}")
                    else:
                        items_per_allele[allele]["variant_phenotype"].add(f"{name}:{source}:{feature}")

            # citations
            pubmed_csq_idx = csq_header_idx["PUBMED"]
            if pubmed_csq_idx < csq_values_len:
                citations = csq_values[pubmed_csq_idx]
                for citation in citations.split("&"):
                    if citation != "":
                        items_per_variant["citation"].add(citation)

            # frequency
            if freq_csq_field:
                af_csq_idx = csq_header_idx[freq_csq_field]
                if af_csq_idx < csq_values_len:
                    frequency = csq_values[af_csq_idx]
                    if frequency != "":
                        items_per_allele[allele]["frequency"] = frequency

        # create summary info for per allele fields
        for field in PER_ALLELE_FIELDS:
            field_nums = []
            for allele in allele_order:
                if allele in items_per_allele and field in items_per_allele[allele]:
                    field_len = len(items_per_allele[allele][field])
                    if field_len > 0:
                        field_nums.append(str(field_len)) 

            if field_nums:
                variant.INFO[PER_ALLELE_FIELDS[field]] = ",".join(field_nums)

        # create summary info for frequency
        field_vals = []
        for allele in allele_order:
            if allele in items_per_allele and "frequency" in items_per_allele[allele]:
                field_vals.append(items_per_allele[allele]["frequency"])
            else:
                field_vals.append(".")

            if not all(freq == "." for freq in field_vals):
                variant.INFO[FREQUENCY_FIELD] = ",".join(field_vals)

        # create summary info for per variant fields
        for field in PER_VARIANT_FIELDS:
            field_len = len(items_per_variant[field])
            if field_len > 0:
                variant.INFO[PER_VARIANT_FIELDS[field]] = str(field_len) 

        output_vcf.write_record(variant)
        
    input_vcf.close()
    output_vcf.close()
    
if __name__ == "__main__":
    sys.exit(main())
