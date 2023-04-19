import sys
from cyvcf2 import VCF, Writer

input_vcf = VCF(sys.argv[1])
id_count = {}
for variant in input_vcf:
    if not variant.ID in id_count:
        id_count[variant.ID] = False
    else:
        id_count[variant.ID] = True
input_vcf.close()

input_vcf = VCF(sys.argv[1])
ouput_file = sys.argv[1].replace("renamed", "processed")
w = Writer(ouput_file, input_vcf)
for variant in input_vcf:
    chr = variant.CHROM
    if ("CTG" in chr) or ("PATCH" in chr) or ("TEST" in chr):
        continue

    if id_count[variant.ID]:
        variant.ID = "."

    w.write_record(variant)
w.close()
input_vcf.close()