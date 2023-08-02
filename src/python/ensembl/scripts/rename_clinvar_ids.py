import sys
from cyvcf2 import VCF, Writer

input_vcf = VCF(sys.argv[1])
ouput_file = sys.argv[2]
w = Writer(ouput_file, input_vcf)
for variant in input_vcf:
    current_id = variant.ID
    
    if not current_id.startswith("VCV"):
        leading_zero = 9 - len(current_id)
        new_id = "VCV" + ("0" * leading_zero) + current_id
        
        variant.ID = new_id       
    
    w.write_record(variant)
w.close()
input_vcf.close()