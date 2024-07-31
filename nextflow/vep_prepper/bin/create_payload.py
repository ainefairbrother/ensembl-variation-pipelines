import json

payload = {}

payload["user"] = "nakib"
payload["genome_uuid"] = "TBD"
payload["name"] = "TBD"
payload["description"] = "GFF file required to run Ensembl VEP"
payload["dataset_type"] = "vep"

payload["dataset_source"] = {}
payload["dataset_source"]["name"] = 
payload["dataset_source"]["type"] = "TBD"

print(json.dumps(payload, indent = 4))