table vep_consequences
"Consequences from VEP"
(
string	chrom;  "Reference sequence chromosome or scaffold"
uint    chromStart;     "Start position of feature on chromosome"
uint    chromEnd;       "End position of feature on chromosome"
string	id;		"ID of the variant"
string	class;		"Class of the variant"
lstring	ref;   		"Reference allele"
lstring	alts;		"Alternative allele(s)"
uint  variantGroup;	"Id of the variant group the variant belongs to"
string	consequence;	"Most severe consequence of the variant"
string  spdi;	"SPDI notation of the variant"
)
