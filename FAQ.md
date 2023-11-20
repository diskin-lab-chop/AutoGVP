## Frequently Asked Questions (FAQ)

### AutoGVP is running too slow. How can a user improve this run-time?
The runtime is machine specific, but there are multiple ways to increase this You can try to add more stringent filters for your inputted vcf when running 01-filter_vcf.sh. Currently, the default is based on the FILTER column. Others criteria such as read depth and allele frequencies could also be used with the filter_criteria parameter.

### Why recommending lift over the gene symbols?
This may impact gene specific variants. We recommend this as to avoid that

### Recommended conflict_res option?
We have these two as an option as to allow the user to tailor to their analysis. Although, we have not seen major differences, there may be some more specific cases where this matter.  As a default we have “latest” based on the assumption that technology and methods get better over time at variant predictions.

### Why do you have a CAVATICA version?
CAVATICA is a public cloud-based tool that uses our workflow. It has its own pre-processing workflow documented here. It is intended for users that want to apply the standard KidsFirst upstream workflows. However if the user wants to use a specific clinVar version and customize their filters, we suggest using the customize workflow option.

 
