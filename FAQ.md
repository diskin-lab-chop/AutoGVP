## Frequently Asked Questions (FAQ)

#### AutoGVP is running too slow. How can a user improve this run-time?
The runtime varies and is machine specific, but there are multiple ways to increase this. You can try to add more stringent filters for your vcf when running ```01-filter_vcf.sh```. Currently, the default is based on the ```FILTER``` column. Others criteria such as read depth and allele frequencies could also be used with the ```filter_criteria``` parameter.

#### Why recommending lift over the gene symbols?
This may impact gene specific variant calling. We recommend this as to avoid that and to make sure you have the most recent gene annotations and symbols.

#### Recommended conflict_res option?
We have ```severe``` and ```last``` as options when running ```run_autogvp.sh``` as to allow the user to tailor for their use-case. Although, we have not seen major differences from our testing, there may be some more specific cases where this matters.  As a default we have “latest” based on the assumption that technology and methods get better over time (re: variant predictions). However, sometimes, especially if the latest entries are deposited close to each other in the clinVar database, it may be more appropriate to use the ```severe```. We recommend manual inspection for genes of interest.

#### Why is there a CAVATICA version?
CAVATICA is a open and public cloud-based tool that uses our workflow based on upstream KidsFirst workflows([Germline SNV Annotations](https://github.com/kids-first/kf-germline-workflow/blob/v0.4.4/docs/GERMLINE_SNV_ANNOT_README.md), [Pathogenicity-Preprocessing](https://github.com/d3b-center/D3b-Pathogenicity-Preprocessing)). It is intended for users that want to apply the standard KidsFirst workflows. However if the user wants to use a specific clinVar version and customize their own filters, we suggest using the customize workflow option.
