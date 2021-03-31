# OTU_Phylogeny
Scripts used for pipeline analyzing effects of species boundaries on community diversity.

`Intraspecific_SeqError_Sim.R` codes the function used to append branches onto existing phylogenies, 
using sim.bd.tree.

`MPDvsDelta.R` codes for function transforming phylogenies based on Pagel's delta and 
(given the corresponding comparative.comm object) calculate the MPD of transformed outputs, 
given a specified vector of values.

`SimulateCommunity.R` is the primary script, and above two commands and `sim.meta.phy.comm` 
to simulate communities/phylogenies, and to transform phylogenies using
`phy.d.transform`. Communities/phylogenies modeling populations (intra) as well as random sequencing error (seq) are
generated from the original community/phylogeny, and are similarly transformed. Also in this script is the 
`abundance.mapping` function, which transfer community abundances from original community to 
communities with populations and sequencing error.

The object generated from a single simulation instance contains 3 community matrices, 3 phylogenies, 3 MPD.matrices (i.e. 
MPD values recorded after delta transforms, for original phylogenies, phylogenies with intraspecific error, and phylogenies 
with sequencing error), and a vector of original MPD values prior to any branch additions or transformations.

`wrapper.R` wraps up these prior scripts and calls them in a simple script (i.e. `source(MPDvsDelta.R)`. It then iterates the 
`SimulateCommunity` function given specified simulation parameters (the params variable) using mcMap 
(to utilize multiple cores). Currently, the only parameters being altered are intraspecific and sequencing error
birth/death rates, and the initial number of species in the community.

`MeasureShiftingDiversity.R` analyzes the output generated from `wrapper.R` (which is stored in an .Rdata file).
It contains the code for calculating 2 different response metrics from simulation data: correlations and rank shifts. 
Linear models are used to measure the effects of delta and other simulation parameters on each response metric. 
This script also contains code for plotting.

