# OTU_Phylogeny
R scripts used to demonstrate delta transformation technique informing OTU boundaries.

`Intraspecific_SeqError_Sim.R` includes code for the function used to append branches onto existing phylogenies, 
using sim.bd.tree.

`MPDvs*.R` scripts includes code for functions used to transform phylogenies based on corresponding Pagel transformation and 
(given the corresponding comparative.comm object) calculate the mpd of transformed outputs, 
given a specified vector of transform values.

`SimulateCommunity.R` is the primary script, and utilizes the commands created in the above two scripts, 
along with `sim.meta.phy.comm`, to simulate a community with a phylogeny, and then transform that phylogeny according to the
`phy.d.transform` function. Communities/phylogenies modeling populations (intra) as well as random sequencing error (seq) are
generated from the original community/phylogeny, and are similarly transformed. The `SimulateCommunity` function is the 
primary function used to generate simulation data.

`wrapper.R` wraps up these prior scripts and calls them in a simple script (i.e. `source(MPDvsDelta.R)`. It then iterates the 
`SimulateCommunity` function given specified simulation parameters (the params variable) using mcMap 
(to utilize multiple cores). 

`MeasureShiftingDiversity.R` analyzes the output generated from `wrapper.R` (which is stored in an .Rdata file).
It contains the code for two separate functions, which are applied over the entirety of one set of MPD matrices contained
in the simulation output (there are MPD matrices for orig, intra, and seq phylogenies). 

