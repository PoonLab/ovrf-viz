
## Major Issues

### Reviewer 1


> 1) The authors rely on NCBI annotations of CDSs. However, these tend to be quite inaccurate (as the authors discovered for bovine adenovirus NC 002513, p20 lines 5-8). Overlapping and other non-standard genes (ribosomal frameshifting, non-AUG initiation, transcriptional slippage, alternative splicing, misannotation of an upstream in-frame AUG when in fact it is not actually used due to being upstream of a transcript start site, etc) are common in viruses but often dealt with poorly by automatic gene prediction software or indeed by manual annotators. Thus, many virus sequences have missing ORFs, while some sequences may have spurious ORFs annotated as if they were real CDSs.
>A few misannotations could seriously bias results - particularly for reverse-frame overlapping genes - for example, the authors 4 reverse-frame overlaps in dsRNA viruses arise from "NC\_003729" Reoviridae, "NC\_042073" Cystoviridae, "NC\_043677" Chrysoviridae, "NC\_043678" Chrysoviridae and are likely all misannotations.
>The authors should carefully curate their input dataset and/or show that their methods and conclusions are robust with respect to misannotations.


Thank you for raising this issue.  It is a true and unfortunate fact that many records in Genbank are misannotated, even for the reference genomes used in our study, which were manually curated as an output of the NCBI Viral genomes Resource.  Due to the scale of our analysis (all available virus genomes), it was not feasible to curate the entire data set of 12,609 genomes (451,228 annotated reading frames).  However, there are some quality control measures that can be automated.

For instance, we have been screening for misannotations by testing whether the length of overlapping regions was consistent with apparent frameshift of the respective reading frames.  For example, an overlap between reading frames that are shifted by `+2` should have a length that is not divisible by 3, with remainder 1, *i.e.*, length mod 3 = 1.  Applying this filter to our data identified 5,241 (~1%) entries with inconsistent overlap lengths - these entries were discarded from further analysis.  We manually determined that these cases were frequently caused by undocumented splicing or slippage.  Thus, the cases of misannotation that the reviewer astutely identified have already been filtered out of our analysis, though those entries remained in the source data.  Due to space constraints, this quality control step was not described in our original manuscript.  We sincerely apologize for this misunderstanding, and we have added this explanation to the revision.

To further evaluate the sensitivity of our results to misannotation, we have performed a simulation experiment in which a random sample of reading frame coordinates were modified.
Specifically, we selected a random sample of 5% of entries in the CSV file containing all coding sequences for viruses in the database, and then shifted the start and end coordinates by an amount drawn from a uniform distribution from -10 to 10.
Next, we re-ran the analysis pipeline on the modified data to assess the impact of this error variation on our results.
**IN PROGRESS**


![](overlap_lengtht_noise.pdf)
Associations between overlap lengths and frame shifts after introduction of random noise to 10\% of the entries of the database


![](treemap_noise.pdf)
Distribution of frameshits across Baltimore classes with the modified database



> 2) The authors include +0 frame overlaps. By this, they mean where the same coding sequence forms part of two different CDS annotations, e.g. in Venezuelan equine encephalitis alphavirus NC\_001449, polyprotein P123 is encoded in the CDS 44-5684 and polyprotein P1234 is encoded in the CDS 44-7526 with readthrough of the stop codon at position 5682-5684. The authors count the region 44-5684 as a +0 frame overlap with the region 44-7526. The authors justification for counting this as an overlap is "Even though the reading frames share codons, they yield different gene products such that the codons are exposed to different selective environments." and "+0 overlaps increase the selective burden of the same nucleotides" (p19 lines 6-10). I find this difficult to swallow. First, the shared domains often have largely similar functions in the two products (e.g. alphavirus P123 versus the P123 domains of P1234) - cf. different-frame overlaps such as orthoreovirus sigma-1 and sigma-1s - where the two proteins have no sequence in common and hence no shared function at all. Second, all codons are subject to multiple selective pressures anyway (translation, mutational, protein function - many proteins are multifunctional); it doesn't make sense to me to single out +0 overlaps as being subject to additional selection pressures.

We concede that sequences belonging to this +0 category of overlaps do not necessarily result in gene products that experience different selective environments.
Hence, we have removed this justification from the manuscript, as well as the corresponding results from figures and text.
We submit, however, that these cases should not be entirely disregarded.
For instance, they represent a convenient control group --- the same sequence being used to produce different gene products, but in the same reading frame --- that provides context for interpreting results for standard overlapping reading frames.
Therefore, we are providing supplementary material that includes the +0 cases so that this information is available, with the caveats raised by the reviewer.
**IN PROGRESS**

> 3) A second problem with including +0 frame overlaps is that their presence in a genome is very much subject to the whim of the annotator. As an example, the alphavirus NC\_001449 has both readthrough and frameshift ORFs annotated:
> 45-5684 and 45-7526 (stop codon readthrough at 5682-5684)
7562-11329 and join(7562-9970,9970-10047) (ribosomal frameshift at 9970)
> leading to $+0$ frame overlaps 45-5684 and 7562-9970. On the other hand, the alphavirus NC\_023812 only has two CDSs annotated 43-7467 and 7541..11269. The stop codon read through and frameshift are there, but extra CDSs are not annotated so the authors find zero +0 frame overlaps.
> Such annotation issues will bias the results as these mainly affect RNA viruses.

We agree that inconsistent annotation of +0 frame overlaps was a significant issue in our original analysis.  As noted above, we have removed the +0 overlaps from our main results in the manuscript, and relegated the previous versions in Supplementary Material.
**DONE**

> 4) There seem to be issues with the annotation scripts. For example in the authors' supplementary overlaps.csv file, we have the following entries for NC\_001449, listing the +0 frame P123/P1234 overlap in the first line, but then the frameshift ORF join(7562-9970,9970-10047), comprising a +0 frame overlap 7562-9970 and a +2 frame overlap 9970-10047, is listed as two +0 frame overlaps:
> "NC\_001449","non-structural polyprotein precursor P123",44,5684,1,"non-structural polyprotein precursor P1234",44,7526,1,5640,7482,5640,+0 \\
"NC\_001449","truncated polyprotein",7561,10047,1,"structural polyprotein precursor",7561,11329,1,2487,3768,2409,+0 \\
"NC\_001449","truncated polyprotein",7561,10047,1,"structural polyprotein precursor",7561,11329,1,2487,3768,78,+0 \\
> Another example is in simian hemorrhagic fever virus NC\_003092, where the 678-nt +1 frame fragment (i.e. 2875-3552) of the join(210-2876,2875-3552) frameshift CDS is labelled as a +0 frame overlap.
"NC\_003092","ORF1aTF polyprotein",209,3552,1,"ORF1a polyprotein"MaSc,209,6527,1,3345,6318,2667,+0 \\
"NC\_003092","ORF1aTF polyprotein",209,3552,1,"ORF1a polyprotein",209,6527,1,3345,6318,678,+0 \\


Thank you for raising this issue.
As noted in our response to an earlier point, cases that involved splices such as the example provided by the reviewer were excluded from our analysis.
**IN PROGRESS?**

This was the correct file. How did we deal with those cases? Removing overlaps product from splicing events.
The input for our calculation of overlaps is based on the file total\_orfs.csv that contains for each protein from each genome the accession number (accno), Product, Strand, Coordinates and Start Codon.
The misscalculation problem seems to be associated with entries that have internal frame shifts.
For example, for the record NC\_003092, the annotations for ORFS1ab protein and ORF1aTF polyprotein have some internal frame shifts that are not even spliced events: 209:6521;6520:10996 and 209:2876;2874:3552 respectively.
Similarly, for NC\_001449 the truncated polyprotein has coordinates: 7561:9970;9969:10047 (again, a frame shift of one nucleotide) where ribosome probably jumps back one nucleotide.
Calculating frame shift on such entries is challenging since there would be two different frame shifts for the same protein.

Following step: Detect how many entries we have with similar characteristics.



### Reviewer #2

> 1. In the graph analysis, homology is based on an alignment free method. how robust are the results of the graph analyses? Will similar conclusions be made with an alignment based homology method? How sensitive are the results and conclusions to the chosen cutoffs? To answer these questions it would be of value to repeat the graph analysis (at least on a small subset) with an alignment based method for homology clustering and possibly also modifying clustering cutoffs on the existing clustering and briefly discuss the impact of these changes on the results and conclusions of the analyses.


The creation of the adjacency graphs using an alignment-free method was motivated for previous difficulties that we encountered when aligning entire virus families.  Such alignments tend to be affected for large periods of evolution that in viruses usually involve events such as recombination, insertion and deletions, and the lack of entire genes in different species.  
An alternative approach would be to use pairwise alignment to generate a similar distance matrix for all coding sequences in a virus family.
Since alignment is more computationally complex than matching k-mers, we wrote a message-passing interface (MPI)-enabled Python script to run the alignment program MAFFT in parallel (15 cores) and compute p-distances (*i.e.*, the proportion of aligned sites that are different and ungapped).
It was then straightforward to graft the resulting distance matrix into our analysis pipeline.
In the revised manuscript, we provide a set of clustering results for Papillomaviridae using pairwise alignment instead of the k-mer-based method.
**IN PROGRESS**



> 2. The graph analysis is lacking information of the overlap frame and overlap length, which are of great interest. Adding edges separate for frame and for length in the supplementary information would be of value (or at least adjusting edge width by a minimum of overlap length and discussing the impact on the resulting graph).

We concur that our current layout algorithm for drawing adjacency graphs does not visualize information in the frame or length of overlaps.  The key challenge is that there are limited channels for visually communicating information at our disposal, such as node size, edge width, node colour.  We have avoided overloading the graph by, for example, varying both node shape, size and colour.  As suggested by the reviewer, we revised our Python script for generating graph layouts, so that we can now provide examples of scaling edge widths to overlap lengths instead of the number of overlaps as supplementary material.

* Discussion - briefly talk about how graphs could be modified to present other information, provide example as supplementary figure (show one)


## Minor Issues

### Reviewer 1

> While the authors have referenced previous studies such as Brandes & Linial, Chirico et al, and Schlub & Holmes, I feel that they should also reference Rancurel et al PMID 19640978 and Pavesi et al PMID 30339683, who used carefully curated datasets of overlaps rather than relying on NCBI annotation.

Thank you for bringing these articles to our attention.  We have incorporated these references in discussing the trade-off between accuracy and scaleability, with respect to the careful manual curation versus automated computational processing of virus genomes.  


> The authors frequently (e.g. p5 line 24, p5 line 25, etc) use the term "reading frame" when the term "open reading frame" (i.e. "ORF") would, I feel, be clearer. A genome has 6 (potential) reading frames (viz. +0,+1,+2,-0,-1,-2) but (for large genomes) can have 100s of ORFs.

change them


> On p6 lines 10-11, I feel it would be useful to state explicitly how (-)strand ORFs (where the biological 5' and 3' ends are in the opposite orientation with respect to the genome coordinates) are dealt with.


> In Fig 4, why are there only 11 clusters in the coronaviridae graph, given that SARS alone has 13 ORFs and there are additional coronavirus accessory genes present in other coronaviruses?

Proteins with high homology to certain clussters can be grouped in that specific cluster and won't form their individual cluster.


> The authors say that simian hemorrhagic fever virus genome (NC 003092) "encodes 33 reading frames" (p9 line 12). I don't know where this number came from because NC\_003092 only has ~15 CDSs.


Double check it


### Reviewer 2


> There seems to be a discrepancy between figure 1c and figure 2a - in overlap of frame +2 the majority of overlap lengths is 1 and 4, and overlap of +2 is the dominant in both dsDNA and dsRNA, however in Figure 1c there are no peaks for dsDNA with length 1 and 4, this does not make sense because there are more dsDNA than dsRNA in the analysis. This apparent discrepancy needs to be fixed or well explained in the text.



> Long overlaps with -0 frame are not explained - is this unique to a single edge between homology clusters or does this occur in several edges, if the latter is true this might have implications for the evolution of the genetic code.



> Although most adjacent nodes don't show overlap, some large nodes have thick edges (e.g., nodes 6 and 7 in the Adenoviridae family) suggesting that in some cases adjacency of specific genes is accompanied by conserved overlap between the two genes - this should be discussed.

Thank you for raising this point.  We have a discussion of this example in the revised manuscript.
