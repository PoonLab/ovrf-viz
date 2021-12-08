
# Major Issues

## Reviewer 1

1. Annotation reliability

>The authors rely on NCBI annotations of CDSs. However, these tend to be quite inaccurate (as the authors discovered for bovine adenovirus NC 002513, p20 lines 5-8). Overlapping and other non-standard genes (ribosomal frameshifting, non-AUG initiation, transcriptional slippage, alternative splicing, misannotation of an upstream in-frame AUG when in fact it is not actually used due to being upstream of a transcript start site, etc) are common in viruses but often dealt with poorly by automatic gene prediction software or indeed by manual annotators. Thus, many virus sequences have missing ORFs, while some sequences may have spurious ORFs annotated as if they were real CDSs.
>A few misannotations could seriously bias results - particularly for reverse-frame overlapping genes - for example, the authors 4 reverse-frame overlaps in dsRNA viruses arise from "NC\_003729" Reoviridae, "NC\_042073" Cystoviridae, "NC\_043677" Chrysoviridae, "NC\_043678" Chrysoviridae and are likely all misannotations.
>The authors should carefully curate their input dataset and/or show that their methods and conclusions are robust with respect to misannotations.


Thank you for raising this issue.  It is a true and unfortunate fact that many records in Genbank are misannotated, even for entries in the RefSeq database that have presumably been manually curated.  Due to the scale of our analysis (all available virus genomes), it was not feasible to curate the entire data set of 12,609 genomes (451,228 annotated reading frames).  However, there are some quality control measures that can be automated.

For instance, we have been screening for misannotations by testing whether the length of overlapping regions was consistent with apparent frameshift of the respective reading frames.  For example, an overlap between reading frames that are shifted by `+2` should have a length that is not divisible by 3, with remainder 1, *i.e.*, length mod 3 = 1.  Applying this filter to our data identified 5,241 (~1%) entries with inconsistent overlap lengths - these entries were discarded from further analysis.  We manually determined that these cases were frequently caused by undocumented splicing or slippage.  Thus, the cases of misannotation that the reviewer astutely identified have already been filtered out of our analysis, though those entries remained in the source data.  Due to space constraints, this quality control step was not described in our original manuscript.  We sincerely apologize for this misunderstanding, and we have added this explanation to the revision.

To further evaluate the sensitivity of our results to misannotation, we have performed a simulation experiment in which a random sample of reading frame coordinates were modified.
Specifically, we selected a random sample of 5% of entries in the CSV file containing all coding sequences for viruses in the database, and then shifted the start and end coordinates by an amount drawn from a uniform distribution from -10 to 10.
Next, we re-ran the analysis pipeline on the modified data to assess the impact of this error variation on our results.



![](overlap_lengtht_noise.pdf)
Associations between overlap lengths and frame shifts after introduction of random noise to 10\% of the entries of the database


![](treemap_noise.pdf)
Distribution of frameshits across Baltimore classes with the modified database


2. Overlaps in the $+0$ frame shift

> The authors include +0 frame overlaps. By this, they mean where the same coding sequence forms part of two different CDS annotations, e.g. in Venezuelan equine encephalitis alphavirus NC\_001449, polyprotein P123 is encoded in the CDS 44-5684 and polyprotein P1234 is encoded in the CDS 44-7526 with readthrough of the stop codon at position 5682-5684. The authors count the region 44-5684 as a +0 frame overlap with the region 44-7526. The authors justification for counting this as an overlap is "Even though the reading frames share codons, they yield different gene products such that the codons are exposed to different selective environments." and "+0 overlaps increase the selective burden of the same nucleotides" (p19 lines 6-10). I find this difficult to swallow. First, the shared domains often have largely similar functions in the two products (e.g. alphavirus P123 versus the P123 domains of P1234) - cf. different-frame overlaps such as orthoreovirus sigma-1 and sigma-1s - where the two proteins have no sequence in common and hence no shared function at all. Second, all codons are subject to multiple selective pressures anyway (translation, mutational, protein function - many proteins are multifunctional); it doesn't make sense to me to single out +0 overlaps as being subject to additional selection pressures.

We concede that sequences belonging to this +0 category of overlaps do not necessarily result in gene products that experience different selective environments.
Hence, we have removed this justification from the manuscript, as well as the corresponding results from figures and text.
We submit, however, that these cases should not be entirely disregarded.
For instance, they represent a convenient control group --- the same sequence being used to produce different gene products, but in the same reading frame --- that provides context for interpreting results for standard overlapping reading frames.
Therefore, we are providing supplementary material that includes the +0 cases so that this information is available, with the caveats raised by the reviewer.


3. Biases on +0 frame shift

> A second problem with including +0 frame overlaps is that their presence in a genome is very much subject to the whim of the annotator. As an example, the alphavirus NC\_001449 has both readthrough and frameshift ORFs annotated:
> 45-5684 and 45-7526 (stop codon readthrough at 5682-5684)
7562-11329 and join(7562-9970,9970-10047) (ribosomal frameshift at 9970)
> leading to +0 frame overlaps 45-5684 and 7562-9970. On the other hand, the alphavirus NC\_023812 only has two CDSs annotated 43-7467 and 7541..11269. The stop codon read through and frameshift are there, but extra CDSs are not annotated so the authors find zero +0 frame overlaps.
> Such annotation issues will bias the results as these mainly affect RNA viruses.

We agree that this is a significant limitation in the analysis (see above, removed +0 from results). As kindly suggested by the revieer, when re-plotting our figures, we see a change in the distribution of overlapping lengths, particularly in RNA viruses. 


4. Annotation bug causing wrong overlap measurement

> There seem to be issues with the annotation scripts. For example in the authors' supplementary overlaps.csv file, we have the following entries for NC\_001449, listing the +0 frame P123/P1234 overlap in the first line, but then the frameshift ORF join(7562-9970,9970-10047), comprising a +0 frame overlap 7562-9970 and a +2 frame overlap 9970-10047, is listed as two +0 frame overlaps:
> "NC\_001449","non-structural polyprotein precursor P123",44,5684,1,"non-structural polyprotein precursor P1234",44,7526,1,5640,7482,5640,+0 \\
"NC\_001449","truncated polyprotein",7561,10047,1,"structural polyprotein precursor",7561,11329,1,2487,3768,2409,+0 \\
"NC\_001449","truncated polyprotein",7561,10047,1,"structural polyprotein precursor",7561,11329,1,2487,3768,78,+0 \\
> Another example is in simian hemorrhagic fever virus NC\_003092, where the 678-nt +1 frame fragment (i.e. 2875-3552) of the join(210-2876,2875-3552) frameshift CDS is labelled as a +0 frame overlap.
"NC\_003092","ORF1aTF polyprotein",209,3552,1,"ORF1a polyprotein"MaSc,209,6527,1,3345,6318,2667,+0 \\
"NC\_003092","ORF1aTF polyprotein",209,3552,1,"ORF1a polyprotein",209,6527,1,3345,6318,678,+0 \\

Thank you for raising this issue. A closer look into it lead us to conclude that such cases appear when we try to calculate overlaps in coding sequences with internal ribosomal frameshifts. 

For example, "ORFS1ab protein" and "ORF1aTF polyprotein" have the following coordinates:  209:6521;6520:10996 and 209:2876;2874:3552 respectively. 

When calculating the overlaps of such reading frames, we are storing different overlap lengths for the different parts of the CDS as independent entries. We recognized a total of 1,426 these cases, and merge them to have one single entry per CDS with the length of the total overlap. 

However, as the reviewer pointed out, many of such cases are present on +0 frameshifts. Therefore, they where removed from the subsequent analysis.

## Reviewer #2

1. Consistency of alignment-free method
> In the graph analysis, homology is based on an alignment free method. how robust are the results of the graph analyses? Will similar conclusions be made with an alignment based homology method? How sensitive are the results and conclusions to the chosen cutoffs? To answer these questions it would be of value to repeat the graph analysis (at least on a small subset) with an alignment based method for homology clustering and possibly also modifying clustering cutoffs on the existing clustering and briefly discuss the impact of these changes on the results and conclusions of the analyses.


The creation of the adjacency graphs using an alignment-free method was motivated for previous difficulties that we encountered when aligning entire virus families. Such alignments tend to be affected for large periods of evolution that in viruses usually involve events such as recombination, insertion and deletions, and the lack of entire genes in different species.
We tried to reproduce the analysis in a small family such as Papillomaviridae but were not able to create a proper alignment from which to estimate alignment scores that would resemble the kmer-counts that we obtained with the clustering method.
* provide more information about problems with alignment


2. Overlap length information on graphs

> The graph analysis is lacking information of the overlap frame and overlap length, which are of great interest. Adding edges separate for frame and for length in the supplementary information would be of value (or at least adjusting edge width by a minimum of overlap length and discussing the impact on the resulting graph).

* we already use edge width for number of proteins (sample size)
* graphs are already saturated with information
* manually annotate graphs with edge labels for overlap lengths
* Discussion - briefly talk about how graphs could be modified to present other information, provide example as supplementary figure (show one)


## Minor Issues

### Reviewer 1

1. Suggested references
> While the authors have referenced previous studies such as Brandes & Linial, Chirico et al, and Schlub & Holmes, I feel that they should also reference Rancurel et al PMID 19640978 and Pavesi et al PMID 30339683, who used carefully curated datasets of overlaps rather than relying on NCBI annotation.

Include them


2. Specify terms
> The authors frequently (e.g. p5 line 24, p5 line 25, etc) use the term "reading frame" when the term "open reading frame" (i.e. "ORF") would, I feel, be clearer. A genome has 6 (potential) reading frames (viz. +0,+1,+2,-0,-1,-2) but (for large genomes) can have 100s of ORFs.

change them


3. Negative strand

> On p6 lines 10-11, I feel it would be useful to state explicitly how (-)strand ORFs (where the biological 5' and 3' ends are in the opposite orientation with respect to the genome coordinates) are dealt with.



4. Coronaviridae genome

> In Fig 4, why are there only 11 clusters in the coronaviridae graph, given that SARS alone has 13 ORFs and there are additional coronavirus accessory genes present in other coronaviruses?

Proteins with high homology to certain clussters can be grouped in that specific cluster and won't form their individual cluster.


5. Simian hemorrhagic fever genome

> The authors say that simian hemorrhagic fever virus genome (NC 003092) "encodes 33 reading frames" (p9 line 12). I don't know where this number came from because NC\_003092 only has ~15 CDSs.


Double check it


### Reviewer 2

1. Figures discrepancy

> There seems to be a discrepancy between figure 1c and figure 2a - in overlap of frame +2 the majority of overlap lengths is 1 and 4, and overlap of +2 is the dominant in both dsDNA and dsRNA, however in Figure 1c there are no peaks for dsDNA with length 1 and 4, this does not make sense because there are more dsDNA than dsRNA in the analysis. This apparent discrepancy needs to be fixed or well explained in the text.



2. Overlaps on -0 frame shift

> Long overlaps with -0 frame are not explained - is this unique to a single edge between homology clusters or does this occur in several edges, if the latter is true this might have implications for the evolution of the genetic code.


3. Conserved overlaps

> Although most adjacent nodes don't show overlap, some large nodes have thick edges (e.g., nodes 6 and 7 in the Adenoviridae family) suggesting that in some cases adjacency of specific genes is accompanied by conserved overlap between the two genes - this should be discussed.

