

# The basic usage for microbiome-related data

The general usage for mining the information related to microbiome-related signatures.


```r
library(wcGeneSummary)
library(ggplot2)
library(ggraph)
library(RColorBrewer)
```

We use BugSigDB, and its R port bugsigdbr to obtain the curated dataset of the relationship with bacterial taxonomy and human diseases ([Geistlinger et al. 2022](https://bugsigdb.org/Main_Page)). Users can query microbiome names, which will be searched for MetaPhlAn taxonomic annotation. If `target="title"`, the title of the corresponding articles will be summarized.


```r
basic <- wcBSDB(c("Veillonella dispar","Neisseria flava"),tag=FALSE,
    curate=TRUE,target="title",pre=TRUE,cl=snow::makeCluster(12),
    pal=RColorBrewer::brewer.pal(10, "Set2"),numWords=80,argList=list(min.freq=1))
#> Input microbes: 2
#>   Found 17 entries for Veillonella dispar
#>   Found 1 entries for Neisseria flava
#> Including 28 entries
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
basic@freqDf |> head(n=20)
#>                    word freq
#> gut                 Gut    8
#> oral               Oral    5
#> patients       patients    5
#> study             study    3
#> arthritis     arthritis    2
#> association association    2
#> bacterial     bacterial    2
#> covid19         COVID19    2
#> diabetes       diabetes    2
#> dysbiosis     Dysbiosis    2
#> infant           infant    2
#> infection     Infection    2
#> obese             obese    2
#> respiratory respiratory    2
#> sequencing   sequencing    2
#> surgery         surgery    2
#> tract             tract    2
#>  features      features    1
#> 16s                 16S    1
#> aerosol         Aerosol    1
basic@wc
```

<img src="02-microbiome_usage_files/figure-html/bsdb_basic-1.png" width="768" />

If `target="abstract"`, the corresponding abstract will be fetched and be summarized.


```r
basic2 <- wcBSDB(c("Veillonella dispar","Neisseria flava"),tag=TRUE,
    curate=TRUE,target="abstract",pre=TRUE,cl=snow::makeCluster(12),
    pal=RColorBrewer::brewer.pal(10, "Set2"),numWords=80)
#> Input microbes: 2
#>   Found 17 entries for Veillonella dispar
#>   Found 1 entries for Neisseria flava
#> Including 28 entries
#> Target is abstract
#>   Querying PubMed for 17 pmids
#>   Querying without API key
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
#> Multiscale bootstrap... Done.
basic2@freqDf |> head()
#>                word freq
#> patients   patients   37
#> gut             Gut   30
#> oral           Oral   27
#> bacterial Bacterial   25
#> species     species   24
#> microbial Microbial   23
basic2@wc
```

<img src="02-microbiome_usage_files/figure-html/bsdb_basic2-1.png" width="768" />

For successful visualization, pre-caculated TF-IDF and frequency data frame is available and one can use them to filter the highly occurring words, or the other prefiltering option used in `wcGeneSummary`.


```r
rmwords <- wcGeneSummary:::allFreqBSDB
filter <- rmwords[rmwords$freq>quantile(rmwords$freq, 0.95),]
filter |> head(n=20)
#>             freq        word
#> microbiota   275  microbiota
#> gut          242         gut
#> microbiome   175  microbiome
#> patients     146    patients
#> study         69       study
#> cancer        64      cancer
#> composition   62 composition
#> oral          55        oral
#> human         49       human
#> intestinal    48  intestinal
#> children      45    children
#> microbial     45   microbial
#> disease       43     disease
#> fecal         41       fecal
#> alterations   36 alterations
#> analysis      34    analysis
#> association   34 association
#> dysbiosis     34   dysbiosis
#> infection     34   infection
#> risk          33        risk
```
The network visualization is possible by enabling `plotType="network"`.
The same parameters that can be passed to `wcGeneSummary` can be used.


```r
net <- wcBSDB(c("Neisseria","Veillonella"),
    curate=TRUE,target="title",pre=TRUE,plotType="network",
    additionalRemove=filter$word, corThresh=0.2, edgeLink=FALSE,
    numWords=60)
#> Input microbes: 2
#>   Found 76 entries for Neisseria
#>   Found 221 entries for Veillonella
#> Including 502 entries
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
net@net
```

<img src="02-microbiome_usage_files/figure-html/bsdb_basic_network-1.png" width="768" />

The words-to-species relationship can be plotted by `mbPlot=TRUE`.


```r
net2 <- wcBSDB(c("Veillonella dispar","Neisseria flava",
    "Veillonella parvula","Akkermansia muciniphila"), mbPlot=TRUE,
    curate=TRUE,target="title",pre=TRUE,plotType="network",
    additionalRemove=filter$word,
    numWords=40, corThresh=0.2, colorText=TRUE)
#> Input microbes: 4
#>   Found 17 entries for Veillonella dispar
#>   Found 1 entries for Neisseria flava
#>   Found 20 entries for Veillonella parvula
#>   Found 21 entries for Akkermansia muciniphila
#> Including 90 entries
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
net2@net
```

<img src="02-microbiome_usage_files/figure-html/bsdb_basic_network_mb-1.png" width="768" />

As the BugSigDB contains the relationship between bacterial taxonomy and disease, disease name can also be plotted. When `disPlot=TRUE`, the `mbPlot`
 will be set to `TRUE` by default.
 

```r
net3 <- wcBSDB(c("Veillonella dispar","Neisseria flava",
    "Veillonella parvula","Akkermansia muciniphila"), mbPlot=TRUE,
    curate=TRUE,target="title",pre=TRUE,plotType="network",
    additionalRemove=filter$word, disPlot=TRUE,
    numWords=40, corThresh=0.2, colorText=TRUE)
#> Input microbes: 4
#>   Found 17 entries for Veillonella dispar
#>   Found 1 entries for Neisseria flava
#>   Found 20 entries for Veillonella parvula
#>   Found 21 entries for Akkermansia muciniphila
#> Including 90 entries
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
net3@net
```

<img src="02-microbiome_usage_files/figure-html/bsdb_basic_network_dis-1.png" width="960" />

Other than curated databases, the PubMed query can also be performed with setting `curate=FALSE`. This way, the text information of the latest literature for the microbes and diseases can be plotted. The options for use in function obtaining PubMed information can be specified to `abstArg` in list format, like `sortOrder="pubdate"`.


```r
net4 <- wcBSDB(c("Veillonella dispar","Neisseria flava",
    "Veillonella parvula","Akkermansia muciniphila"), mbPlot=TRUE,
    curate=FALSE,target="title",pre=TRUE,plotType="network",
    additionalRemove=filter$word, disPlot=TRUE,
    numWords=40, corThresh=0.2, colorText=TRUE,
    abstArg = list(retMax=80, sortOrder="pubdate"))
#> Input microbes: 4
#>   Found 17 entries for Veillonella dispar
#>   Found 1 entries for Neisseria flava
#>   Found 20 entries for Veillonella parvula
#>   Found 21 entries for Akkermansia muciniphila
#> Including 90 entries
#> Proceeding without API key
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
net4@net
```

<img src="02-microbiome_usage_files/figure-html/no_curate-1.png" width="960" />

For microbiome analysis, it is often the case that investigating coded enzymes is important. Using `wcEC` function and `getUPtax` function, the queried species or genus can be linked to possible interaction with enzymes using following databases. The downloaded file path should be specified to the function like below to link the queried taxonomy and enzymes. Specifically, enzymes listed in `enzyme.dat` are searched, and corresponding UniProt identifiers are obtained, followed by mapping using `speclist.txt`. This way, the links to microbe - textual information - enzyme can be plotted. 

- [enzyme.dat - from Expasy](https://enzyme.expasy.org/)
- [speclist.txt - UniProt Controlled vocabulary of species](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt)


```r
# Not run
vp <- wcBSDB(c("Veillonella parvula"),
                 plotType="network",
                curate=TRUE, target="title", edgeLink=FALSE,
                mbPlot = TRUE, ecPlot=TRUE, disPlot=TRUE, tag=TRUE,
                cl=snow::makeCluster(10),colorText=TRUE, pre=TRUE, numWords=30,
                nodePal=RColorBrewer::brewer.pal(10, "Set2"),
                ecFile="../enzyme.dat",
                upTaxFile = "../speclist.txt")
#> Input microbes: 1
#>   Found 20 entries for Veillonella parvula
#> Including 31 entries
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
#> Processing EC file
#>   Linking taxonomy to EC
#> Multiscale bootstrap... Done.
vp@net
```

<img src="02-microbiome_usage_files/figure-html/enzyme-1.png" width="960" />

Further, the relationship between metabolites and microbiome is of interest. Recent studies have revealed various associations in gut microbiome composition and human plasma metabolites, as well as in the other environments.

- Wishart DS, Oler E, Peters H, et al. MiMeDB: the Human Microbial Metabolome Database. Nucleic Acids Res. 2023;51(D1):D611-D620. doi:10.1093/nar/gkac868
- Muller E, Algavi YM, Borenstein E. The gut microbiome-metabolome dataset collection: a curated resource for integrative meta-analysis. npj Biofilms and Microbiomes. 2022;8(1):1-7. doi:10.1038/s41522-022-00345-5
- Dekkers KF, Sayols-Baixeras S, Baldanzi G, et al. An online atlas of human plasma metabolite signatures of gut microbiome composition. Nat Commun. 2022;13(1):5370. doi:10.1038/s41467-022-33050-0

We now use data obtained in [Dekkers et al.](https://www.nature.com/articles/s41467-022-33050-0).
First, read the downloaded file using `readxl`.


```r
metab <- readxl::read_excel(
  "../41467_2022_33050_MOESM8_ESM.xlsx",skip = 7)
```

Pass this tibble, as well as the columns to represent `taxonomy`, `metabolites`, and `quantitative values to threshold` to `metCol`. In this case, `metCol <- c("Metagenomic species", "Metabolite", "Spearman's ρ")`.



```r

metabEx <- wcBSDB(c("Akkermansia muciniphila"),
                edgeLink=FALSE,
                curate=TRUE,
                corThresh=0.3,
                nodePal=RColorBrewer::brewer.pal(10, "Dark2"),
                pre=TRUE, tag=TRUE,
                additionalRemove = filter$word,
                target="abstract", colorText=TRUE,
                plotType="network", numWords=50,
                mbPlot=TRUE,
                metCol=c("Metagenomic species", "Metabolite", "Spearman's ρ"),
                metab =metab, metThresh=0.15,
                preserve = TRUE,
                cl=snow::makeCluster(10),
                abstArg = list(retMax=80,
                               sortOrder="relevance"))
#> Input microbes: 1
#>   Found 21 entries for Akkermansia muciniphila
#> Including 31 entries
#> Target is abstract
#>   Querying PubMed for 20 pmids
#>   Querying without API key
#> Filter based on BugSigDB
#> Filtering 0 words (frequency and/or tfidf)
#> Checking metabolites
#> Multiscale bootstrap... Done.
metabEx@net
```

<img src="02-microbiome_usage_files/figure-html/metabex-1.png" width="960" />

In this way, we can plot links between microbes - metabolites - textual information. For all the information combined, one can plot textual information - metabolites - coded enzymes - diseases - microbes link in one query.

For the complex network, the resulting image might be unreadable.
`exportCyjsWithoutImage` function can be used to export the graph to readily interactive interface using `Cytoscape.js`. The below chunk shows the output produced by the function, hosted by GitHub Pages.



```r
## exportCyjsWithoutImage(metabEx@igraph, rootDir=".", netDir="",nodeColorDiscretePal = "Pastel1")
knitr::include_url("https://noriakis.github.io/cyjs_test/complex")
```

<iframe src="https://noriakis.github.io/cyjs_test/complex" width="100%" height="400px" data-external="1" style="border: none;"></iframe>

## Annotating the dendrogram of taxonomy

Annotating the dendrogram (or cladogram) of taxonomy is possible by the `plotEigengeneNetworksWithWords` function. The used need to provide customized data.frame depicting query (taxonomy) and its corresponding text, like MetaCyc.



```r

## Read pathway description
file <- "../../metacyc/24.5/data/pathways.dat"
input <- parseMetaCycPathway(file,
                              candSp="all",
                              withTax=TRUE,
                              clear=TRUE,
                              noComma = TRUE)


## Convert the species name by taxonomizr
input$spConverted <- convertMetaCyc(input$species)
input$spConverted |> head()
#> [1] "Bacteria;Cyanobacteria;NA;Synechococcales;Prochloraceae;Prochloron;Prochloron didemni"                                 
#> [2] "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli"          
#> [3] "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;Salmonella enterica"        
#> [4] "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli"          
#> [5] "Bacteria;Pseudomonadota;Betaproteobacteria;Burkholderiales;Burkholderiaceae;Paraburkholderia;Paraburkholderia phymatum"
#> [6] "Bacteria;Pseudomonadota;Betaproteobacteria;Burkholderiales;Burkholderiaceae;Cupriavidus;Cupriavidus basilensis"
```

Now we fetch species under `Bacillota`.


```r
bc <- input[grepl("Bacillota", input$spConverted),]
bc <- bc[!duplicated(bc$pathwayID),]
bc |> head()
#>      pathwayID
#> 12  PCEDEG-PWY
#> 173   PWY-7872
#> 349   PWY-7013
#> 748   PWY-7886
#> 749   PWY-8146
#> 799   PWY-7942
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               text
#> 12                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                Another group of organisms, called the halorespirators, has been shown to metabolize tetrachloroethene in a process known as dehalorespiration, in which halogenated organic compounds are used as electron acceptors in an anaerobic respiratory process . Halorespiration can result in the complete dechlorination of PCE and TRICHLOROETHENE (TCE) to the benign end products ETHYLENE-CMPD and CPD-9312. In most cases different organisms specialize in different steps of this pathway. Thus far, only TAX-243164 has been shown to completely dechlorinate TETRACHLOROETHENE to ETHYLENE-CMPD through anaerobic reductive dechlorination (dehalorespiration) .  TAX-243164 achieves this utilizing only two enzymes: EC-1.21.99.5 (PCE-RDase), which degrades PCE to TCE, and EC-1.21.99.M1 (TCE-RDase), which performs all the other steps of the pathway. However, TCE-RDase has a very low activity using DICHLOROETHENE and VINYL-CHLORIDE (vinyl chloride, or VC) as substrates, and is apparently catalyzing these reductions as a slow, cometabolic reaction (the reaction rate in &mu;molminmg is 5.2 for TCE but only 0.036 for VC) . This phenomenon poses a problem for the treatment of contaminated sites, as it results in the accumulation of VINYL-CHLORIDE, which is the most toxic compound of all the chloroethenes. The problem can be alleviated by other organisms, such as TAX-311424, which have a dedicated EC-1.21.99.M2, which reduces VC and all dichloroethene (DCE) isomers at a high rate .  All of these enzymes are inhibited by the addition of iodoalkanes, and this inhibition is reversed by exposure to light, indicating the presence of corrinoid cofactors. In addition, most (but not all) of these enzymes have been shown to contain at least one  cluster. 
#> 173                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             The pathway for the degradation of CPD-19877 was found by a bioinformatic analysis of genes clustered with a D-threonate-binding protein (SBP), which is a part of a TRAP (tripartite ATP-independent permease) transporter for the compound .  Following import into the cell, CPD-19877 is phosphorylated to form ERYTHRONATE-4P by EC-2.7.1.220. The enzyme, encoded by the G-55076 gene, is one of the first members of a large family of enzymes, known as DUF (domain of unknown function) 1537, to which a function has been assigned .  Following phosphorylation, the compound is dehygrogenated at position 3 by EC-1.1.1.409. This dehydrogenase, with has some similarity to EC-1.1.1.262, catalyzes a dehydrogenation, generating an unstable product that undergoes a decarboxylation to yield the central metabolite DIHYDROXY-ACETONE-PHOSPHATE . 
#> 349  PROPANE-1-2-DIOL (propylene glycol) is produced from LACTALD during the bacterial degradation of L-rhamnose and L-Fucopyranoses (see PWY0-1315 and the associated pathway links) (in ). TAX-28901 can utilize PROPANE-1-2-DIOL as a carbon source and its metabolism may be a virulence factor . Under fermentative conditions, CPD-665 (propionaldehyde) is formed and converted to PROPANOL and PROPIONATE (propionate), ATP is produced, and NAD is regenerated (in ). The reactions forming CPD-665 and PROPIONYL-COA occur within a specialized bacterial microcompartment known as a metabolosome (see the About This Pathway section below).   The first enzyme in this pathway is CPLX1R65-171, a diol dehydratase composed of polypeptides PduCDE. This enzyme is dependent upon ADENOSYLCOBALAMIN (vitamin B12), which requires anaerobic conditions for its de novo biosynthesis (see pathway PWY-5507). Initially it was thought that PROPANE-1-2-DIOL utilization required oxygen, leading to an apparent ADENOSYLCOBALAMIN paradox. This paradox was resolved by the discovery that TAX-28901 can use CPD-14 as an alternative electron acceptor during anaerobic degradation of both PROPANE-1-2-DIOL and ETHANOL-AMINE .  The adenosyl group of the ADENOSYLCOBALAMIN cofactor of CPLX1R65-171 is unstable in vivo and can be converted to CPD0-1256 in a reaction that inactivates the enzyme . Active enzyme can be regenerated within the microcompartment by the combined action of a reactivase PduGH  (biochemically characterized in TAX-571 ), a dedicated STM2050-MONOMER , and a cobalamin reductase PduS, that together regenerate ADENOSYLCOBALAMIN. Thus, the microcompartment contains a complete ADENOSYLCOBALAMIN recycling system  (and in ).  Under aerobic conditions, exogenously added ADENOSYLCOBALAMIN (or its precursor COBINAMIDE) can be imported into the cell to allow use of PROPANE-1-2-DIOL as a carbon and energy source via the PWY0-42 as shown in the pathway link  and in .  There is also evidence for the production of PROPANOL and PROPIONATE via this pathway in TAX-1637 species, although the PWY0-42 may not be present in these organisms .  In addition to the organisms shown here, genes for ADENOSYLCOBALAMIN-dependent propanediol utilization have also been identified in species of FRAME:TAX-620, FRAME:TAX-629, FRAME:TAX-1578 and FRAME:TAX-1357, as well as FRAME:TAX-562 E24377A (reviewed in ).  About This Pathway  The first two reactions of this pathway occur within the lumen of a bacterial microcompartment (BMC), a capsid-like shell composed of multiple proteins that is analogous to the carboxysome of cyanobacteria. The ADENOSYLCOBALAMIN-dependent PROPANE-1-2-DIOL utilization (Pdu) microcompartment consists of major polypeptides PduABB'CDEGHJKOPTU which includes enzymes PduCDE, PduGH, PduO and PduP, as well as BMC-domain containing proteins PduABB'JKNTU . A structural proein (PduM) has also been described . The proposed function of this microcompartment is to sequester toxic CPD-665 (propionaldehyde).  The remaining reactions occur in the cytosol (reviewed in  and ).  The conversion of CPD-665 to PROPANOL is catalyzed within the compartment by an alcohol dehydrogenasealdehyde reductase encoded by the STM2052 gene, which provides ample supply of NAD+ to EC-1.2.1.87 (STM2051 PduP) .  The PROPIONATE formed can be activated to PROPIONYL-COA by the reversible reactions catalyzed by PduW and PduL. PROPIONYL-COA Propanoyl-CoA can also be formed from PROPIONATE by the Acs and PrpE enzymes in pathway PWY0-42 (see EC-6.2.1.17). PROPIONYL-COA Propanoyl-CoA is thus an intermediate in the catabolism of both PROPANE-1-2-DIOL (this pathway) and PROPIONATE (in pathway PWY0-42) and is necessary for expression of the prpBCDE genes of the PWY0-42 during growth on PROPANE-1-2-DIOL .  The genes of this pathway are clustered and are located adjacent to the cbicob gene cluster needed for ADENOSYLCOBALAMIN biosynthesis. Both clusters are regulated by the pocR gene ( and in ).  The 21 gene regulon encoding the pdu organelle and propanediol degradation enzymes from TAX-546 was cloned and expressed successfully in TAX-562, allowing PROPANE-1-2-DIOL utilization . 
#> 748                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   Bacteria do not possess the PWY3O-450 eukaryote-like CDP-choline-dependent pathway for phosphatidylcholine synthesis, but it is now well-documented that a wide range of Gram-positive and Gram-negative species that colonize mucosal tissue utilize a pathway very similar to the eukaryotic pathway as a means to decorate cell wall glycoconjugates with phosphocholine. While the addition of phosphocholine to these glycoconjugates is not as common as acetylation or methylation, it has far more potent biological consequences due to the multiple interactions between phosphocholine and host proteins, which are important at all stages of the infection process by pathogens .   These phosphocholine glycoconjugate modifications are found on various surface structures depending on the species. Examples include Lipopolysaccharides lipopolysaccharides (TAX-727) , pili (TAX-487 and TAX-485) , Lipoteichoic-Acids lipoteichoic acid (TAX-1313) , a 43-kDa surface protein (TAX-287) , and surface components of diverse commensal oral species .  The first two enzyme activities in the pathway, EC-2.7.1.32, and EC-2.7.7.15, are essentially identical to those of the eukaryotic pathway. These enzymes are encoded by homologues of the eukaryotic genes (G-44724 and G-44725, respectively) that have been detected in many species that colonize mucosal tissue.  The third and last enzyme, diphosphonucleoside choline phosphotransferase, catalyzes the transfer of the phosphoclholine moiety to a cell surface glycoconjugate. The gene encoding this enzyme has not been definitively identified in many of these organisms , but is proposed to be G-43708 in TAX-727  and TAX-1313 . 
#> 749                                                                                                                           C-type cytochromes (cytC) are electron-transfer proteins that have one or more HEME_C groups. They are characterized by the covalent attachment of the heme to the polypeptide chain via two (or rarely one) thioether bonds formed between thiol groups of cysteine residues and vinyl groups of a PROTOHEME molecule. The two cysteine residues almost always occur in the amino-acid sequence CXXCH .  C-type cytochromes possess a wide range of properties and function in a large number of different redox processes. In general, the differences in c-type cytochromes among bacteria are much larger than those among animal species. The use of the term term c-type cytochromes was introduced to distinguish this family of diverse proteins from the well-characterized mitochondrial protein, which is often referred to as cytochrome c.  Bacterial c-type cytochromes function in the electron transport chains of bacteria with many different types of energy metabolism, including phototrophs, methylotrophs, denitrifiers, sulfate reducers and the nitrogen-fixers.   The biogenesis of c-type cytochromes involves covalent attachment of PROTOHEME to two cysteines at a conserved CXXCH sequence in the apocytochrome. The histidine residue acts as an axial ligand to the heme iron. Three unique cytochrome c assembly pathways, known as systems I, II, and III, have ben described. This pathway describes system III, which is found in the mitochondria of fungi, vertebrates, invertebrates, and sme protozoa.  About This Pathway  The first indication for the existence of the system II pathway was the discovery that the ccsA gene (named for cytochrome c synthesis) in chloroplasts contains a WWD domain (a tryptophan-rich domain that is known to maintain a heme molecule in the reduced state by two external histidines) that is different from that found in the bacterial CcmF genes of system I . The gene was later shown to be required for cytochrome f biosynthesis (cytochrome f, a part of the cytochrome b6f complex, is analogous to cytochrome c1 of the bc1 complex and contains the CXXCH motif), and led to the discovery of the rest of the involved genes . Analogous genes were discovered in the Gram-positive bacterium TAX-1423  and the &beta;-proteobacterium TAX-520 .  The system is now known to operate in a wide range of organisms, including chloroplasts, Gram-positive bacteria, cyanobacteria, most &beta;-proteobacteria, some &delta;-proteobacteria, and all &epsilon;-proteobacteria .  In organisms that utilize system II the apocytochrome polypeptide, upon translocation from the cytosol, is exposed to the activity of EC-1.8.4.15, which forms a disulfide bridge between the two cysteines in the CXXCH motif. This bridge must be reduced before the heme molecule can be attached to the polypeptide. System II uses reducing equivalents from inside the cell (likely thioredoxin) to reduce the cysteine residues. Dedicated thioredoxin-related proteins called ResA (in TAX-1423) and CcsX (in TAX-520) can reduce the cysteines directly . These proteins are kept in their reduced active form by the DsbD protein, or by the similar, smaller protein CcdA .  The other components of system II are the CcsA and CcsB proteins, which form a tight complex with the activity of EC-4.4.1.17. In a few bacterial species (such as TAX-209, TAX-816, and TAX-843) the two proteins are fused into one large ORF called CcsBA. The two proteins are membrane-bound. CcsA, which has either eight or six transmembrane domains , contains the WWD domain that binds the heme molecule, as well as conserved histidine residues that are believed to help keep the heme iron in its reduced form. The reduced heme molecule moves from the cytosol through CcsA until it reaches the WWD domain, where it is liganded in place. CcsB has four transmembrane domains, with the third and fourth separated by a large periplasmic domain that includes the apocytochrome c CXXCH recognition function. CcsB thus recruites the apocytochrome polypeptide and presents it to the heme bound to CcsA . 
#> 799                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     The cyclization is not limited to free amino acids. N-terminal glutaminyl and glutamyl residues of polypeptides can also spontaneously cyclize to 5-OXOPROLINE residues. The former transformation is also enzymically catalyzed by EC-2.3.2.5 . Once formed, the abberant residue is removed from the peptide by EC-3.4.19.3, yielding free 5-OXOPROLINE .  5-OXOPROLINE is also formed by additional routes. One such route is from GLUTATHIONE, in a reaction catalyzed by EC-4.3.2.7 . In eukaryotic organisms it is also formed by EC-4.3.2.9 as part of the PWY-4041 .   Accumulation of 5-OXOPROLINE is deleterious to the organism. Its reported effects include growth inhibition in prokaryotes and plants  and interference with energy production, lipid synthesis, and antioxidant defenses in mammalian brain . Human inborn errors of glutathione metabolism that lead to 5-OXOPROLINE buildup result in metabolic acidosis, hemolytic anemia, neurological problems, and massive urinary excretion of 5-OXOPROLINE .  An ATP-dependent enzyme that hydrolyzes 5-OXOPROLINE to GLT, EC-3.5.2.9, has been known for a while, as it functions in the PWY-4041. Until recently, the enzyme, which is encoded by the HS11319 gene, was thought to be a eukaryotic enzyme exclusively . However, in 2017 a widespread prokaryotic enzyme that catalyzes the same reaction was discovered . The trimeric enzyme is encoded by the G6382, G6380, and G6381 genes. Inactivation of any of these genes in TAX-1423 slowed growth and caused 5-OXOPROLINE accumulation in the cells. 
#>                                                         commonName
#> 12                                   tetrachloroethene degradation
#> 173                                    D-erythronate degradation I
#> 349                        (<i>S</i>)-propane-1,2-diol degradation
#> 748 cell-surface glycoconjugate-linked phosphocholine biosynthesis
#> 749                cytochrome <i>c</i> biogenesis (system II type)
#> 799                                     5-oxo-L-proline metabolism
#>        species taxonomicRange
#> 12   TAX-55583          TAX-2
#> 173 TAX-498761          TAX-2
#> 349   TAX-1642          TAX-2
#> 748   TAX-1313          TAX-2
#> 749 TAX-224308      TAX-28216
#> 799   TAX-1423       TAX-2759
#>                                                                                                   spConverted
#> 12     Bacteria;Bacillota;Clostridia;Eubacteriales;Desulfitobacteriaceae;Dehalobacter;Dehalobacter restrictus
#> 173 Bacteria;Bacillota;Clostridia;Eubacteriales;Heliobacteriaceae;Heliomicrobium;Heliomicrobium modesticaldum
#> 349                              Bacteria;Bacillota;Bacilli;Bacillales;Listeriaceae;Listeria;Listeria innocua
#> 748        Bacteria;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;Streptococcus pneumoniae
#> 749                              Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus subtilis
#> 799                              Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus subtilis
```

This time we use a random dendrogram and rename the rows with taxonomic name included in MetaCyc.


```r

while (TRUE) {
  uniq <- sample(bc$species, 10)
  if (length(uniq)==unique(length(uniq))) {
    if (sum(is.na(sapply(strsplit(bc[bc$species %in% uniq,]$spConverted, ";"), "[", 7)))==0) {
      break
    }
  }
}
data <- matrix(sample(seq(1,100),100), ncol = 10)
rownames(data) <- uniq

## Dendrogram
dhc <- as.dendrogram(hclust(dist(data)))
plot(dhc)
```

<img src="02-microbiome_usage_files/figure-html/example-1.png" width="672" />
Using the dendrogram (`dhc` argument) and input data.frame with the column name `query`, we can plot the dendrogram with pathway information. Note we need to provide the function named vector of nodes (corresponding to gene clusters when the gene as input). This time, bi-gram wordclouds are to be plotted, by specifying in `argList`.


```r
# Some filtering
load("../allFreqMetaCyc.rda")
deleter <- (allFreqMetaCyc |> dplyr::filter(freq>5000) |> dplyr::select(word))$word
# Named vector
sampled <- row.names(data)
names(sampled) <- sampled
# Query column (including node names) is specified as query
input$query <- input$species

micro <- plotEigengeneNetworksWithWords(NA, sampled,
                               useWC = TRUE, # Use wordcloud
                               useFunc = "wcMan", # Use manual function (as the input is custom data.frame)
                               useDf=input,dendPlot="ggplot",dhc=dhc,
                               argList=list(additionalRemove=deleter,
                                ngram=1),
                               useWGCNA=FALSE, spacer=1,
                               horiz=FALSE, wcScale = 10)
micro + scale_y_continuous(expand=c(0,20))
```

<img src="02-microbiome_usage_files/figure-html/plotDendroWord-1.png" width="1440" />


```r
sessionInfo()
#> R version 4.2.2 (2022-10-31 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 22621)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=Japanese_Japan.utf8 
#> [2] LC_CTYPE=Japanese_Japan.utf8   
#> [3] LC_MONETARY=Japanese_Japan.utf8
#> [4] LC_NUMERIC=C                   
#> [5] LC_TIME=Japanese_Japan.utf8    
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets 
#> [6] methods   base     
#> 
#> other attached packages:
#> [1] RColorBrewer_1.1-3   ggraph_2.1.0        
#> [3] ggplot2_3.4.0        wcGeneSummary_0.99.0
#> 
#> loaded via a namespace (and not attached):
#>   [1] GeneSummary_0.99.4     colorspace_2.0-3      
#>   [3] rjson_0.2.21           ellipsis_0.3.2        
#>   [5] XVector_0.38.0         GlobalOptions_0.1.2   
#>   [7] base64enc_0.1-3        ggdendro_0.1.23       
#>   [9] fs_1.5.2               rstudioapi_0.14       
#>  [11] farver_2.1.1           graphlayouts_0.8.4    
#>  [13] ggrepel_0.9.2          bit64_4.0.5           
#>  [15] AnnotationDbi_1.60.0   fansi_1.0.3           
#>  [17] xml2_1.3.3             downlit_0.4.2         
#>  [19] cachem_1.0.6           knitr_1.41            
#>  [21] polyclip_1.10-4        jsonlite_1.8.4        
#>  [23] png_0.1-8              graph_1.76.0          
#>  [25] ggforce_0.4.1          shiny_1.7.4           
#>  [27] bugsigdbr_1.5.3        rentrez_1.2.3         
#>  [29] compiler_4.2.2         httr_1.4.4            
#>  [31] assertthat_0.2.1       fastmap_1.1.0         
#>  [33] cli_3.6.0              later_1.3.0           
#>  [35] tweenr_2.0.2           htmltools_0.5.4       
#>  [37] tools_4.2.2            igraph_1.3.5          
#>  [39] NLP_0.2-1              gtable_0.3.1          
#>  [41] glue_1.6.2             GenomeInfoDbData_1.2.9
#>  [43] dplyr_1.0.10           Rcpp_1.0.9            
#>  [45] slam_0.1-50            Biobase_2.58.0        
#>  [47] jquerylib_0.1.4        vctrs_0.5.1           
#>  [49] Biostrings_2.66.0      xfun_0.36             
#>  [51] stringr_1.5.0          mime_0.12             
#>  [53] lifecycle_1.0.3        pvclust_2.2-0         
#>  [55] XML_3.99-0.13          dendextend_1.16.0     
#>  [57] org.Hs.eg.db_3.16.0    zlibbioc_1.44.0       
#>  [59] MASS_7.3-58.1          scales_1.2.1          
#>  [61] tidygraph_1.2.2        promises_1.2.0.1      
#>  [63] parallel_4.2.2         cyjShiny_1.0.34       
#>  [65] yaml_2.3.6             memoise_2.0.1         
#>  [67] gridExtra_2.3          yulab.utils_0.0.6     
#>  [69] sass_0.4.4             stringi_1.7.12        
#>  [71] RSQLite_2.2.20         highr_0.10            
#>  [73] S4Vectors_0.36.1       BiocGenerics_0.44.0   
#>  [75] GenomeInfoDb_1.34.6    rlang_1.0.6           
#>  [77] pkgconfig_2.0.3        bitops_1.0-7          
#>  [79] evaluate_0.19          purrr_1.0.1           
#>  [81] patchwork_1.1.2        htmlwidgets_1.6.1     
#>  [83] cowplot_1.1.1          bit_4.0.5             
#>  [85] tidyselect_1.2.0       magrittr_2.0.3        
#>  [87] bookdown_0.31          R6_2.5.1              
#>  [89] IRanges_2.32.0         generics_0.1.3        
#>  [91] DBI_1.1.3              pillar_1.8.1          
#>  [93] withr_2.5.0            KEGGREST_1.38.0       
#>  [95] RCurl_1.98-1.9         tibble_3.1.8          
#>  [97] crayon_1.5.2           wordcloud_2.6         
#>  [99] utf8_1.2.2             rmarkdown_2.19        
#> [101] viridis_0.6.2          GetoptLong_1.0.5      
#> [103] grid_4.2.2             blob_1.2.3            
#> [105] digest_0.6.31          xtable_1.8-4          
#> [107] tm_0.7-10              tidyr_1.2.1           
#> [109] httpuv_1.6.8           gridGraphics_0.5-1    
#> [111] stats4_4.2.2           munsell_0.5.0         
#> [113] ggplotify_0.1.0        viridisLite_0.4.1     
#> [115] bslib_0.4.2
```
