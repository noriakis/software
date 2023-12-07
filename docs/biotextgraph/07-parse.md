

# Parsing

Some functions can be used to the parsing purpose.

## Parsing Enzyme Comission numbers, its description and comment.


```r
library(biotextgraph)
#>  要求されたパッケージ ggplot2 をロード中です
#> 
#> Registered S3 method overwritten by 'pvclust':
#>   method       from      
#>   text.pvclust dendextend
library(data.table)
ecdf <- enzyme("enzyme.dat",ecnum="all", onlyDf=TRUE)
#> Processing EC file
data.table(ecdf) |> head()
#>     number                                 desc
#> 1: 1.1.1.1                alcohol dehydrogenase
#> 2: 1.1.1.2      alcohol dehydrogenase (NADP(+))
#> 3: 1.1.1.3             homoserine dehydrogenase
#> 4: 1.1.1.4       (R,R)-butanediol dehydrogenase
#> 5: 1.1.1.5 Transferred entry: 111303 and 111304
#> 6: 1.1.1.6               glycerol dehydrogenase
#>                                                                                                                                                                                                                                           comment
#> 1: Acts on primary or secondary alcohols or hemi-acetals with very broad     specificity; however the enzyme oxidizes methanol much more poorly     than ethanol The animal, but not the yeast, enzyme acts also on cyclic secondary     alcohols
#> 2:                                                    Some members of this group oxidize only primary alcohols; others act     also on secondary alcohols May be identical with EC 11119, EC 11133 and EC 11155 Re-specific with respect to NADPH
#> 3:                 The enzyme from Saccharomyces cerevisiae acts most rapidly with     NAD(+); the Neurospora enzyme with NADP(+) The enzyme from Escherichia coli is a multifunctional protein, which     also catalyzes the reaction of EC 2724
#> 4:                                                                                                                                                                                     Also converts diacetyl into acetoin with NADH as reductant
#> 5:                                                                                                                                                                                                                                               
#> 6:                                                                                                                                                                                                                   Also acts on 1,2-propanediol
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     DRs
#> 1: P07327,ADH1A_HUMAN;P28469,ADH1A_MACMU;Q5RBP7,ADH1A_PONAB;P25405,ADH1A_SAAHA;P25406,ADH1B_SAAHA;P00327,ADH1E_HORSE;P00326,ADH1G_HUMAN;O97959,ADH1G_PAPHA;P00328,ADH1S_HORSE;P80222,ADH1_ALLMI;P30350,ADH1_ANAPL;P49645,ADH1_APTAU;P06525,ADH1_ARATH;P41747,ADH1_ASPFN;Q17334,ADH1_CAEEL;P43067,ADH1_CANAX;P85440,ADH1_CATRO;P14219,ADH1_CENAM;P48814,ADH1_CERCA;Q70UN9,ADH1_CERCO;P23991,ADH1_CHICK;P86883,ADH1_COLLI;P19631,ADH1_COTJA;P23236,ADH1_OHY;P48586,ADH1_OMN;P09370,ADH1_OMO;P22246,ADH1_OMT;P07161,ADH1_OMU;P12854,ADH1_ONA;P08843,ADH1_EMENI;A0A165U5V5,ADH1_EUPLT;P26325,ADH1_GADMC;Q9Z2M2,ADH1_GEOAT;Q64413,ADH1_GEOBU;Q64415,ADH1_GEOKN;P12311,ADH1_GEOSE;P05336,ADH1_HORVU;P20369,ADH1_KLULA;Q07288,ADH1_KLUMA;P00333,ADH1_MAIZE;P86885,ADH1_MESAU;P00329,ADH1_MOUSE;P80512,ADH1_NAJNA;Q9P6C8,ADH1_NEUCR;Q75ZX4,ADH1_ORYSI;Q2R8Z5,ADH1_ORYSJ;P12886,ADH1_PEA;P22797,ADH1_PELPE;P41680,ADH1_PERMA;P25141,ADH1_PETHY;O00097,ADH1_PICST;Q03505,ADH1_RABIT;P06757,ADH1_RAT;P14673,ADH1_SOLTU;P80338,ADH1_STRCA;P13603,ADH1_TRIRP;P00330,ADH1_YEAST;Q07264,ADH1_ZEALU;P20368,ADH1_ZYMMO;O45687,ADH2_CAEEL;O94038,ADH2_CANAL;P48815,ADH2_CERCA;Q70UP5,ADH2_CERCO;Q70UP6,ADH2_CERRO;P27581,ADH2_OAR;P25720,ADH2_OBU;P23237,ADH2_OHY;P48587,ADH2_OMN;P09369,ADH2_OMO;P07160,ADH2_OMU;P24267,ADH2_OWH;P37686,ADH2_ECOLI;P54202,ADH2_EMENI;Q24803,ADH2_ENTHI;P42327,ADH2_GEOSE;P10847,ADH2_HORVU;P49383,ADH2_KLULA;Q9P4C2,ADH2_KLUMA;P04707,ADH2_MAIZE;Q4R1E8,ADH2_ORYSI;Q0ITW7,ADH2_ORYSJ;O13309,ADH2_PICST;P28032,ADH2_SOLLC;P14674,ADH2_SOLTU;F2Z678,ADH2_YARLI;P00331,ADH2_YEAST;F8DVL8,ADH2_ZYMMA;P0DJA2,ADH2_ZYMMO;P07754,ADH3_EMENI;P42328,ADH3_GEOSE;P10848,ADH3_HORVU;P49384,ADH3_KLULA;P14675,ADH3_SOLTU;P07246,ADH3_YEAST;P49385,ADH4_KLULA;Q09669,ADH4_SCHPO;A6ZTT5,ADH4_YEAS7;P10127,ADH4_YEAST;Q6XQ67,ADH5_SACPS;P38113,ADH5_YEAST;P28332,ADH6_HUMAN;P41681,ADH6_PERMA;Q5R7Z8,ADH6_PONAB;Q5XI95,ADH6_RAT;P40394,ADH7_HUMAN;Q64437,ADH7_MOUSE;P41682,ADH7_RAT;P9WQC0,ADHA_MYCTO;P9WQC1,ADHA_MYCTU;O31186,ADHA_RHIME;Q7U1B9,ADHB_MYCBO;P9WQC6,ADHB_MYCTO;P9WQC7,ADHB_MYCTU;P9WQB8,ADHD_MYCTO;P9WQB9,ADHD_MYCTU;P33744,ADHE_CLOAB;P0A9Q8,ADHE_ECO57;P0A9Q7,ADHE_ECOLI;P81600,ADHH_GADMO;P72324,ADHI_CERS4;Q9SK86,ADHL1_ARATH;Q9SK87,ADHL2_ARATH;A1L4Y2,ADHL3_ARATH;Q8VZ49,ADHL4_ARATH;Q0V7W6,ADHL5_ARATH;Q8LEB2,ADHL6_ARATH;Q9FH04,ADHL7_ARATH;P81601,ADHL_GADMO;P39451,ADHP_ECOLI;O46649,ADHP_RABIT;O46650,ADHQ_RABIT;Q96533,ADHX_ARATH;Q3ZC42,ADHX_BOVIN;Q17335,ADHX_CAEEL;Q54TC2,ADHX_DICDI;P46415,ADHX_OME;P19854,ADHX_HORSE;P11766,ADHX_HUMAN;P93629,ADHX_MAIZE;P28474,ADHX_MOUSE;P80360,ADHX_MYXGL;P81431,ADHX_OCTVU;A2XAZ3,ADHX_ORYSI;Q0DWH1,ADHX_ORYSJ;P80572,ADHX_PEA;O19053,ADHX_RABIT;P12711,ADHX_RAT;P80467,ADHX_SAAHA;P86884,ADHX_SCYCA;P79896,ADHX_SPAAU;Q9NAR7,ADH_BACOL;P14940,ADH_CUPNE;Q0KDL6,ADH_CUPNH;Q00669,ADH_OAD;P21518,ADH_OAF;P25139,ADH_OAM;Q50L96,ADH_OAN;P48584,ADH_OBO;P22245,ADH_ODI;Q9NG42,ADH_OEQ;P28483,ADH_OER;P48585,ADH_OFL;P51551,ADH_OGR;Q09009,ADH_OGU;P51549,ADH_OHA;P21898,ADH_OHE;Q07588,ADH_OIM;Q9NG40,ADH_OIN;Q27404,ADH_OLA;P10807,ADH_OLE;P07162,ADH_OMA;Q09010,ADH_OMD;P00334,ADH_OME;Q00671,ADH_OMM;P25721,ADH_OMY;Q00672,ADH_ONI;P07159,ADH_OOR;P84328,ADH_OPB;P37473,ADH_OPE;P23361,ADH_OPI;P23277,ADH_OPL;Q6LCE4,ADH_OPS;Q9U8S9,ADH_OPU;Q9GN94,ADH_OSE;Q24641,ADH_OSI;P23278,ADH_OSL;Q03384,ADH_OSU;P28484,ADH_OTE;P51550,ADH_OTS;B4M8Y0,ADH_OVI;Q05114,ADH_OWI;P26719,ADH_OYA;P17648,ADH_FRAAN;P48977,ADH_MALDO;Q8GIX7,ADH_MORSE;P9WQC2,ADH_MYCTO;P9WQC3,ADH_MYCTU;P39462,ADH_SACS2;P25988,ADH_SCAAL;Q00670,ADH_SCACA;P00332,ADH_SCHPO;Q2FJ31,ADH_STAA3;Q2G0G1,ADH_STAA8;Q2YSX0,ADH_STAAB;Q5HI63,ADH_STAAC;Q99W07,ADH_STAAM;Q7A742,ADH_STAAN;Q6GJ63,ADH_STAAR;Q6GBM4,ADH_STAAS;Q8NXU1,ADH_STAAW;Q5HRD6,ADH_STAEQ;Q8CQ56,ADH_STAES;Q4J781,ADH_SULAC;P50381,ADH_SULSR;Q96XE0,ADH_SULTO;A0A0S1X9S7,ADH_THEBA;P51552,ADH_ZAPTU;Q5AR48,ASQE_EMENI;A5JYX5,DHS3_CAEEL;P76553,EUTG_ECOLI;P41795,EUTG_SALTY;P32771,FADH_YEAST;A7ZIA4,FRMA_ECO24;Q8X5J4,FRMA_ECO57;A7ZX04,FRMA_ECOHS;A1A835,FRMA_ECOK1;Q0TKS7,FRMA_ECOL5;Q8FKG1,FRMA_ECOL6;B1J085,FRMA_ECOLC;P25437,FRMA_ECOLI;B1LIP1,FRMA_ECOSM;Q1RFI7,FRMA_ECOUT;P44557,FRMA_HAEIN;P39450,FRMA_PHODP;Q3Z550,FRMA_SHISS;P73138,FRMA_SYNY3;E1ACQ9,NOTN_ASPSM;N4WE73,OXI1_COCH4;N4WE43,RED2_COCH4;N4WW42,RED3_COCH4;P33010,TERPD_PSESP;O07737,Y1895_MYCTU
#> 2:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     Q6AZW2,A1A1A_DANRE;Q568L5,A1A1B_DANRE;Q24857,ADH3_ENTHI;Q04894,ADH6_YEAST;P25377,ADH7_YEAST;O57380,ADH8_PELPE;P74721,ADHA_SYNY3;Q9F282,ADHA_THEET;P0CH36,ADHC1_MYCS2;P0CH37,ADHC2_MYCS2;P0A4X1,ADHC_MYCBO;P9WQC4,ADHC_MYCTO;P9WQC5,ADHC_MYCTU;P27250,AHR_ECOLI;Q3ZCJ2,AK1A1_BOVIN;Q5ZK84,AK1A1_CHICK;O70473,AK1A1_CRIGR;P14550,AK1A1_HUMAN;Q9JII6,AK1A1_MOUSE;P50578,AK1A1_PIG;Q5R5D5,AK1A1_PONAB;P51635,AK1A1_RAT;Q6GMC7,AK1A1_XENLA;Q28FD1,AK1A1_XENTR;Q9UUN9,ALD2_SPOSA;P27800,ALDX_SPOSA;P75691,YAHK_ECOLI;Q46856,YQHD_ECOLI
#> 3:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               P00561,AK1H_ECOLI;P27725,AK1H_SERMA;P00562,AK2H_ECOLI;Q9SA18,AKH1_ARATH;P49079,AKH1_MAIZE;O81852,AKH2_ARATH;P49080,AKH2_MAIZE;P57290,AKH_BUCAI;Q8K9U9,AKH_BUCAP;Q89AR4,AKH_BUCBP;P37142,AKH_DAUCA;P44505,AKH_HAEIN;P19582,DHOM_BACSU;P08499,DHOM_CORGL;Q5B998,DHOM_EMENI;Q9ZL20,DHOM_HELPJ;P56429,DHOM_HELPY;Q9CGD8,DHOM_LACLA;P52985,DHOM_LACLC;P37143,DHOM_METGL;Q58997,DHOM_METJA;P63630,DHOM_MYCBO;P46806,DHOM_MYCLE;P9WPX0,DHOM_MYCTO;P9WPX1,DHOM_MYCTU;P29365,DHOM_PSEAE;O94671,DHOM_SCHPO;P52986,DHOM_SYNY3;P31116,DHOM_YEAST;P37144,DHON_METGL
#> 4:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              P14940,ADH_CUPNE;Q0KDL6,ADH_CUPNH;P39714,BDH1_YEAST;O34788,BDHA_BACSU;Q00796,DHSO_HUMAN
#> 5:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
#> 6:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    A4IP64,ADH1_GEOTN;O13702,GLD1_SCHPO;P45511,GLDA_CITFR;P0A9S6,GLDA_ECOL6;P0A9S5,GLDA_ECOLI;P32816,GLDA_GEOSE;P50173,GLDA_PSEPU;Q9WYQ4,GLDA_THEMA;Q92EU6,GOLD_LISIN
```

## Parsing UniProt organism identifiers.


```r
updf <- getUPtax("speclist.txt",candUP="all")
data.table(updf) |> head()
#>     UPID                                           Taxonomy
#> 1: AADNV Aedes albopictus densovirus (isolate Boublik/1994)
#> 2:  AAV2                           Adeno-associated virus 2
#> 3: AAV2S Adeno-associated virus 2 (isolate Srivastava/1982)
#> 4: ABAMA                                     Abacion magnum
#> 5: ABANI                                     Abaeis nicippe
#> 6: ABAPA                              Abax parallelepipedus
```


## Parsing MetaCyc pathway text information


```r
candidateSpecies <- c("Escherichia coli")
file <- "metacyc/24.5/data/pathways.dat"
metacyc <- parseMetaCycPathway(file, candidateSpecies, withTax=TRUE, clear=TRUE)
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Throughout evolution, bioluminescence
#> has been reinvented many times; some 30 different
#> independent systems are still extant. The responsible
#> enzymes are unrelated in bacteria, unicellular algae,
#> coelenterates, beetles, fishes, and other organisms, al...'
#> to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Throughout evolution, bioluminescence
#> has been reinvented many times; some 30 different
#> independent systems are still extant. The responsible
#> enzymes are unrelated in bacteria, unicellular algae,
#> coelenterates, beetles, fishes, and other organisms, al...'
#> to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'However, only a handful of
#> biologically synthesized carbon-fluorine bonds are known
#> |CITS: [15209126]|. The most common natural organofluorine
#> species is |FRAME: CPD-12709|, a toxin found in the leaves
#> and seeds of a variety of tropical plants, often i...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Nitrification, the oxidation of
#> |FRAME: AMMONIA| to |FRAME: NITRATE| by microorganisms, is
#> a key process in the global nitrogen cycle, resulting in
#> nitrogen loss from ecosystems, eutrophication of surface
#> and ground waters, and the production of atmos...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'The most heavily modified class of RNA
#> molecules is the tRNA - an average of 11%, 14%, and 17% of
#> the tRNA residues are modified in |FRAME: TAX-562|, |FRAME:
#> TAX-4932|, and mammalian cytoplasm, respectively |CITS:
#> [9399820]|. These modifications expand...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' 2-Thiouridine derivatives are present
#> at position 34 (the <91>wobble<92> base) in anticodons of
#> tRNA for |FRAME: LYS|, |FRAME: GLT|, and |FRAME: GLN| in
#> all domains of life (see |FRAME: PWY-7892| for an example),
#> but also at position 54 in the T-loop ...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'The most heavily modified class of RNA
#> molecules is the tRNA - an average of 11%, 14%, and 17% of
#> the tRNA residues are modified in |FRAME: TAX-562|, |FRAME:
#> TAX-4932|, and mammalian cytoplasm, respectively |CITS:
#> [9399820]|. These modifications expand...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' In many bacteria the enzyme is
#> encoded by the |FRAME: EG11528| gene. Triclosan forms a
#> tight, non-covalent ternary complex with the oxidized NADH
#> cofactor in the active site of FabI, thus inactivating the
#> enzyme |CITS: [10196195]|.  However, some bact...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Erythromycin A is used in clinical
#> medicine against infections caused by Gram-positive
#> bacteria. It is also used for many pulmonary infections
#> such as Legionnaire<92>s disease and as an alternative for
#> patients allergic to penicillins |CITS: [11851474...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' The predominant circulating
#> pyrimidine in many species, including humans is |FRAME:
#> URIDINE| |CITS: [8155674][16769123]|. In mammals |FRAME:
#> URIDINE Uridine| is produced <i>de novo<i> in the liver and
#> kidney, and is circulated to other organs, where i...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: CPD0-1308 Glyphosate| is a
#> broad-spectrum systemic herbicide used to kill weeds. It
#> was initially patented and sold by the Monsanto Company in
#> the 1970s under the tradename Roundup, and has become the
#> most used herbicide in the world |CITS: [2...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Queuosine and its derivatives occur
#> exclusively at position 34 (the wobble position) in the
#> anticodons of tRNAs coding for the amino acids |FRAME:
#> HIS|, |FRAME:L-ASPARTATE|, |FRAME: ASN| and |FRAME: TYR|
#> |CITS: [ 4550561]|. Each of these tRNAs posses ...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' The processing of tRNA in |FRAME:
#> TAX-2157| and |FRAME: TAX-2759| involves the splicing of
#> introns by a mechanism that involves protein enzymes
#> (unlike mRNA splicing, which is RNA catalyzed). While there
#> are some examples of tRNA intron splicing in pr...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' The processing of tRNA in |FRAME:
#> TAX-2157| and |FRAME: TAX-2759| involves the splicing of
#> introns by a mechanism that involves three protein enzymes
#> (unlike mRNA splicing, which is RNA catalyzed). While there
#> are some examples of tRNA intron splicing...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Erythromycin A is used in clinical
#> medicine against infections caused by Gram-positive
#> bacteria. It is also used for many pulmonary infections
#> such as Legionnaire<92>s disease and as an alternative for
#> patients allergic to penicillins |CITS: [11851474...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' The ABH and the related Lewis blood
#> group determinants are formed from |FRAME:
#> Blood-Group-Antigen-Precursor "precursor disaccharides"|
#> that decorate glycosphingolipids or protein glycans by
#> concerted actions of several glycosyltransferases.  The ABH
#> ...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'The glucosinolates play an important
#> role in plant defense. Upon herbivore attack glucosinolates
#> are hydrolyzed and broken down to toxic compounds in a
#> process that starts with the action of |FRAME: EC-3.2.1.147
#> "EC 3.2.1.147, myrosinase|. The myrosina...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Erythromycin A is used in clinical
#> medicine against infections caused by Gram-positive
#> bacteria. It is also used for many pulmonary infections
#> such as Legionnaire<92>s disease and as an alternative for
#> patients allergic to penicillins |CITS: [11851474...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Erythromycin A is used in clinical
#> medicine against infections caused by Gram-positive
#> bacteria. It is also used for many pulmonary infections
#> such as Legionnaire<92>s disease and as an alternative for
#> patients allergic to penicillins |CITS: [11851474...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'The most heavily modified class of RNA
#> molecules is the tRNA - an average of 11%, 14%, and 17% of
#> the tRNA residues are modified in |FRAME: TAX-562|, |FRAME:
#> TAX-4932|, and mammalian cytoplasm, respectively |CITS:
#> [9399820]|. These modifications expand...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' The term <91><91>indirect defense" is
#> used to describe a defense mechanism used by certain
#> plants, in which herbivore damage induces the emission of
#> volatile organic compounds that attract natural enemies of
#> the herbivores |CITS:[7753779][11251117]|. ...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Methylotrophs are organisms that are
#> capable of growing on C1 compounds, like methanol. These
#> organisms can derive all their energy and carbon needs from
#> reduced molecules that have no C<96>C bond. Methylotrophy
#> is found only in a few prokaryotic and ...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Deoxyribonucleotides are synthesised
#> <i>de novo<i> at the diphosphate level through reduction of
#> the 2<92>-hydroxyl group of the corresponding
#> ribonucleotides (see |FRAME: PWY-7184|). This reduction is
#> mediated by the key enzyme |FRAME: EC-1.17.4.1 "r...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Most pathogenic bacteria have an
#> absolute requirement for iron, which is scarce within their
#> host environment. In vertebrates the majority of the iron
#> is found in the form of heme. Several heme-iron acquisition
#> systems have been characterized from pat...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Most pathogenic bacteria have an
#> absolute requirement for iron, which is scarce within their
#> host environment. In vertebrates the majority of the iron
#> is found in the form of heme. Several heme-iron acquisition
#> systems have been characterized from pat...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate '|FRAME:URIDINE| to |FRAME: CPD-497|,
#> which is carried out by the enzyme pseudouridine synthase
#> in a reaction that does not require any cofactors. In
#> pseudouridine uracil is bound to the to ribose through C5
#> rather than through N1 as is the case for uri...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: D-APIOSE| is a unique
#> branched-chain pentose found principally in plants. It is a
#> key component of structurally complex cell wall
#> polysaccharides such as |FRAME: CPD-14581
#> "rhamnogalacturonan-II"| (RG-II) in the cell walls of
#> higher plants and...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: D-APIOSE| is a unique
#> branched-chain pentose found principally in plants. It is a
#> key component of structurally complex cell wall
#> polysaccharides such as |FRAME: CPD-14581
#> "rhamnogalacturonan-II"| (RG-II) in the cell walls of
#> higher plants and...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: D-APIOSE| is a unique
#> branched-chain pentose found principally in plants. It is a
#> key component of structurally complex cell wall
#> polysaccharides such as |FRAME: CPD-14581
#> "rhamnogalacturonan-II"| (RG-II) in the cell walls of
#> higher plants and...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: Glucuronoxylans
#> Glucuronoxylans| are a component of |FRAME:
#> Hemicelluloses|. They are linear polymers of
#> &beta;-D-xylopyranosyl units linked by (1&rarr;4)
#> glycosidic bonds, with many of the xylose units substituted
#> with 2, 3 or 2,3-linked gluc...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: ASCORBATE|, also known as
#> vitamin C, fulfils multiple essential roles in both plants
#> and animals. Being a strong reducing agent, it functions as
#> an antioxidant and a redox buffer. It is also a cofactor
#> for several enzymes, which are involved i...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: ASCORBATE|, also known as
#> vitamin C, fulfils multiple essential roles in both plants
#> and animals. Being a strong reducing agent, it functions as
#> an antioxidant and a redox buffer. It is also a cofactor
#> for several enzymes, which are involved i...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: ASCORBATE|, also known as
#> vitamin C, fulfils multiple essential roles in both plants
#> and animals. Being a strong reducing agent, it functions as
#> an antioxidant and a redox buffer. It is also a cofactor
#> for several enzymes, which are involved i...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: ASCORBATE|, also known as
#> vitamin C, fulfils multiple essential roles in both plants
#> and animals. Being a strong reducing agent, it functions as
#> an antioxidant and a redox buffer. It is also a cofactor
#> for several enzymes, which are involved i...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Oxalate degradation by ruminal
#> microbes was recognized as an important function quite
#> early |CITS: [Talapatra48][Morris55]|. Increased dietary
#> oxalate induces selection of oxalate-degrading ruminal
#> bacteria and make it possible for the animal to toler...'
#> to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Fungi are believed to utilize oxalate
#> in lignin degradation, nutrient availability, pathogenesis,
#> and competition. There is good correlation between
#> pathogenesis, virulence, and oxalate secretion. Solubility
#> of soil nutrients is achieved by soil-livin...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Deoxyribonucleotides are synthesised
#> <i>de novo<i> at the diphosphate level through reduction of
#> the 2<92>-hydroxyl group of the corresponding
#> ribonucleotides (see |FRAME: PWY-7184|). This reduction is
#> mediated by the key enzyme |FRAME: EC-1.17.4.1 "r...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'the diversity amongst glucosinolates
#> arises from the addition of different sized alkyl groups to
#> the side chain of the amino acids, principally valine,
#> phenylalanine and methionine.  Glucosinolates are found
#> prominently in the order |FRAME: TAX-3700|,...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Protein substrates of the |FRAME:
#> Ubiquitins ubiquitin| (Ub) system, which controls the
#> levels of many intracellular proteins, are conjugated to
#> ubiquitin through the action of |FRAME:
#> E3-independent-Ubiquitin-carrier-E2 "E2
#> ubiquitin-conjugating enzy...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Deoxyribonucleotides are synthesised
#> <i>de novo<i> at the diphosphate level through reduction of
#> the 2<92>-hydroxyl group of the corresponding
#> ribonucleotides (see |FRAME:PWY-7220| and
#> |FRAME:PWY-7222|). This reduction is mediated by the key
#> enzyme |F...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: E3-ubiquitin-carrier-proteins
#> "E3 ubiquitin transferases"| recognize their substrates by
#> specific degradation signals, known as degrons |CITS:
#> [11017125][14744430][15688063][15950869]|. Once recognized,
#> the protein substrate is ubiquitylated a...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Glycosylated microbial natural
#> products possess a wide range of pharmaceutical properties
#> such as anticancer, antiviral, antifungal or antibiotic.
#> The sugar moieties often contribute significantly to the
#> biological activities of the these compounds.  ...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'of sulfur. It is formed in megatonne
#> quantities in the atmosphere from the chemical oxidation of
#> atmospheric |FRAME: CPD-7670|. Depending on the
#> geographical latitude, 25<96>70% of the flux of
#> dimethylsulfide is oxidized to methanesulfonate, which is
#> d...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Manganese is the second most abundant
#> transition metal found in the Earth<92>s crust. In the
#> environment, manganese is mostly found in three different
#> oxidation states: II, III, and IV. Mn(II), primarily
#> occurring as the soluble |FRAME: MN+2| species,...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Manganese is the second most abundant
#> transition metal found in the Earth<92>s crust. In the
#> environment, manganese is mostly found in three different
#> oxidation states: II, III, and IV. Mn(II), primarily
#> occurring as the soluble |FRAME: MN+2| species,...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: TAX-920 | is an obligate
#> chemolithotroph that lives at extremely low pH, and is the
#> main bacterium used in the industrial extraction of copper
#> and uranium from ores, using the microbial leaching
#> technique. The sole energy-producing process tha...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: CPD-15423| is a prenylated
#> hydroxybenzoate unit that is attached to the 2-amino group
#> of the aminocoumarin core of the two aminocoumarin
#> antibiotics |FRAME: CPD-15417| and |FRAME: CPD-15421|,
#> which are produced by the bacterial species |FRAME:...' to
#> a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate '|FRAME: CPD-4886| (L-limonene) is
#> found in a variety of trees and herbs such as Mentha spp.,
#> while |FRAME: CPD-8785| (D-limonene) is the major component
#> of peel oil from oranges and lemons, and the essential oil
#> of caraway. The natural functions of lim...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' The mycobacterial cell envelope
#> contains an usually rich variety of complex free lipids and
#> glycolipids associated with the mycolic acid layer of the
#> cell wall. These lipids and glycolipids are located in the
#> membrane outer leaflet. Together with the ...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'There are two types of eumelanin,
#> black and brown, which differ by their pattern of polymer
#> bonds. A small amount of black eumelanin in the absence of
#> other pigments causes grey hair, while a small amount of
#> brown eumelanin in the absence of other pigm...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' An important and common pathway for
#> catabolism of amino acids by yeast is the Ehrlich pathway
#> |CITS: [Ehrlich07]|. In this pathway, following
#> transamination of an amino acid into the corresponding
#> 2-oxo acid, the 2-oxo acid is decarboxylated to an ald...'
#> to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: Lignins "Lignin"| is the most
#> abundant aromatic compound in nature, and its
#> mineralization is an important step in the terrestrial
#> carbon cycle. Lignin is composed of various intermolecular
#> linkages between phenylpropanes, including guaiacyl, ...'
#> to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Numerous filamentous fungi, including
#> the food biotechnology fungus |FRAME: TAX-5061|, the
#> opportunistic human pathogen |FRAME: TAX-746128|, the
#> phytopathogenic fungi |FRAME: TAX-40559| and |FRAME:
#> TAX-5180|, as well as many brown-rot and white-rot ba...'
#> to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'it can obtain its nitrogen from many
#> different sources, including most amino acids |CITS: [Large
#> 68]|. The most common mechanism for retreiving the nitrogen
#> from amino acids is transamination, using |FRAME:
#> 2-KETOGLUTARATE| or other 2-oxo acids as amin...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'it can obtain its nitrogen from many
#> different sources, including most amino acids |CITS: [Large
#> 68]|. The most common mechanism for retreiving the nitrogen
#> from amino acids is transamination, using |FRAME:
#> 2-KETOGLUTARATE| or other 2-oxo acids as amin...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'Even though the commercial production
#> of this chemical is very limited compared to |FRAME:
#> O-DICHLOROBENZENE| and |FRAME: 14-DICHLOROBENZENE|, its
#> concentration in the environment can be almost as high as
#> the other isomers (for example, 11 <b1> 13 ngli...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' They are both |FRAME:
#> Glycosaminoglycans "glycosaminoglycan"| polymers composed
#> of repeating disaccharide units, each containing |FRAME:
#> CPD-12557| linked by a (1,3) or (1,4) linkage to a uronic
#> acid, |FRAME: CPD-12521| or |FRAME: CPD-12515|, with e...'
#> to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: HEPARIN "Heparin"| is a
#> glycosaminoglycan polymer with a molecular weight ranging
#> from 3 to 30 kDa. It is composed of a variable sulfated
#> repeating disaccharide units composed of hexuronic (either
#> iduronic or glucuronic) acid and glucosamine r...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: Heparan-Sulfate "Heparan
#> sulfate"| (also known as heparitin sulfate) binds to
#> specific proteins such as antithrombin and several growth
#> factors, thereby regulating various biological processes
#> including anticoagulation and angiogenesis |CITS: ...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' They are both |FRAME:
#> Glycosaminoglycans "glycosaminoglycan"| polymers composed
#> of repeating disaccharide units, each containing |FRAME:
#> CPD-12557| linked by a (1,3) or (1,4) linkage to a uronic
#> acid, |FRAME: CPD-12521| or |FRAME: CPD-12515|, with e...'
#> to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: Hyaluronan Hyaluronan| is an
#> anionic, nonsulfated |FRAME: Glycosaminoglycans
#> "glycosaminoglycan"| distributed widely throughout
#> connective, epithelial, and neural tissues, including the
#> vitreous humor of the eye, the aorta, blood, brain,
#> liver...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Of the ten existing hexitols (|FRAME:
#> CPD-15703|, |FRAME: CPD-12811|, |FRAME: CPD-15702|, |FRAME:
#> SORBITOL|, |FRAME: CPD-12810|, |FRAME: MANNITOL|, |FRAME:
#> CPD-15701|, |FRAME: CPD-357|, |FRAME: CPD-369| and |FRAME:
#> GALACTITOL|) only three occur natura...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: D-APIOSE| is a unique
#> branched-chain pentose found principally in plants. It is a
#> key component of structurally complex cell wall
#> polysaccharides such as |FRAME: CPD-14581
#> "rhamnogalacturonan-II"| (RG-II) in the cell walls of
#> higher plants and...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: D-APIOSE| is a unique
#> branched-chain pentose found principally in plants. It is a
#> key component of structurally complex cell wall
#> polysaccharides such as |FRAME: CPD-14581
#> "rhamnogalacturonan-II"| (RG-II) in the cell walls of
#> higher plants and...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' The predominant circulating
#> pyrimidine in humans is |FRAME: URIDINE| |CITS:
#> [8155674][16769123]|. Among different species, including
#> man, its plasma level is strictly maintained at 3<96>5
#> &mu;M, a concentration higher than that of other
#> nucleosides |C...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' A higher-titer of |FRAME: BUTANOL|
#> production (15<96>30 gL) was achieved by replacing |FRAME:
#> G-93| with the |FRAME:MONOMER-16778| from |FRAME: TAX-3039|
#> (|FRAME:G-14619|), which directly utilizes NADH as the
#> reducing equivalent, and by artificially b...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Deoxyribonucleotides are synthesised
#> <i>de novo<i> at the diphosphate level through reduction of
#> the 2<92>-hydroxyl group of the corresponding
#> ribonucleotides (see |FRAME: PWY-7184|). This reduction is
#> mediated by the key enzyme |FRAME: EC-1.17.4.1 "r...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: ASCORBATE|, also known as
#> vitamin C, fulfils multiple essential roles in both plants
#> and animals. Being a strong reducing agent, it functions as
#> the primary water soluble antioxidant in cells, interacting
#> with reactive oxygen species generated...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Polyunsaturated long- and
#> very-long-chain fatty acids (of 20-22 carbons) accumulate
#> in many algae, mosses and ferns. They are found in all
#> major glycerol lipids of these organisms and represent
#> major constituents of membranes, although their function
#> ...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Teichoic acids were first reported in
#> 1958 |CITS: [Armstrong58]|, after being discovered
#> serendipitously as "contaminants" in bacterial extracts
#> from |FRAME: TAX-1590| |CITS: [2662967]|. As they make up a
#> significant portion of the cell walls (about 6...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Teichoic acids were first reported in
#> 1958 |CITS: [Armstrong58]|, after being discovered
#> serendipitously as "contaminants" in bacterial extracts
#> from |FRAME: TAX-1590| |CITS: [2662967]|. As they make up a
#> significant portion of the cell walls (about 6...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Saffron is is commonly considered the
#> most expensive spice on Earth. It consists of the dried red
#> stigmas of |FRAME: TAX-82528|, a perennial, sterile,
#> vegetatively propagated small plant from the |FRAME:
#> TAX-26339| family. The high price is a result ...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: CPD-13015| (MGG) is a
#> compatible solute identified in the slightly thermophilic
#> species |FRAME: TAX-28237| (optimum growth temperature of
#> 55 to 60<b0>C), a member of the order Thermotogales. Under
#> optimum growth conditions MGG is the major sol...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: CPD-13015|
#> (mannosylglucosylglycerate, MGG) is a compatible solute
#> identified in the slightly thermophilic species |FRAME:
#> TAX-28237| (optimum growth temperature of 55 to 60<b0>C), a
#> member of the order Thermotogales. Under optimum growth
#> cond...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Glycogen is a highly branched glucose
#> polymer that serves as a form of energy storage in animals
#> and fungi. It resembles plant-produced |FRAME: Starch
#> starch| and is sometimes called "animal starch".  Glycogen
#> is formed of small chains of 8 to 12 gluc...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Glycogen is a highly branched glucose
#> polymer that serves as a form of energy storage in animals
#> and fungi. It resembles plant-produced |FRAME: Starch
#> starch| and is sometimes called "animal starch".  Glycogen
#> is formed of small chains of 8 to 12 gluc...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: TAX-1773| and other
#> pathogenic mycobacterial species produce a capsular layer
#> outside the bacterial wall and plasma membrane that is
#> thought to be involved in immune evasion and virulence
#> |CITS: [13630872][13711244]|. The capsular layer is com...'
#> to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: ALGINATE Alginate| is a major
#> component of the cell walls of |FRAME: TAX-2870| (brown
#> algae), and can make up to 45% of the total algal material.
#> It is a linear polysaccharide of (1-4)-linked |FRAME:
#> CPD-15205|, with variable amounts of its C-...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: ALGINATE Alginate| is a major
#> component of the cell walls of |FRAME: TAX-2870| (brown
#> algae), and can make up to 45% of the total algal material.
#> It is a linear polysaccharide of (1-4)-linked |FRAME:
#> CPD-15205|, with variable amounts of its C-...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Cell-cell communication in bacteria
#> is accomplished through the exchange of extracellular
#> signalling molecules called autoinducers. This process,
#> termed quorum sensing, allows bacterial populations to
#> coordinate gene expression as a function of cell d...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Cell-cell communication in bacteria
#> is accomplished through the exchange of extracellular
#> signalling molecules called autoinducers. This process,
#> termed quorum sensing, allows bacterial populations to
#> coordinate gene expression as a function of cell d...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Cell-cell communication in bacteria
#> is accomplished through the exchange of extracellular
#> signalling molecules called autoinducers. This process,
#> termed quorum sensing, allows bacterial populations to
#> coordinate gene expression as a function of cell d...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Cell-cell communication in bacteria
#> is accomplished through the exchange of extracellular
#> signalling molecules called autoinducers. This process,
#> termed quorum sensing, allows bacterial populations to
#> coordinate gene expression as a function of cell d...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: PHOSPHATIDYLCHOLINE "A
#> phosphatidylcholine"| is |FRAME: Phosphoglycerides|
#> composed of a |FRAME: GLYCEROL| backbone esterified to
#> |FRAME: PHOSPHORYL-CHOLINE| and two |FRAME: Fatty-Acids
#> "fatty acids"|. Phosphatidylcholine is a major component
#> ...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: PHOSPHATIDYLCHOLINE "A
#> phosphatidylcholine"| is |FRAME: Phosphoglycerides|
#> composed of a |FRAME: GLYCEROL| backbone esterified to
#> |FRAME: PHOSPHORYL-CHOLINE| and two |FRAME: Fatty-Acids
#> "fatty acids"|. Phosphatidylcholine is a major component
#> ...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: PHOSPHATIDYLCHOLINE "A
#> phosphatidylcholine"| is |FRAME: Phosphoglycerides|
#> composed of a |FRAME: GLYCEROL| backbone esterified to
#> |FRAME: PHOSPHORYL-CHOLINE| and two |FRAME: Fatty-Acids
#> "fatty acids"|. Phosphatidylcholine is a major component
#> ...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: PHOSPHATIDYLCHOLINE "A
#> phosphatidylcholine"| is |FRAME: Phosphoglycerides|
#> composed of a |FRAME: GLYCEROL| backbone esterified to
#> |FRAME: PHOSPHORYL-CHOLINE| and two |FRAME: Fatty-Acids
#> "fatty acids"|. Phosphatidylcholine is a major component
#> ...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: PHOSPHATIDYLCHOLINE "A
#> phosphatidylcholine"| is |FRAME: Phosphoglycerides|
#> composed of a |FRAME: GLYCEROL| backbone esterified to
#> |FRAME: PHOSPHORYL-CHOLINE| and two |FRAME: Fatty-Acids
#> "fatty acids"|. Phosphatidylcholine is a major component
#> ...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: PHOSPHATIDYLCHOLINE "A
#> phosphatidylcholine"| is |FRAME: Phosphoglycerides|
#> composed of a |FRAME: GLYCEROL| backbone esterified to
#> |FRAME: PHOSPHORYL-CHOLINE| and two |FRAME: Fatty-Acids
#> "fatty acids"|. Phosphatidylcholine is a major component
#> ...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: PHOSPHATIDYLCHOLINE "A
#> phosphatidylcholine"| is |FRAME: Phosphoglycerides|
#> composed of a |FRAME: GLYCEROL| backbone esterified to
#> |FRAME: PHOSPHORYL-CHOLINE| and two |FRAME: Fatty-Acids
#> "fatty acids"|. Phosphatidylcholine is a major component
#> ...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'The difference among these compounds
#> stem from structural differences of an acyl group attached
#> to C2 of the &gamma;-lactam ring.  Slinosporamides A, C, D
#> and E contain a chloroethyl, ethyl, methyl and propyl
#> groups, respectively. A similar compound, |...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Bile (or gall) is a dark green to
#> yellowish brown fluid, produced by the liver of most
#> vertebrates, that aids the digestion of lipids in the small
#> intestine. In humans bile is produced continuously by the
#> liver and stored and concentrated in the gall ...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Bile (or gall) is a dark green to
#> yellowish brown fluid, produced by the liver of most
#> vertebrates, that aids the digestion of lipids in the small
#> intestine. In humans bile is produced continuously by the
#> liver and stored and concentrated in the gall ...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'worldwide. The pungency of the fruits
#> of red pepper is caused by a lipophilic alkaloid |FRAME:
#> CAPSAICIN| and its analogs, termed capsaicinoids |CITS:
#> [18190936]|. The fundamental structure of capsaicinoids is
#> that of an acid amide of |FRAME: CPD-932...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Natural rubber is produced from latex
#> - a material that makes up the cytoplasmic content of the
#> latecifers or latex vessels of |FRAME: TAX-3981| (the
#> Brazilian rubber tree) |CITS: [17545224]|. Although over
#> 2000 plant species capable of producing rubb...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Fagopyritols are galactosyl cyclitols
#> found primarily in seeds of buckwheat (|FRAME:TAX-3617|).
#> Whereas maturing embryos of many plant seeds accumulate
#> |FRAME:CPD-1099|, |FRAME:CPD-170| and |FRAME:CPD-8065| (in
#> addition to sucrose) as the predominan...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Fagopyritols are galactosyl cyclitols
#> found primarily in seeds of buckwheat (|FRAME:TAX-3617|).
#> Whereas maturing embryos of many plant seeds accumulate
#> |FRAME:CPD-1099|, |FRAME:CPD-170| and |FRAME:CPD-8065| (in
#> addition to sucrose) as the predominan...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' DAT and PAT contain one middle-chain
#> saturated fatty acid (|FRAME: PALMITATE| or |FRAME:
#> STEARIC_ACID|) at the 2-position and one or four
#> polymethyl-branched long-chain fatty acids at the 3-, 6-,
#> 2'-, 4'-, or 6'-position. The polymethyl-branched fatty...'
#> to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' The name ganglioside is derived from
#> the ganglion cells of the brain, from which the first
#> gangliosides were isolated by the German scientist Ernst
#> Klenk in 1942 |CITS: [Klenk42]|. Indeed, gangliosides
#> dominate the composition of the glycocalyces of m...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Organisms produce hydrocarbons of
#> different types by different mechanisms. Several mechanisms
#> have been described for the production of hydrocarbons from
#> fatty acids or their intermediates, including synthesis of
#> alkanes from fatty aldehydes by decarb...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Since acetaldehyde is reactive in the
#> atmosphere, its emisions is of interest to atmospheric
#> chemists, who are not able to account for large amounts of
#> biogenic |FRAME: ACETALD| |CITS: [19538397]|.  One source
#> of |FRAME: ACETALD| is |FRAME: ETOH|. In...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate 'which usually occurs as a glycoside
#> with a deoxy(amino)sugar.  Elloramycin is produced by
#> |FRAME: TAX-47716| T<fc>2353, and was found to possess a
#> weak activity against a variety of Gram-positive bacteria
#> (especially streptomycetes) and against stem ce...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: STAUROSPORINE
#> "Staurosporine"| was isolated in 1977 (under the name
#> AM-2282) from a culture of |FRAME: TAX-65499 "<i>Lentzea
#> albida<i> (then known as <i>Streptomyces stauroporeus<i>
#> AM-2282)"| that was collected from a soil sample in Japan
#> |CI...' to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: Pseudomonic-Acids
#> "Pseudomonic acids"| are polyketide antibiotics produced by
#> |FRAME: ORG-6360| |CITS: [5003547]|. The antibiotic
#> mupirocin is a mixture of four pseudomonic acids: |FRAME:
#> CPD0-2581| (about 90%), |FRAME: CPD-21448| (about 8%), ...'
#> to a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Green sulfur bacteria (|FRAME:
#> TAX-191412|) and green flimentous bacteria (|FRAME:
#> TAX-1106|) contain chlorosomes - large organelles that are
#> believed to be the most efficient light-harvesting antennae
#> in nature |CITS: [20130996][Saga10]|. Green sulfu...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Green sulfur bacteria (|FRAME:
#> TAX-191412|) and green flimentous bacteria (|FRAME:
#> TAX-1106|) contain chlorosomes - large organelles that are
#> believed to be the most efficient light-harvesting antennae
#> in nature |CITS: [20130996][Saga10]|. Green sulfu...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Green sulfur bacteria (|FRAME:
#> TAX-191412|) and green flimentous bacteria (|FRAME:
#> TAX-1106|) contain chlorosomes - large organelles that are
#> believed to be the most efficient light-harvesting antennae
#> in nature |CITS: [20130996][Saga10]|. Green sulfu...' to a
#> wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: Cobalamins Cobalamin|-derived
#> cofactors are required by two important human enzymes -
#> |FRAME: EC-2.1.1.13|, and |FRAME: EC-5.4.99.2|. The former
#> utilizes |FRAME: CPD-9037|, while the latter requires
#> |FRAME: ADENOSYLCOBALAMIN|.  While humans ca...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Disproportionation is a process in
#> which a reactant is both oxidized and reduced in the same
#> chemical reaction, forming two separate compounds. The
#> disproportionation of |FRAME: S2O3|, an important step in
#> sulfur transformation in aquatic systems, was...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' |FRAME: Ulvan Ulvans| are cell wall
#> matrix polysaccharides in |FRAME: TAX-3041 "green algae"|
#> belonging to the genus |FRAME: TAX-3118|. They represent
#> 8<96>30% of the dry weight of Ulva species.  Ulvans have a
#> complex structure that has not been compl...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' Arabinoxylans have been identified in
#> wheat, rye, barley, oat, rice, and sorghum, as well as in
#> some less common plants, including pangola grass, bamboo
#> shoots and rye grass. |FRAME: Glucuronoxylans
#> Glucuronoxylans| and |FRAME: Glucuronoarabinoxylans ...' to
#> a wide string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> unable to translate ' A bacterial catabolic pathway of DCA
#> was proposed in |FRAME: TAX-13689| TMY1009 |CITS:
#> [Habu88]|. Several mutant strains were created and shown to
#> accumulated several compounds including |FRAME: CPD-17077|
#> (DCA-C), |FRAME: CPD-17080| (DCA-S), |FRAME: ...' to a wide
#> string
#> Warning in grepl(paste(candSp, collapse = "|"), coms):
#> input string 1 is invalid
metacyc |> head()
#>   pathwayID
#> 1  PWY-7622
#> 2  PWY-7591
#> 3  PWY-7613
#> 4  PWY-7529
#> 5  PWY-7599
#> 6  PWY-7536
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          text
#> 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  D-galactofuransoe is the five-carbon ring form of D-galactopyranose. In solution both forms exist in equilibrium with D-galactopyranose formation favored. Although D-galactopyranose is ubiquitous in cellular organisms, D-galactofuransoe is found only in some non-mammalian eukaryotes and in some bacteria including mycobacteria and Escherichia coli (as indicated in the pathway links and E. coli enzyme EG11983-MONOMER).  UDP-D-GALACTO-14-FURANOSE is a nucleotide-activated form of D-galactofuransoe that is used by some organisms in the biosynthesis of polysaccharides and glycoconjugates. Although D-galactofuransoe residues are not found in mammals, higher plants and yeasts, they are found in the glycans of some bacteria and lower eukaryotes including trypanosomatids, nematodes, the free-living alga TAX-3055 and filamentous fungi. These groups include known pathogens. D-galactofuransoe residues are antigenic in humans and are therefore of interest as therapeutic targets ( and reviewed in ).  About This Pathway  In filamentous fungi UDP-D-GALACTO-14-FURANOSE is biosynthesized in the cytosol starting with CPD-12575 which is derived from GLC-6-P as indicated in the pathway link. Following the epimerization of CPD-12575 to CPD-14553, EC-5.4.99.9 converts CPD-14553 to its five-carbon ring form UDP-D-GALACTO-14-FURANOSE. The substrate of the mutase, CPD-14553, may also be formed by the Leloir pathway enzyme EC-2.7.7.12. In protozoan parasites such as TAX-5664 which appear to lack an ortholog encoding the Leloir pathway enzyme, CPD-14553 can be formed by EC-2.7.7.64 . The pathway product UDP-D-GALACTO-14-FURANOSE is then transported from the cytosol to the Golgi lumen for glycoconjugate biosynthesis (reviewed in ) (see pathway PWY-6317).   In filamentous fungi galactomannan metabolism has been studied in both the non-pathogenic model organism TAX-162425 and in the opportunistic pathogen TAX-746128, and relevant orthologs have been identified in these organisms. In TAX-162425 cytosolic UDP-glucose 4-epimerase (synonym UDP-galactose 4-epimerase) is encoded by ugeA . In TAX-746128 uge5 encodes the dominant cytosolic UDP-glucose 4-epimerase that is essential for growth on D-Galactose and the synthesis of D-galactofuransoe. In TAX-746128 a second gene uge3 with no identified ortholog in TAX-162425 has also been characterized. Both Uge5 and Uge3 are required for galactosaminogalactan synthesis in TAX-746128 . In TAX-162425 cytosolic UDP-galactose mutase is encoded by ugmA, and in TAX-746128 by ugm1 (glfA) .   UDP-D-GALACTO-14-FURANOSE is then transported from the cytosol to the Golgi lumen for galactomannoprotein biosynthesis (not shown). In TAX-162425 the UDP-galactofuranose transporter is encoded by ugtA , and in TAX-746128 by ugt1 (glfB) ( and reviewed in ). In both organisms a novel, G-55162-encoded, Golgi-localized, EC-2.4.1.M27, has been characterized that is involved in synthesis of the galactofurnaose antigen of O-linked glycans. The O-glycans are then transported in vesicles to the hyphal cell surface . 
#> 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             CPD-17304 Okenone is a red carotenoid pigment. Structurally this ketocarotenoid has a &chi;-ring at one end and an open chain &psi;-end that is methoxylated at the C-1' position and contains a keto group at the C-4' position. It is found in some purple sulfur bacterial members of the family TAX-1046. These photosynthetic organisms reside in illuminated anoxic zones of aquatic habitats. Purple sulfur bacteria are anaerobic or microaerophilic, oxidizing HS to produce granules of Elemental-Sulfur (elemental sulfur) (see pathway PWY-5274). The light absorption properties of CPD-17304 allow these bacteria to exist in deeper water layers. In addition to its role in bacterial ecology, CPD-17304 also becomes an important geochemical biomarker via its diagenesis product CPD-17311, which is found in ancient rock formations (  and in .  TAX-572262 is a purple sulfur bacterium that synthesizes CPD-17304 as its only carotenoid when grown anoxically under chlorophototrophic conditions. A gene cluster encoding the enzymes for CPD-17304 biosynthesis (genes crtE, crtB,  crtI, crtC, crtU and crtY and an unlinked gene crtF) was identified and analyzed. The biosynthetic pathway for CPD-17304 was elucidated based on heterologous expression of recombinant enzymes in carotenoid-producing hosts, followed by carotenoid analysis .  About This Pathway  In this pathway CrtY was shown to catalyze the cyclization of one &psi;-end of CPD1F-114 producing the &beta;-ring of CPD1F-126. The CrtC hydratase hydrates the double bond of the &psi;-end group of CPD1F-126 producing CPD-11446. CruO is a unique C-44' ketolase required for CPD-17304 biosynthesis that introduces a keto group at the C-4' position of CPD-11446. CrtF is an O-methyltransferase that further modifies the &psi;-end group by transfer of a methyl group to the hydroxyl formed by CrtC, producing the methoxy group of CPD-17302. Finally, the &beta;-ring of CPD-17302 is converted to a &chi;-ring and the molecule is further desaturated by the carotene desaturasemethyltransferase CrtU to produce CPD-17304. The key enzymes of the pathway are CrtY and CrtU .  In addition, a putative OXYGEN-MOLECULE-dependent pathway for CPD-17307 biosynthesis involving genes crtC, cruS and crtF was also proposed (not shown), based on the demonstration of these enzymatic functions in Escherichia coli. However, no physiological evidence for the pathway could be found in either TAX-572262 or TAX-244573. The novel gene cruS was shown to encode a unique 2-ketolase3,4-desaturase that could participate in the hypothetical pathway . 
#> 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     Heptose sugars, mostly glycero-manno-heptoses, are found in the cell surface polymers of many bacteria. Prior to their incorporation, the heptose residues are activated by attachment to a nucleotide (see pathways PWY0-1241 and PWY-6478). These activated heptoses can be uniquely modified by enzymes before incorporation into the polymers. In bacterial pathogens heptose derivatives such as 6-deoxyheptoses that are not found in mammals can play an important role in virulence, as demonstrated in TAX-633 . Heptose modifying enzymes are therefore of interest as potential drug targets .  The modified heptose CPD-17411 is found in the O-Antigens O-antigen of ORG-6272. The O-antigen along with core oligosaccharide and Lipid-A lipid A comprise Lipopolysaccharides lipopolysaccharide (LPS), a major component of the outer membrane in Gram-negative bacteria. LPS is one of several virulence factors in TAX-633. It has been shown that 6-deoxy-D-manno-heptose affects the barrier function of LPS and the overall virulence of ORG-6272 .  About This Pathway  DmhA along with its corresponding reductase DmhB catalyze the formation of 6-deoxyheptose, a modified heptose present in the O-antigen of ORG-6272. In this pathway both DmhA and DmhB utilize and release sugars that are in the D-manno configuration. Unlike PWY-7610 in TAX-197, no epimerization steps are involved .  Both dmhA and dmhB mutants were analyzed for the composition and structure of their LPS and the virulence-related properties of these mutants and their complemented counterparts were assessed in vitro . Recombinant proteins were expressed in Escherichia coli, purified and biochemically characterized . 
#> 4  Sialic acids are a family of polyhydroxylated &alpha;-keto acids that contain nine carbon atoms. Most sialic acids are derivatives of N-ACETYLNEURAMINATE, or CPD-10734 (KDN). N-ACETYLNEURAMINATE is the most common sialic acid in mammals (see pathway PWY-6138), while KDN is abundant in lower vertebrates (see pathway PWY-6140). Their core structures can be modified at the hydroxyl groups, lactonized, or hydroxylated at the acetamido group, generating many derivatives. CPD-262 is a derivative of CMP-N-ACETYL-NEURAMINATE (see pathway PWY-6144). Reviewed in    and .  Sialic acids are found mainly in vertebrates and a few higher invertebrates (ascidians and echinoderms). These acidic sugars are usually the terminal sugar residue in the glycan chains of vertebrate glycoconjugates (mostly glycoproteins and glycolipids, but also proteoglycans and glycosylphosphatidylinositol anchors). They function in mediating cellular recognition and adhesion events for many important processes such as development, the immune and inflammatory responses, and oncogenesis. Sialic acid occurs rarely in invertebrates. Endogenous sialylation has been shown to occur in TAX-7215, but details of sialic acid biosynthesis in this organism remain to be determined . However, it is possible the sialic acids might be biosynthesized by other eukaryotes in a species andor life cycle-dependent manner. In  and reviewed in    and .  Most bacteria do not biosynthesize sialic acids, but some pathogenic, or symbiotic bacteria biosynthesize sialic acids as a means of evading a host's immune system (this pathway). The sialic acid is displayed on the bacterial cell surface (in capsular polysaccharides) in order to mimic vertebrate cells. Pathogens that biosynthesize sialic acids include TAX-487, TAX-1392869 and TAX-197. In addition, the human gut symbiont TAX-818 has been shown to synthesize CPD-10734  (see pathway PWY-6140). Whether or not archaea contain sialic acids remains to be determined (reviewed in    and ). Other sialic acid-like sugars biosynthesized by bacteria include the nonulosonic acids CPD-10754  (see pathway PWY-6143) and legionaminic acid see pathway class CMP-Legionaminate-Biosynthesis .  Protists are thought to lack the ability to biosynthesize sialic acids although more genome data are needed to confirm this. Sialic acids have been thought to be absent in plants but some studies raise the possibility . Fungi appear to lack any known sialic acid biosynthetic pathway, although strain-specific, or novel pathways could exist. Reviewed in    and .  Also see PWY-6145.  About This Pathway  TAX-316275 is a Gram-negative, psychrophilic fish pathogen whose genome sequence has been determined . Little is known about its mechanism of virulence. This organism contains three copies of a gene cluster that is homologous to the Escherichia coli neu gene cluster for the synthesis of capsular sialic acids (see PWY-6139). In TAX-316275 the first copy of the gene cluster, neu1, is also likely involved in sialic acids biosynthesis. The second copy neu2, which is less homologous, is likely involved in the synthesis of alternative compounds such as legionaminate (see pathway class CMP-Legionaminate-Biosynthesis). The third copy is an exact duplicate of the neu2 gene cluster .  Data suggest that both N-ACETYLNEURAMINATE and CPD-1803 (7-O-acetyl-N-acetylneuraminate) are present in this organism and that their synthesis is catalyzed by the sialic acid synthase product of gene neuB1 which is present in the neu1 gene cluster. Based on genome analysis, enzyme kinetics, and structural analysis, the putative substrate for NeuB1 is  CPD-16880 and its product is CPD-1803. In addition to TAX-316275, this gene was also predicted to be present in other species including Escherichia coli and TAX-1311.   A pathway for the biosynthesis of CPD-1803 has been proposed and is shown here. In addition to the experimentally determined NeuB1, predicted enzymes of the pathway include the O-acetyl transferase NeuD1, the hydrolyzing 2-epimerase NeuC1, and NeuA1 which activates CPD-1803 to a CMP derivative .  
#> 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          CPD-17343 Anditomin is CPD-17353. These complex organo-oxygen natural products are produced by fungi from polyketide and terpenoid precursors. They have unique, highly oxygenated structures containing a complex bridged-ring system. The elucidation of the biosynthetic pathway for the meroterpenoid CPD-17343 in TAX-1549217 provides an opportunity for future construction of novel scaffolds for use in drug discovery .  About This Pathway  The CPD-17343 biosynthetic gene cluster of TAX-1549217 was identified by bioinformatic analysis and the and gene products were functionally characterized (note that NCBI Taxonomy merged the Emericella variecolor and Aspergillus stellifer entries). The early-stage and late-stage biosynthetic steps were determined by expression of recombinant enzymes in a strain of TAX-5062 followed by analysis of products by HPLC, 1H NMR, 13C NMR and mass spectrometry. The mid-stage biosynthetic steps were elucidated using feeding experiments to predict the sequence of tailoring enzymes. In the case of AndA and AndF, recombinant enzyme was expressed in Escherichia coli, purified and characterized .  The pathway begins with a polyketide synthase encoded by gene andM that produces CPD-17316 as its final product. AndK is a bifunctional P450 monooxygenase and hydrolase fusion protein that produces the phthalide compound CPD-17317. AndD is a prenyltransferase that incorporates the farnesyl moiety. AndE is an epoxidase that forms the (S)-epoxide CPD-17325. The terpene cyclase AndB then forms CPD-17333 .   The mid-stage biosynthetic steps include AndA, AndJ and AndI. The nonheme iron-dependent dioxygenase AndA was shown to catalyze two reactions, the dehydrogenation of CPD-17334 to produce the enone CPD-17335 that contains a &Delta;1,2-conjugated double bond, and an isomerization involving an unprecedented skeletal rearrangement that results in the bridged-ring of CPD-17336. This is in contrast to a previously proposed DielsAlder reaction. AndA thus generates the scaffold of the andilesins class of secondary metabolites. AndJ is a FAD-dependent Baeyer-Villiger monooxygenase that generates a seven-membered lactone ring from CPD-17336. AndI is a short-chain dehydrogenasereductase (SDR) that appears to reduce CPD-17337 to CPD-17339 .   In the late-stage biosynthetic steps acetyltransferase AndG attaches an acetyl group to the hydroxyl group of CPD-17339. This acetyl group is then lost in a spontaneous reaction. AndH is another SDR reductase that catalyzes reduction of the C-6 double bond of CPD-17341 to produce CPD-17342. Finally AndF, another nonheme iron-dependent dioxygenase like AndA, oxidizes CPD-17342 to CPD-17343. Based on the cofactor requirements for the production of both CPD-17336 by AndA and CPD-17343 by AndF, as well as the analogous FtmF reaction (see EC-1.14.11.38), detailed reaction mechanisms for AndA and AndF were proposed . 
#> 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    This pathway depicts the biosynthesis of the CPD-16935 moiety of the natural product CPD-16932. The reaction link shows its amide linkage to the synthetic partial polyenoate scaffold CPD-16939 which was derived from the structure of CPD-16932 .  The antifungal agent CPD-16932 is a linear polyketide antibiotic produced by ORG-6265. In the CPD-16932 biosynthetic gene cluster of this organism three tandem ORFs, ORF33, ORF34 and ORF35 catalyze the formation of the CPD-16935 moiety of CPD-16932. Other ORFs in the cluster encode a large type I polyketide synthetase 15844935. CPD-16935 2-Amino-3-hydroxycyclopent-2-enone (an enol tautomer of CPD-16936) is a C5N unit that is present in many members of the manumycin antibiotic family as well as in other bioactive metabolites. In CPD-16932 the amino group of CPD-16935 is amide-bonded to the polyketide-derived polyenoic acid component, whereas in MOENOMYCIN it is amide-bonded to Hexuronates group .  About This Pathway  In the first reaction GLY is condensed with SUC-COA which produces 5-AMINO-LEVULINATE, catalyzed by  MONOMER-18786 encoded by ORF34. 5-AMINO-LEVULINATE 5-Aminolevulinate is also an intermediate in tetrapyrrole biosynthesis in some species of TAX-1883  (see pathway PWY-5189). In the next reaction an ATP-dependent MONOMER-18787 encoded by ORF35 catalyzes the ligation of 5-AMINO-LEVULINATE to CO-A producing CPD-16937. The latter compound is unstable and can either undergo spontaneous intramolecular cyclization to produce the off-pathway shunt product CPD-16938, or be converted enzymatically to CPD-16935 by the cyclase activity of the ORF34 product .  The reaction link shows a condensation step catalyzed by the MONOMER-18788 product of ORF33 that is proposed to incorporate CPD-16935 into the polyketide chain of CPD-16932 which constitutes the final chain termination step. This reaction was experimentally demonstrated using the synthetic partial polyenoate scaffold CPD-16939 which was derived from the structure of CPD-16932 .  Recombinant proteins encoded by ORF33, ORF34 and ORF35 were overexpressed in Escherichia coli, purified and biochemically characterized. In addition, the entire CPD-16935 biosynthetic pathway was reconstituted in vitro. The three enzymes were incubated with MG+2, ATP, CO-A, GLY and SUC-COA (or 5-AMINO-LEVULINATE) and CPD-16939. The reaction product CPD-16941 was identified by LC-MS. Gel filtration chromatography demonstrated that none of the three enzymes formed complexes, indicating freely diffusible intermediates .    
#>                                                      commonName
#> 1                    UDP-&alpha;-D-galactofuranose biosynthesis
#> 2                                          okenone biosynthesis
#> 3               GDP-6-deoxy-D-<i>manno</i>-heptose biosynthesis
#> 4 CMP-<i>N</i>-acetyl-7-<i>O</i>-acetylneuraminate biosynthesis
#> 5                                        anditomin biosynthesis
#> 6               2-amino-3-hydroxycyclopent-2-enone biosynthesis
#>                            species
#> 1 TAX-746128,TAX-330879,TAX-162425
#> 2  TAX-572262,TAX-37487,TAX-244573
#> 3                         ORG-6272
#> 4                       TAX-316275
#> 5                      TAX-1549217
#> 6                         ORG-6265
#>                                taxonomicRange
#> 1 TAX-6231,TAX-3052,TAX-147538,TAX-5654,TAX-2
#> 2                                    TAX-1046
#> 3                                       TAX-2
#> 4                                       TAX-2
#> 5                                    TAX-4751
#> 6                                  TAX-201174
#>              query
#> 1 Escherichia coli
#> 2 Escherichia coli
#> 3 Escherichia coli
#> 4 Escherichia coli
#> 5 Escherichia coli
#> 6 Escherichia coli
```
