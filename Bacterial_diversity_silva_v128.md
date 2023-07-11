# Analyses of bacterial diversity


```R
#Load libraries
library(scales)
library(phyloseq)
library(phylogeo)
library(ggplot2)
library(ape)
library(vegan)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(forcats)
```

    Warning message:
    “replacing previous import ‘dplyr::combine’ by ‘gridExtra::combine’ when loading ‘phylogeo’”
    Loading required package: permute
    
    Loading required package: lattice
    
    This is vegan 2.5-7
    



```R
# Load OTU table
otu <- as.matrix(read.table("biofert_16S_otu.tsv", header=T, row.names=1))
colnames(otu) <- c("root","substrate")
otu <- otu[, c("substrate", "root")]
OTU <- otu_table(otu, taxa_are_rows=T)
```


```R
# Load taxonomy table
taxa <- as.matrix(read.table("biofert_16S_tax.tsv", row.names = 1))
TAXA <- tax_table(taxa)
colnames(TAXA) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
```


```R
# Create phyloseq object with taxonomy and OTU tables 
bact <-phyloseq(OTU,TAXA)
bact
```


    phyloseq-class experiment-level object
    otu_table()   OTU Table:         [ 12249 taxa and 2 samples ]
    tax_table()   Taxonomy Table:    [ 12249 taxa by 7 taxonomic ranks ]



```R
# Rarefaction curves
pdf("rarefaction.pdf")
rarefaction <- rarecurve(t(otu), step=1000, cex=0.5, col="red")
dev.off()

rarecurve(t(otu), step=1000, cex=0.5, col="red")
```


<strong>png:</strong> 2



    
![png](output_5_1.png)
    



```R
# Shannon diversity estimations
estimate_richness(bact)
write.table(estimate_richness(bact), file="bact_diversity.tsv",append = FALSE, quote = TRUE, sep = "\t")
```


<table class="dataframe">
<caption>A data.frame: 2 × 9</caption>
<thead>
	<tr><th></th><th scope=col>Observed</th><th scope=col>Chao1</th><th scope=col>se.chao1</th><th scope=col>ACE</th><th scope=col>se.ACE</th><th scope=col>Shannon</th><th scope=col>Simpson</th><th scope=col>InvSimpson</th><th scope=col>Fisher</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>substrate</th><td>6971</td><td>7359.233</td><td>29.58714</td><td>7871.580</td><td>39.12867</td><td>6.080002</td><td>0.9895416</td><td> 95.61659</td><td>1634.206</td></tr>
	<tr><th scope=row>root</th><td>8819</td><td>8982.350</td><td>16.48938</td><td>9393.791</td><td>39.29644</td><td>7.617671</td><td>0.9981706</td><td>546.63185</td><td>2426.967</td></tr>
</tbody>
</table>




```R
div <- read.table("bact_diversity.tsv", header=TRUE, row.names = NULL)
colnames(div) <- c("Sample", "Observed", "Chao1", "SE.Chao1", 
                    "ACE", "SE.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")

observed_plot <- ggplot(div, aes(reorder(Sample, Observed), y=Observed))  + geom_point() + 
                        theme_light() +  ggtitle("Observed diversity") 
ggsave("bact_observed_plot.pdf", width=10, height=5, units="cm")
observed_plot


shannon_plot <- ggplot(div, aes(reorder(Sample, Shannon), y=Shannon))  + geom_point() + 
                        theme_light() + ggtitle("Shannon diversity")
ggsave("bact_shannon_plot.pdf", width=10, height=5, units="cm")
shannon_plot

simpson_plot <- ggplot(div, aes(reorder(Sample, Simpson), y=Simpson))  + geom_point() + 
                        theme_light() + ggtitle("Simpson diversity")
ggsave("bact_simpson_plot.pdf", width=10, height=5, units="cm")
simpson_plot
```


    
![png](output_7_0.png)
    



    
![png](output_7_1.png)
    



    
![png](output_7_2.png)
    



```R
#Estimate relative abundance
rel_bact <- transform_sample_counts(bact, function(x) x / sum(x))
```


```R
#Most abundant phylum
#Tax glom at Phylum level
phy_bact <- tax_glom(rel_bact, "Phylum")
phy_bact
```


    phyloseq-class experiment-level object
    otu_table()   OTU Table:         [ 41 taxa and 2 samples ]
    tax_table()   Taxonomy Table:    [ 41 taxa by 7 taxonomic ranks ]



```R
#Export the file with taxa names and taxa counts
bact_otu_phy <- otu_table(phy_bact)
bact_tax_phy <- tax_table(phy_bact)
phy_bact_tab <- cbind(bact_otu_phy, bact_tax_phy)
phy_bact_tab <- phy_bact_tab[, 1:4]
write.table(phy_bact_tab, "phy_bact_tab.tsv", sep = "\t")
phy_bact_tab <- read.table("phy_bact_tab.tsv", header=TRUE, row.names=1, stringsAsFactors = FALSE)
m_phy_bact_tab <- melt(phy_bact_tab)
colnames(m_phy_bact_tab) <- c("Kingdom", "Phylum", "Sample", "Relative_abundance")
m_phy_bact_tab <- m_phy_bact_tab[, 2:4]

#Collapse Phylum with low abundance
m_phy_bact_tab$Phylum[m_phy_bact_tab$Relative_abundance <= 0.015] <- "Low_abundance"
cm_phy_bact_tab <- aggregate(m_phy_bact_tab$Relative_abundance 
                                          ,by=list(m_phy_bact_tab$Phylum, 
                                           m_phy_bact_tab$Sample),sum)
colnames(cm_phy_bact_tab) <- c("Phylum", "Sample", "Relative_abundance")
cm_phy_bact_tab

#Plot phylum
print("Most abundant Phylum")
phylum_plot <- ggplot(cm_phy_bact_tab, aes(x=Sample, y=Relative_abundance,
                                          fill=fct_reorder(Phylum, Relative_abundance))) + 
       scale_fill_manual(values=rev(colorRampPalette(brewer.pal(n = 8, name = "Set1"))(14))) +
       geom_bar(stat="identity", color="black", width=0.8) + scale_y_continuous(expand = c(0 ,0)) + 
       theme_light(base_size = 7) + theme(axis.text.y = element_text(size = 15),
                                          axis.text.x = element_text(size= 15)) 
ggsave("bact_phylum_plot.pdf", width=20, height=20, units="cm")

phylum_plot
```

    Using Kingdom, Phylum as id variables
    



<table class="dataframe">
<caption>A data.frame: 19 × 3</caption>
<thead>
	<tr><th scope=col>Phylum</th><th scope=col>Sample</th><th scope=col>Relative_abundance</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Acidobacteria   </td><td>substrate</td><td>0.02130849</td></tr>
	<tr><td>Actinobacteria  </td><td>substrate</td><td>0.01619271</td></tr>
	<tr><td>Bacteroidetes   </td><td>substrate</td><td>0.27979049</td></tr>
	<tr><td>Low_abundance   </td><td>substrate</td><td>0.04783734</td></tr>
	<tr><td>Proteobacteria  </td><td>substrate</td><td>0.60578859</td></tr>
	<tr><td>Saccharibacteria</td><td>substrate</td><td>0.02908238</td></tr>
	<tr><td>Acidobacteria   </td><td>root     </td><td>0.10147018</td></tr>
	<tr><td>Actinobacteria  </td><td>root     </td><td>0.03670412</td></tr>
	<tr><td>Bacteroidetes   </td><td>root     </td><td>0.21991168</td></tr>
	<tr><td>Chloroflexi     </td><td>root     </td><td>0.01659120</td></tr>
	<tr><td>Cyanobacteria   </td><td>root     </td><td>0.01711666</td></tr>
	<tr><td>Firmicutes      </td><td>root     </td><td>0.02219241</td></tr>
	<tr><td>Gemmatimonadetes</td><td>root     </td><td>0.02157751</td></tr>
	<tr><td>Low_abundance   </td><td>root     </td><td>0.03665940</td></tr>
	<tr><td>Parcubacteria   </td><td>root     </td><td>0.01547320</td></tr>
	<tr><td>Planctomycetes  </td><td>root     </td><td>0.03976745</td></tr>
	<tr><td>Proteobacteria  </td><td>root     </td><td>0.34503885</td></tr>
	<tr><td>Saccharibacteria</td><td>root     </td><td>0.05405556</td></tr>
	<tr><td>Verrucomicrobia </td><td>root     </td><td>0.07344178</td></tr>
</tbody>
</table>



    [1] "Most abundant Phylum"



    
![png](output_10_3.png)
    



```R
#Most abundant phylum
#Tax glom at Phylum level
class_bact <- tax_glom(rel_bact, "Class")
class_bact
```


    phyloseq-class experiment-level object
    otu_table()   OTU Table:         [ 139 taxa and 2 samples ]
    tax_table()   Taxonomy Table:    [ 139 taxa by 7 taxonomic ranks ]



```R
#Export the file with taxa names and taxa counts
bact_otu_class <- otu_table(class_bact)
bact_tax_class <- tax_table(class_bact)
class_bact_tab <- cbind(bact_otu_class, bact_tax_class)
head(class_bact_tab)
```


<table class="dataframe">
<caption>A matrix: 6 × 9 of type chr</caption>
<thead>
	<tr><th></th><th scope=col>substrate</th><th scope=col>root</th><th scope=col>Kingdom</th><th scope=col>Phylum</th><th scope=col>Class</th><th scope=col>Order</th><th scope=col>Family</th><th scope=col>Genus</th><th scope=col>Species</th></tr>
</thead>
<tbody>
	<tr><th scope=row>25</th><td>6.10059001420566e-05</td><td>0.00202358991559059 </td><td>Bacteria</td><td>Planctomycetes  </td><td>vadinHA49              </td><td>NA</td><td>NA</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>27</th><td>0.000130726928875836</td><td>0.0063055508971994  </td><td>Bacteria</td><td>Gemmatimonadetes</td><td>S0134_terrestrial_group</td><td>NA</td><td>NA</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>76</th><td>0.00725970211690473 </td><td>0.00684219352674828 </td><td>Bacteria</td><td>Verrucomicrobia </td><td>Spartobacteria         </td><td>NA</td><td>NA</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>129</th><td>0.00158615340369347 </td><td>0.00470680306333501 </td><td>Bacteria</td><td>Gemmatimonadetes</td><td>Longimicrobia          </td><td>NA</td><td>NA</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>237</th><td>0                   </td><td>7.8260383475879e-05 </td><td>Bacteria</td><td>Verrucomicrobia </td><td>WCHB1-41               </td><td>NA</td><td>NA</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>264</th><td>0                   </td><td>8.94404382581475e-05</td><td>Bacteria</td><td>Latescibacteria </td><td>Ambiguous_taxa         </td><td>NA</td><td>NA</td><td>NA</td><td>NA</td></tr>
</tbody>
</table>




```R
#Export the file with taxa names and taxa counts
bact_otu_class <- otu_table(class_bact)
bact_tax_class <- tax_table(class_bact)
class_bact_tab <- cbind(bact_otu_class, bact_tax_class)
class_bact_tab <- class_bact_tab[, 1:5]
write.table(class_bact_tab, "class_bact_tab.tsv", sep = "\t")
class_bact_tab <- read.table("class_bact_tab.tsv", header=TRUE, row.names=1, stringsAsFactors = FALSE)
m_class_bact_tab <- melt(class_bact_tab)
colnames(m_class_bact_tab) <- c("Kingdom", "Phylum", "Class", "Sample", "Relative_abundance")

#Collapse Class with low abundance
m_class_bact_tab$Class[m_class_bact_tab$Relative_abundance <= 0.015] <- "Low_abundance"
cm_class_bact_tab <- aggregate(m_class_bact_tab$Relative_abundance 
                                          ,by=list(m_class_bact_tab$Class, 
                                           m_class_bact_tab$Sample),sum)
colnames(cm_class_bact_tab) <- c("Class", "Sample", "Relative_abundance")
cm_class_bact_tab

#Plot Class
print("Most abundant Class")
Class_plot <- ggplot(cm_class_bact_tab, aes(x=Sample, y=Relative_abundance,
                                          fill=fct_reorder(Class, Relative_abundance))) + 
       scale_fill_manual(values=rev(colorRampPalette(brewer.pal(n = 8, name = "Spectral"))(17))) +
       geom_bar(stat="identity", color="black", width=0.8) + scale_y_continuous(expand = c(0 ,0)) + 
       theme_light(base_size = 7) + theme(axis.text.y = element_text(size = 15),
                                          axis.text.x = element_text(size= 15)) 
ggsave("bact_Class_plot.pdf", width=20, height=20, units="cm")

Class_plot
```

    Using Kingdom, Phylum, Class as id variables
    



<table class="dataframe">
<caption>A data.frame: 24 × 3</caption>
<thead>
	<tr><th scope=col>Class</th><th scope=col>Sample</th><th scope=col>Relative_abundance</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Alphaproteobacteria </td><td>substrate</td><td>0.24133934</td></tr>
	<tr><td>Betaproteobacteria  </td><td>substrate</td><td>0.30849812</td></tr>
	<tr><td>Cytophagia          </td><td>substrate</td><td>0.01980077</td></tr>
	<tr><td>Gammaproteobacteria </td><td>substrate</td><td>0.04393296</td></tr>
	<tr><td>Low_abundance       </td><td>substrate</td><td>0.11293064</td></tr>
	<tr><td>Sphingobacteriia    </td><td>substrate</td><td>0.25574545</td></tr>
	<tr><td>uncultured_bacterium</td><td>substrate</td><td>0.01775272</td></tr>
	<tr><td>Actinobacteria      </td><td>root     </td><td>0.01779865</td></tr>
	<tr><td>Alphaproteobacteria </td><td>root     </td><td>0.12926379</td></tr>
	<tr><td>Betaproteobacteria  </td><td>root     </td><td>0.09999441</td></tr>
	<tr><td>Blastocatellia      </td><td>root     </td><td>0.03447929</td></tr>
	<tr><td>Cyanobacteria       </td><td>root     </td><td>0.01560736</td></tr>
	<tr><td>Cytophagia          </td><td>root     </td><td>0.03977863</td></tr>
	<tr><td>Deltaproteobacteria </td><td>root     </td><td>0.05544189</td></tr>
	<tr><td>Gammaproteobacteria </td><td>root     </td><td>0.06009279</td></tr>
	<tr><td>Low_abundance       </td><td>root     </td><td>0.17924982</td></tr>
	<tr><td>OPB35_soil_group    </td><td>root     </td><td>0.04117614</td></tr>
	<tr><td>Opitutae            </td><td>root     </td><td>0.02327687</td></tr>
	<tr><td>Phycisphaerae       </td><td>root     </td><td>0.01706076</td></tr>
	<tr><td>Planctomycetacia    </td><td>root     </td><td>0.01947566</td></tr>
	<tr><td>Solibacteres        </td><td>root     </td><td>0.01770921</td></tr>
	<tr><td>Sphingobacteriia    </td><td>root     </td><td>0.17481134</td></tr>
	<tr><td>Subgroup_6          </td><td>root     </td><td>0.02532282</td></tr>
	<tr><td>uncultured_bacterium</td><td>root     </td><td>0.04946056</td></tr>
</tbody>
</table>



    [1] "Most abundant Class"



    
![png](output_13_3.png)
    



```R
#Most abundant genus
#Tax glom at genus level
gen_bact <- tax_glom(bact, "Genus")
gen_bact
```


    phyloseq-class experiment-level object
    otu_table()   OTU Table:         [ 968 taxa and 2 samples ]
    tax_table()   Taxonomy Table:    [ 968 taxa by 7 taxonomic ranks ]



```R
#Heat map of top abundant genus in both samples
top50 <- prune_taxa(names(sort(taxa_sums(gen_bact),TRUE)[1:50]), gen_bact)
gen_heat_plot <- plot_heatmap(top50, taxa.label = "Genus", method = NULL, low = "#f9f8f9ff",  trans = log_trans(2),
                              high ="black", na.value = "white", 
                              sample.order = c("substrate", "root"),
                              taxa.order  = names(sort(taxa_sums(top50))))+ theme_light()
ggsave("gen_heat_plot.pdf", width=20, height=22, units="cm")
gen_heat_plot
```


    
![png](output_15_0.png)
    



```R
#Heat map of top abundant genus in both samples
top50 <- prune_taxa(names(sort(taxa_sums(gen_bact),TRUE)[1:50]), gen_bact)
gen_heat_plot <- plot_heatmap(top50, method = NULL, low = "#f9f8f9ff",  trans = log_trans(2),
                              high ="black", na.value = "white", 
                              sample.order = c("substrate", "root"),
                              taxa.order  = names(sort(taxa_sums(top50))))+ theme_light()
ggsave("gen_heat_plot_NL.pdf", width=20, height=22, units="cm")
gen_heat_plot
```


    
![png](output_16_0.png)
    



```R
#Most abundant genus
#Tax glom at genus level
gen_bact <- tax_glom(rel_bact, "Genus")
gen_bact
```


    phyloseq-class experiment-level object
    otu_table()   OTU Table:         [ 968 taxa and 2 samples ]
    tax_table()   Taxonomy Table:    [ 968 taxa by 7 taxonomic ranks ]



```R
options(scipen = 999)
gen_otu <- otu_table(gen_bact)
gen_tax <- tax_table(gen_bact)
gen_fung_tab <- cbind(gen_otu, gen_tax)
gen_fung_tab
write.table(gen_fung_tab, "gen_fung_tab.tsv", sep = "\t")
```


<table class="dataframe">
<caption>A matrix: 485 × 8 of type chr</caption>
<thead>
	<tr><th></th><th scope=col>substrate</th><th scope=col>root</th><th scope=col>Kingdom</th><th scope=col>Phylum</th><th scope=col>Class</th><th scope=col>Order</th><th scope=col>Family</th><th scope=col>Genus</th></tr>
</thead>
<tbody>
	<tr><th scope=row>19</th><td>0.0000871543241617933 </td><td>0.00489609765367375  </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Nitrosomonadaceae     </td><td>IS-44                                     </td></tr>
	<tr><th scope=row>75</th><td>0.0098222923330341    </td><td>0.000603628477850188 </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Methylophilaceae      </td><td>MM1                                       </td></tr>
	<tr><th scope=row>76</th><td>0.00163850129424171   </td><td>0.00049184542639645  </td><td>Bacteria</td><td>Verrucomicrobiota</td><td>Verrucomicrobiae   </td><td>Chthoniobacterales   </td><td>Chthoniobacteraceae   </td><td>LD29                                      </td></tr>
	<tr><th scope=row>95</th><td>0                     </td><td>0.0000223566102907477</td><td>Bacteria</td><td>Gemmatimonadota  </td><td>Longimicrobia      </td><td>Longimicrobiales     </td><td>Longimicrobiaceae     </td><td>Longimicrobium                            </td></tr>
	<tr><th scope=row>129</th><td>0.00154263153766374   </td><td>0.00435953900669581  </td><td>Bacteria</td><td>Gemmatimonadota  </td><td>Longimicrobia      </td><td>Longimicrobiales     </td><td>Longimicrobiaceae     </td><td>YC-ZSS-LKJ147                             </td></tr>
	<tr><th scope=row>233</th><td>0.000331186431814815  </td><td>0.00163203255122458  </td><td>Bacteria</td><td>Acidobacteriota  </td><td>Thermoanaerobaculia</td><td>Thermoanaerobaculales</td><td>Thermoanaerobaculaceae</td><td>Subgroup                                  </td></tr>
	<tr><th scope=row>340</th><td>0.00134217659209162   </td><td>0.00266043662459898  </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Comamonadaceae        </td><td>Kinneretia                                </td></tr>
	<tr><th scope=row>379</th><td>0.000252747540069201  </td><td>0.000704233224158553 </td><td>Bacteria</td><td>Acidobacteriota  </td><td>Vicinamibacteria   </td><td>Vicinamibacterales   </td><td>Vicinamibacteraceae   </td><td>Vicinamibacter                            </td></tr>
	<tr><th scope=row>417</th><td>0.000052292594497076  </td><td>0.000156496272035234 </td><td>Bacteria</td><td>Bdellovibrionota </td><td>Oligoflexia        </td><td>Silvanigrellales     </td><td>Silvanigrellaceae     </td><td>Silvanigrella                             </td></tr>
	<tr><th scope=row>463</th><td>0.000026146297248538  </td><td>0.000771303055030796 </td><td>Bacteria</td><td>Myxococcota      </td><td>Myxococcia         </td><td>Myxococcales         </td><td>Myxococcaceae         </td><td>Archangium                                </td></tr>
	<tr><th scope=row>499</th><td>0.218182135106633     </td><td>0.00568975731899529  </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Burkholderiaceae      </td><td>Burkholderia-Caballeronia-Paraburkholderia</td></tr>
	<tr><th scope=row>566</th><td>0.000104585188994152  </td><td>0.0016655674666607   </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Steroidobacterales   </td><td>Steroidobacteraceae   </td><td>Povalibacter                              </td></tr>
	<tr><th scope=row>587</th><td>0.0000871543241617933 </td><td>0.000245922713198225 </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Oxalobacteraceae      </td><td>Paraherbaspirillum                        </td></tr>
	<tr><th scope=row>645</th><td>0.000296324702150097  </td><td>0.000257101018343599 </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Comamonadaceae        </td><td>Leptothrix                                </td></tr>
	<tr><th scope=row>694</th><td>0.00000871543241617933</td><td>0.000558915257268693 </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Nitrosomonadaceae     </td><td>oc32                                      </td></tr>
	<tr><th scope=row>751</th><td>0.0000610080269132553 </td><td>0.000514202036687198 </td><td>Bacteria</td><td>Firmicutes       </td><td>Bacilli            </td><td>Bacillales           </td><td>Planococcaceae        </td><td>Filibacter                                </td></tr>
	<tr><th scope=row>765</th><td>0.000078438891745614  </td><td>0.000201209492616729 </td><td>Bacteria</td><td>Myxococcota      </td><td>Myxococcia         </td><td>Myxococcales         </td><td>Myxococcaceae         </td><td>KD3-10                                    </td></tr>
	<tr><th scope=row>779</th><td>0.00275407664351267   </td><td>0.00732178987021988  </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Steroidobacterales   </td><td>Steroidobacteraceae   </td><td>Steroidobacter                            </td></tr>
	<tr><th scope=row>825</th><td>0                     </td><td>0.0000223566102907477</td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Comamonadaceae        </td><td>Comamonas                                 </td></tr>
	<tr><th scope=row>871</th><td>0.000191739513155945  </td><td>0.00137493153288098  </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Comamonadaceae        </td><td>Caenimonas                                </td></tr>
	<tr><th scope=row>874</th><td>0.00779159658006432   </td><td>0.00985926513821974  </td><td>Bacteria</td><td>Firmicutes       </td><td>Bacilli            </td><td>Bacillales           </td><td>Bacillaceae           </td><td>Bacillus                                  </td></tr>
	<tr><th scope=row>919</th><td>0.00161235499699318   </td><td>0.00401301154718922  </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Oxalobacteraceae      </td><td>Noviherbaspirillum                        </td></tr>
	<tr><th scope=row>1066</th><td>0.000026146297248538  </td><td>0.000134139661744486 </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Xanthomonadales      </td><td>Xanthomonadaceae      </td><td>Luteimonas                                </td></tr>
	<tr><th scope=row>1081</th><td>0.0000958697565779726 </td><td>0.00198973831587655  </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Comamonadaceae        </td><td>Ideonella                                 </td></tr>
	<tr><th scope=row>1125</th><td>0.0000958697565779726 </td><td>0.000178852882325982 </td><td>Bacteria</td><td>Acidobacteriota  </td><td>Vicinamibacteria   </td><td>Vicinamibacterales   </td><td>Vicinamibacteraceae   </td><td>Luteitalea                                </td></tr>
	<tr><th scope=row>1215</th><td>0.00230958959028752   </td><td>0.00497434578969137  </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Comamonadaceae        </td><td>Rhizobacter                               </td></tr>
	<tr><th scope=row>1224</th><td>0.0113213467086169    </td><td>0.00220212611363865  </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Oxalobacteraceae      </td><td>Massilia                                  </td></tr>
	<tr><th scope=row>1269</th><td>0.00026146297248538   </td><td>0.00738885970109212  </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Nitrosomonadaceae     </td><td>MND1                                      </td></tr>
	<tr><th scope=row>1314</th><td>0.001263737700346     </td><td>0.0187795526442281   </td><td>Bacteria</td><td>Verrucomicrobiota</td><td>Verrucomicrobiae   </td><td>Opitutales           </td><td>Opitutaceae           </td><td>Opitutus                                  </td></tr>
	<tr><th scope=row>1370</th><td>0.000200454945572125  </td><td>0.0010284040733744   </td><td>Bacteria</td><td>Proteobacteria   </td><td>Gammaproteobacteria</td><td>Burkholderiales      </td><td>Nitrosomonadaceae     </td><td>mle1-7                                    </td></tr>
	<tr><th scope=row>⋮</th><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td><td>⋮</td></tr>
	<tr><th scope=row>51822</th><td>0.00000871543241617933</td><td>0.0000111783051453739</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Caulobacterales  </td><td>Caulobacteraceae  </td><td>PMMR1           </td></tr>
	<tr><th scope=row>51932</th><td>0.000366048161479532  </td><td>0.000100604746308365 </td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Caulobacterales  </td><td>Caulobacteraceae  </td><td>Asticcacaulis   </td></tr>
	<tr><th scope=row>52058</th><td>0.000113300621410331  </td><td>0.000167674577180608 </td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Acetobacterales  </td><td>Acetobacteraceae  </td><td>Rhodopila       </td></tr>
	<tr><th scope=row>52437</th><td>0                     </td><td>0.0000223566102907477</td><td>Bacteria</td><td>Cyanobacteria </td><td>Cyanobacteriia     </td><td>Cyanobacteriales </td><td>Nostocaceae       </td><td>Anabaena        </td></tr>
	<tr><th scope=row>52478</th><td>0.0000174308648323587 </td><td>0.0000111783051453739</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Rhizobiales      </td><td>Beijerinckiaceae  </td><td>Roseiarcus      </td></tr>
	<tr><th scope=row>52588</th><td>0.0000174308648323587 </td><td>0.0000111783051453739</td><td>Bacteria</td><td>Chloroflexi   </td><td>Chloroflexia       </td><td>Chloroflexales   </td><td>Herpetosiphonaceae</td><td>Herpetosiphon   </td></tr>
	<tr><th scope=row>52650</th><td>0                     </td><td>0.0000223566102907477</td><td>Bacteria</td><td>Chloroflexi   </td><td>Chloroflexia       </td><td>Chloroflexales   </td><td>Chloroflexaceae   </td><td>Chloronema      </td></tr>
	<tr><th scope=row>52793</th><td>0.00000871543241617933</td><td>0.0000111783051453739</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Rhodobacterales  </td><td>Rhodobacteraceae  </td><td>Falsirhodobacter</td></tr>
	<tr><th scope=row>52801</th><td>0                     </td><td>0.0000223566102907477</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Defluviicoccales </td><td>Defluviicoccaceae </td><td>Defluviicoccus  </td></tr>
	<tr><th scope=row>52823</th><td>0.00000871543241617933</td><td>0.0000111783051453739</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Sphingomonadales </td><td>Sphingomonadaceae </td><td>DSSF69          </td></tr>
	<tr><th scope=row>53151</th><td>0.0000435771620808966 </td><td>0.0000223566102907477</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Acetobacterales  </td><td>Acetobacteraceae  </td><td>Roseomonas      </td></tr>
	<tr><th scope=row>53276</th><td>0.000104585188994152  </td><td>0.0000111783051453739</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Caulobacterales  </td><td>Hyphomonadaceae   </td><td>Hirschia        </td></tr>
	<tr><th scope=row>53393</th><td>0.00013073148624269   </td><td>0.0000894264411629909</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Rhizobiales      </td><td>Devosiaceae       </td><td>Arsenicitalea   </td></tr>
	<tr><th scope=row>53642</th><td>0.0000174308648323587 </td><td>0.0000223566102907477</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Sphingomonadales </td><td>Sphingomonadaceae </td><td>Parablastomonas </td></tr>
	<tr><th scope=row>53692</th><td>0.00000871543241617933</td><td>0.0000111783051453739</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Rhizobiales      </td><td>Rhizobiaceae      </td><td>Ochrobactrum    </td></tr>
	<tr><th scope=row>53820</th><td>0.000226601242820663  </td><td>0.0000111783051453739</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Caedibacterales  </td><td>Caedibacteraceae  </td><td>Caedibacter     </td></tr>
	<tr><th scope=row>54317</th><td>0.000156877783491228  </td><td>0.000078248136017617 </td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Rhizobiales      </td><td>Xanthobacteraceae </td><td>Variibacter     </td></tr>
	<tr><th scope=row>54326</th><td>0.000174308648323587  </td><td>0.0000111783051453739</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Rhizobiales      </td><td>Rhizobiaceae      </td><td>[Rhizobium]     </td></tr>
	<tr><th scope=row>54340</th><td>0.000357332729063352  </td><td>0.0000335349154361216</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Sphingomonadales </td><td>Sphingomonadaceae </td><td>Stakelama       </td></tr>
	<tr><th scope=row>54410</th><td>0.0000871543241617933 </td><td>0.0000111783051453739</td><td>Bacteria</td><td>Proteobacteria</td><td>Gammaproteobacteria</td><td>Burkholderiales  </td><td>Oxalobacteraceae  </td><td>s3t2d-1089      </td></tr>
	<tr><th scope=row>54426</th><td>0.000113300621410331  </td><td>0.0000111783051453739</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Rhizobiales      </td><td>Rhizobiaceae      </td><td>Ciceribacter    </td></tr>
	<tr><th scope=row>54465</th><td>0.0000871543241617933 </td><td>0                    </td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Rhizobiales      </td><td>Kaistiaceae       </td><td>Kaistia         </td></tr>
	<tr><th scope=row>54655</th><td>0.000418340755976608  </td><td>0.0000335349154361216</td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Rhizobiales      </td><td>Xanthobacteraceae </td><td>Oligotropha     </td></tr>
	<tr><th scope=row>54675</th><td>0.000026146297248538  </td><td>0                    </td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Holosporales     </td><td>Holosporaceae     </td><td>Candidatus      </td></tr>
	<tr><th scope=row>54761</th><td>0.0000348617296647173 </td><td>0                    </td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Acetobacterales  </td><td>Acetobacteraceae  </td><td>Roseococcus     </td></tr>
	<tr><th scope=row>54770</th><td>0.0000174308648323587 </td><td>0                    </td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Rhizobiales      </td><td>Xanthobacteraceae </td><td>Rhodopseudomonas</td></tr>
	<tr><th scope=row>54948</th><td>0.0000871543241617933 </td><td>0.0000111783051453739</td><td>Bacteria</td><td>Proteobacteria</td><td>Gammaproteobacteria</td><td>Burkholderiales  </td><td>Alcaligenaceae    </td><td>Ampullimonas    </td></tr>
	<tr><th scope=row>55028</th><td>0.0000174308648323587 </td><td>0                    </td><td>Bacteria</td><td>Proteobacteria</td><td>Alphaproteobacteria</td><td>Tistrellales     </td><td>Geminicoccaceae   </td><td>Candidatus      </td></tr>
	<tr><th scope=row>56529</th><td>0                     </td><td>0.0000894264411629909</td><td>Bacteria</td><td>Deinococcota  </td><td>Deinococci         </td><td>Deinococcales    </td><td>Trueperaceae      </td><td>Truepera        </td></tr>
	<tr><th scope=row>56887</th><td>0                     </td><td>0.0000335349154361216</td><td>Archaea </td><td>Crenarchaeota </td><td>Nitrososphaeria    </td><td>Nitrososphaerales</td><td>Nitrososphaeraceae</td><td>Candidatus      </td></tr>
</tbody>
</table>




```R
#Venn diagram
otu_list <- otu_table(gen_bact)

substrate <- otu_list[, 1]
substrate[substrate==0] <- NA
substrate2 <- substrate[complete.cases(substrate),]
write.table(rownames(substrate2) ,"bac_substrate_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

root <- otu_list[, 2]
root[root==0] <- NA
root2 <- root[complete.cases(root),]
write.table(rownames(root2),"bac_root_list.txt", quote = FALSE,  row.names= FALSE, col.names = FALSE)


head(substrate2)
head(root2)

#The files are used as input list to make Venn diagrams with the online tool:
#http://bioinformatics.psb.ugent.be/webtools/Venn/

```


<table class="dataframe">
<caption>A otu_table: 6 × 1 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>substrate</th></tr>
</thead>
<tbody>
	<tr><th scope=row>19</th><td>0.00008715432</td></tr>
	<tr><th scope=row>75</th><td>0.00982229233</td></tr>
	<tr><th scope=row>76</th><td>0.00163850129</td></tr>
	<tr><th scope=row>129</th><td>0.00154263154</td></tr>
	<tr><th scope=row>233</th><td>0.00033118643</td></tr>
	<tr><th scope=row>340</th><td>0.00134217659</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A otu_table: 6 × 1 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>root</th></tr>
</thead>
<tbody>
	<tr><th scope=row>19</th><td>0.00489609765</td></tr>
	<tr><th scope=row>75</th><td>0.00060362848</td></tr>
	<tr><th scope=row>76</th><td>0.00049184543</td></tr>
	<tr><th scope=row>95</th><td>0.00002235661</td></tr>
	<tr><th scope=row>129</th><td>0.00435953901</td></tr>
	<tr><th scope=row>233</th><td>0.00163203255</td></tr>
</tbody>
</table>


