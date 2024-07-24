---
title: "Network_Analysis_by_netcomi"
output: 
  html_document: 
    keep_md: true
date: "2024-07-24"
---



  Networks are used all around: Political science, games, ecology, metabolics, etc. Co-occurrence networks assume that when two subjects, such as animals, people, or virtual connections, are frequently appearing together that they are associating in some way. For example, two people frequently messaging on social media are likely friends and an individual who is often messaged or retweeted is important to the structure of the community they've surrounded themselves with, just like an influencer.
  We can do similar things with microbial communities. By looking at Amplicon sequencing data from environmental samples. we can see which species or taxa often occur together. This can imply potential interactions (such as predation or symbiosis) or a reliance on similar conditions for survival. One of the most important parts about networks is that the statistics are all dependent, and it is this dependency that interests us.
  To that end, we use microeco (https://github.com/ChiLiubio/microeco) to explore networks as it allows for the conversion of phylseq objects to microeco objects, the generation of networks, and various statistical tests to explore how network metrics correlate to environmental factors. Though microeco does allow for visualization, we use Gephi which is a program meant for the visualization of networks. Instructions for download can be found online.
  
  Load the following packages:

```r
library(readxl)
library(openxlsx2)
```

```
## 
## Attaching package: 'openxlsx2'
```

```
## The following object is masked from 'package:readxl':
## 
##     read_xlsx
```

```r
library(BiocManager)
```

```
## Bioconductor version '3.18' is out-of-date; the current release version '3.19'
##   is available with R version '4.4'; see https://bioconductor.org/install
```

```r
library(phyloseq)
library(dada2)
```

```
## Loading required package: Rcpp
```

```r
library(microeco)
library(file2meco)
library(magrittr)
library(ggh4x)
```

```
## Loading required package: ggplot2
```

```
## 
## Attaching package: 'ggh4x'
```

```
## The following object is masked from 'package:ggplot2':
## 
##     guide_axis_logticks
```

```r
library(rgexf)
library(paletteer)
library(igraph)
```

```
## 
## Attaching package: 'igraph'
```

```
## The following objects are masked from 'package:stats':
## 
##     decompose, spectrum
```

```
## The following object is masked from 'package:base':
## 
##     union
```

```r
library(pulsar)
library(CoDaSeq)
```

```
## Loading required package: ALDEx2
```

```
## Loading required package: zCompositions
```

```
## Warning: package 'zCompositions' was built under R version 4.3.3
```

```
## Loading required package: MASS
```

```
## Loading required package: NADA
```

```
## Loading required package: survival
```

```
## 
## Attaching package: 'NADA'
```

```
## The following object is masked from 'package:stats':
## 
##     cor
```

```
## Loading required package: truncnorm
```

```
## Loading required package: lattice
```

```
## Loading required package: latticeExtra
```

```
## 
## Attaching package: 'latticeExtra'
```

```
## The following object is masked from 'package:ggplot2':
## 
##     layer
```

```
## Loading required package: car
```

```
## Loading required package: carData
```

Step one is to convert our phyloseq object into a microeco object. This will be using functions from the file2meco package.


```r
#Loading Phyloseq object from library
data(GlobalPatterns)
meco<-phyloseq2meco(GlobalPatterns)
```

```
## 228 taxa with 0 abundance are removed from the otu_table ...
```

```r
#when converting to a meco object, it will "clean up" the OTU and Taxa table by removing those that have 0 abundances across all samples. In truth, this causes issues down the line so we have to tidy the dataset before moving on. It is also important to run this function anytime you modify the microeco object by filtering, subgrouping, etc.
meco$tidy_dataset()
```
  Meco objects are a little special. If you inspect it directly you will be able to see various objects such as cal_abund, otu_tables, sample_table, and so much more. Meco objects stores the taxa table, otu table, and metadata file from phyloseq but it also stores the functions which the package calls on (such as cal_abund). This can make troubleshooting or modifying steps in the process a bit difficult unless you are familiar with messing around with functions in R.
  Now that we have this object, we can begin exploring networks!


```r
#Dataset should be the microeco object, cor_methods are elaborated on the microeco page. Taxa level can be specified to any level. Lets look at this network from a Family level.
network<-trans_network$new(dataset = meco, cor_method = "spearman",filter_thres = .0001,
                          taxa_level = "Family")
```

```
## After filtering, 232 features are remained ...
```

```
## The correlation result list is stored in object$res_cor_p ...
```

```r
#With our network created, we can now calculate it. In this case we have our significance set to 95% confidence (p=.5) and have a correlation cutoff of .6.Additionally, we use the add_taxa_name function as it facilatates exporting this important information to gephi.
network$cal_network(COR_p_thres = .05,COR_cut=.6,add_taxa_name = c("Kingdom",
  "Phylum","Class","Order","Family","Genus","Species"))
```

```
## ---------------- 2024-07-24 16:35:27.159755 : Start ----------------
```

```
## Perform p value adjustment with fdr method ...
```

```
## ---------------- 2024-07-24 16:35:27.500901 : Finish ----------------
```

```
## Skip adding taxonomy: Genus to node as it is not in colnames of tax_table ...
```

```
## Skip adding taxonomy: Species to node as it is not in colnames of tax_table ...
```

```
## The result network is stored in object$res_network ...
```

```r
#Now that is has calculated which members of the network are significant, we calculate how our modules are shaped. Again, the microeco page has more information on the different methods.
network$cal_module(method = "cluster_fast_greedy")
```

```
## Use cluster_fast_greedy function to partition modules ...
```

```
## Totally, 4 modules are idenfified ...
```

```
## Modules are assigned in network with attribute name -- module ...
```

```r
#Finally, we export the object as a gephi file.
network$save_network(filepath = "network.gexf")
```

  Great, we have a network! The fun part of visualizing it occurs in Gephi. First, we have to open gephi and open our new Gephi file which should be saved in the same place as your R repository:

![Opening Gephi File](C:/Users/eckha/OneDrive/Documents/Amplicon-Analysis-Pipeline/Gephi Reference Images/Screenshot 2024-07-24 154218.png)


  When you open it, a screen like this should appear. Simply make sure that graph type is "undirected" and hit ok.

![Finishing importing](C:/Users/eckha/OneDrive/Documents/Amplicon-Analysis-Pipeline/Gephi Reference Images/Screenshot 2024-07-24 154241.png)
  
  
  Now your network is in! And looking like a weird orange glob. Since this isn't the most imformative, lets use Gephi's tool to improve the appearance of this network. First, lets work in the appearance tab. From here we can modify the color, size, labels, and text of nodes and edges. Lets change the color of our nodes so that they represent their Family.
  To do so, we click on the paint pallete symbol and make sure node is selected. Click on the partition tab and select Family. Hit apply.
  
  ![Partitioning Color](C:/Users/eckha/OneDrive/Documents/Amplicon-Analysis-Pipeline/Gephi Reference Images/Screenshot 2024-07-24 155239.png)
  
  
  Well, it seems that there aren't enough colors in the palette. Lets generate a better palette. We can do so by clicking Palette... above apply. From here, we click on generate. Turn of "limit number of colors" and hit generate. You are also able to change the range of colors through the "presets" dropdown.
  
  ![Clicking on palette](C:/Users/eckha/OneDrive/Documents/Amplicon-Analysis-Pipeline/Gephi Reference Images/Screenshot 2024-07-24 155736.png)
  ![Generating a palette](C:/Users/eckha/OneDrive/Documents/Amplicon-Analysis-Pipeline/Gephi Reference Images/Screenshot 2024-07-24 155814.png)
  

  Now we can hit Ok and when we check our palette options you should see the new palette. Apply it. Now that our blob is colorful, we can modify our nodes further. We can make the size of each node correlate to some helpful information, such as the relative abundance of our taxa. To do so click on the size tab (stacked circles). From here make sure you are on "Nodes" and click on the "Ranking" tab. click on RelativeAbundance in the dropdown menu. From here you'll get the option to change the minimum and maximum size. change these to what you like and hit apply.
  
  ![Node Size](C:/Users/eckha/OneDrive/Documents/Amplicon-Analysis-Pipeline/Gephi Reference Images/Screenshot 2024-07-24 161316.png)
  
  
  Now we can de-blob our network. The Layout pane allows us to change how our network is visually structured. There are a lot of options ranging from simple displays to force directed algorithms such as the Yifan Hu layout. Lets try a simple display first.
  click on the "choose a layout dropdown" and select "Fruchterman Reingold". Once you click on it, you will be able to modify the graph area, the gravity, and the speed. Mess with these as you like! If you ever need to "reset it", select "Random Layout" and hit run. Once you have a setting you like, run the Fruchterman Reingold layout.
  
  ![Layout](C:/Users/eckha/OneDrive/Documents/Amplicon-Analysis-Pipeline/Gephi Reference Images/Screenshot 2024-07-24 161823.png)
  ![Frucht Layout](C:/Users/eckha/OneDrive/Documents/Amplicon-Analysis-Pipeline/Gephi Reference Images/Screenshot 2024-07-24 162103.png)
  
  
  Now that the layout is built, we can continue modifying to make it visually more appealing. In this case, I would change the size of the nodes so that they are more visible and maybe even reconsider some of the abhorrent color choices I made...
  Regardless, this is the basics of getting a network visualized in Gephi. there are plenty of other things you can do, and Gephi provides some tutorial on its website (https://gephi.org/users/).
  
  
  
