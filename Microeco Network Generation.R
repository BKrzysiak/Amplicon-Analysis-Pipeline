#Load Necessary Packages
library(microeco); packageVersion("microeco")
library(file2meco); packageVersion("file2meco")
library(magrittr); packageVersion("magrittr")
library(ggh4x); packageVersion("ggh4x")
library(rgexf); packageVersion("rgexf")
library(paletteer); packageVersion("paletteer")
library(microViz); packageVersion("microViz")
library(igraph); packageVersion("igraph")
library(pulsar); packageVersion("pulsar")
library(CoDaSeq); packageVersion("CoDaSeq")
library(tibble); packageVersion("tibble")
library(tidyverse); packageVersion("tidyverse")

#Phyloseq to Microeco Handoff
meco<-phyloseq2meco(psr)
meco$tidy_dataset()
