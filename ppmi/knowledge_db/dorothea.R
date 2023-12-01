## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load, message=FALSE------------------------------------------------------
library(dorothea)
library(decoupleR)
library(ggplot2)
library(dplyr)

## ----model--------------------------------------------------------------------
net <- decoupleR::get_dorothea(levels = c('A', 'B', 'C', 'D'))
head(net)

## ----n_genes------------------------------------------------------------------
n_genes <- net %>%
  group_by(source) %>%
  summarize(n = n())

ggplot(data=n_genes, aes(x=n)) +
  geom_density() +
  theme(text = element_text(size=12)) +
  xlab('Number of target genes') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")


## ----n_edges------------------------------------------------------------------
n_edges <- net %>%
  group_by(confidence) %>%
  summarize(n = n())

ggplot(data=n_edges, aes(x=confidence, y=log10(n), color=confidence, fill=confidence)) +
  geom_bar(stat="identity") +
  theme(text = element_text(size=12)) +
  xlab('Confidence') +
  ylab('log10(Number of edges)') +
  theme_bw() +
  theme(legend.position = "none")


## ----prop---------------------------------------------------------------------
prop <- net %>%
  mutate(mor = case_when(mor < 0 ~ -1, mor > 0 ~ 1)) %>%
  group_by(source, mor) %>%
  summarize(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(mor == 1)

ggplot(data=prop, aes(x=freq)) +
  geom_density() +
  theme(text = element_text(size=12)) +
  xlab('% of positive edges') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")


## ----collectri----------------------------------------------------------------
net <- decoupleR::get_collectri(split_complexes = FALSE)
head(net)

## ----n_genes_collectri--------------------------------------------------------
n_genes <- net %>%
  group_by(source) %>%
  summarize(n = n())

ggplot(data=n_genes, aes(x=n)) +
  geom_density() +
  theme(text = element_text(size=12)) +
  xlab('Number of target genes') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")


## ----prop_collectri-----------------------------------------------------------
prop <- net %>%
  mutate(mor = case_when(mor < 0 ~ -1, mor > 0 ~ 1)) %>%
  group_by(source, mor) %>%
  summarize(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(mor == 1)

ggplot(data=prop, aes(x=freq)) +
  geom_density() +
  theme(text = element_text(size=12)) +
  xlab('% of positive edges') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")

