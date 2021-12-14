# Author ------------------------------------------------------------------

# @author:   Axel Kunstner
# @date:    2021-12-13

# libraries ---------------------------------------------------------------

# for data handling
library(tidyverse)
library(phyloseq)
library(ape)

# for estimates and analysis
library(DivNet)
library(estimatr)
library(selbal)
library(ppcor)
library(ALDEx2)
library(corncob)

# for multi panel plots
library(cowplot)
library(gridExtra)
library(ggpubr)

# Set colors and seed -----------------------------------------------------

# set seed
seedID <- 246
set.seed(seedID)

# general color schemes
colv2 <- c("steelblue", "orange")
colv4 <- c("steelblue1", "steelblue4", "orange1", "orange4")

# Colors for taxonomic abundance plots
taxNumberPlot <- 19
cols <- c(pals::tableau20(taxNumberPlot+1), "grey45")

# Retrieve final data set -------------------------------------------------

ps <- readRDS(file = "data/phyloseq.ps.RDS")

ps_ra  <- transform_sample_counts(ps, function(x){x / sum(x)})
ps_clr <- microbiome::transform(ps, "clr")

metadata <- sample_data(ps) %>% data.frame()

ps_kegg <- readRDS(file = "data/phyloseq.KEGG.RDS")

# Figure 1 ----------------------------------------------------------------

# alpha diversity
# Species richness using breakaway estimator
ba <- breakaway::breakaway(input_data = ps )

ba.richness <- summary(ba) %>%
    add_column("SampleNames" = ps %>% otu_table %>% sample_names) %>%
    add_column("Depth"       = ps %>% sample_sums)
# check whether richness and sequencing depth are associated
lm1 <- estimatr::lm_robust( estimate ~ Depth, se_type = "HC1", data = ba.richness )
summary(lm1) # no sign. association with sequencing depth

ba.richness <- merge(x = metadata,
                     y = ba.richness %>% dplyr::select(sample_names, Richness = estimate, Richness.error = error),
                     by.x = 'id', by.y = 'sample_names')

fig_1a <- ba.richness %>%
    dplyr::mutate(Var = factor(paste(Var, Strain))) %>%
    ggplot(., aes(x = Var, y = Richness, color = Var)) +
    geom_violin(color = 'black' ) +
    ggforce::geom_sina() +
    ylim(20, 100) +
    xlab("") + ylab("breakaway estimate of species richness") +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          # axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12)) +
    ggpubr::stat_compare_means(label.x = 1.5, label.y = 100) +
    ggpubr::stat_compare_means(comparisons = mycomp, label.y = c(80, 80, 80, 85, 90),
                               method = "wilcox.test") +
    scale_color_manual(values=c("grey45", "grey45", colv4[c(1, 1, 3, 3)])) +
    ggtitle("")
fig_1a

# beta diversity -> Aitchison distance
otu.table_clr <- otu_table(ps_clr)
ps_clr_dist <- dist(otu.table_clr, method="euclidean")
ps_clr_ord <- phyloseq::ordinate(ps_clr, "RDA", distance = "euclidean")
fig_1b <- plot_ordination( physeq = ps_clr, ordination = ps_clr_ord, color='Var', shape = "Strain")  +
    scale_color_manual(values=c("grey45", colv4[c(1,3)])) +
    geom_point(size=4) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12)) +
    ggtitle("")
fig_1b

# PERMANOVA on Aitchison distance
vegan::adonis(formula = ps_clr_dist ~ Strain + Week + Strain:Diet2,
               data = metadata,
               permutations = 99999)

# Supplementary Fig 1 -----------------------------------------------------

#
# Alpha diversity was estimated using DivNet;
# precalculated estimates are loaded below
#

# dv_ind <- DivNet::divnet(ps, base = "ASV1", ncores = 10,
#                          network = "diagonal",
#                          tuning = 'fast')
# saveRDS(object = dv_ind, file = "data/divnet.ind.ASV.RDS")

dv_ind <- readRDS(file = "data/divnet.ind.ASV.RDS")

combined_shannon <- ps %>% sample_data %>% data.frame %>%
    mutate(sample_names = rownames(.)) %>%
    left_join(dv_ind$shannon %>% summary,
              by = "sample_names") %>%
    dplyr::filter( !is.na(estimate) )

mycomp <- list(
    c("w0 CD ALR", "w0 CD BPL"),
    c("w8 CD ALR", "w8 CD BPL"),
    c("w8 HFD ALR", "w8 HFD BPL"),
    c("w0 CD ALR", "w8 CD ALR"),
    c("w0 CD BPL", "w8 CD BPL")
)

suppl_fig1 <- combined_shannon %>%
    dplyr::mutate(Var = factor(paste(Var, Strain))) %>%
    ggplot(., aes(x = Var, y = estimate, color = Var)) +
    geom_violin(color = 'black' ) +
    ggforce::geom_sina() +
    ylim(0, 6) +
    xlab("") + ylab("DivNet estimate of Shannon") +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text.y = element_text(size=12)) +
    ggpubr::stat_compare_means(label.x = 1.5, label.y = 5.5) +
    ggpubr::stat_compare_means(comparisons = mycomp, label.y = c(4.25, 4.25, 4.25, 4.5, 4.75),
                               method = "wilcox.test") +
    scale_color_manual(values=c("grey45", "grey45", colv4[c(1, 1, 3, 3)]))
suppl_fig1

# Figure 2 ----------------------------------------------------------------

# Phylum and Genus abundances (proportions)

# Phyla
fig_2a <- ps %>%
    tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} ) %>%
    psmelt() %>%                                         # Melt to long format
    #filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
    arrange(Phylum) %>%
    dplyr::select( Phylum, Strain, Var, Abundance) %>%
    dplyr::group_by(Phylum, Strain, Var) %>%
    dplyr::summarise(mean = mean(Abundance), groups = c(Phylum)) %>%
    unique() %>%
    dplyr::mutate(Phylum = gsub(pattern = "_", replacement = " ", x = Phylum)) %>%
    ggplot(., aes(x = Var, y = mean, fill = Phylum)) +
    facet_wrap(Strain~., scales = "free_x") +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cols) +
    # Remove x axis title
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
    ylab("Relative Abundance \n") +
    ggtitle("") +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid")) +
    ylim(0,1) +
    guides(fill = guide_legend(title = "Phylum", ncol = 1))
fig_2a

# Genera
means_genus <- ps %>%
    tax_glom(taxrank = "Genus") %>% # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} ) %>%
    psmelt() %>% # Melt to long format
    dplyr::arrange(Genus) %>%
    dplyr::select( Genus, Strain, Var, Abundance) %>%
    dplyr::group_by(Genus, Strain, Var) %>%
    dplyr::summarise(mean = mean(Abundance), groups = c(Genus)) %>%
    unique() %>%
    dplyr::mutate(Genus = gsub(pattern = "_", replacement = " ", x = Genus))
# filter(mean >= 0.02) # Filter out low abundance taxa

# summarise lowly abundant genera
means_genus$Genus[means_genus$mean < 0.02] <- "Other"
means_genus <- means_genus %>%
    group_by(Genus, Strain, Var) %>%
    summarise(mean = sum(mean), groups = c(Genus)) %>%
    unique()

fig_2b <- ggplot(means_genus, aes(x = Var, y = mean, fill = Genus)) +
    facet_wrap(Strain~., scales = "free_x") +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cols) +
    # Remove x axis title
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
    ylab("Relative Abundance (Genus abundance > 2%)") +
    ggtitle("") +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(colour="white", fill="grey95", size=1.5, linetype="solid")) +
    ylim(0,1) +
    guides(fill = guide_legend(title = "Genus", ncol = 1))
fig_2b

# Figure 3 ----------------------------------------------------------------

# selbal

# Week 0
TAXRANK <- "Phylum"
selbal_dat <- ps %>%
    subset_samples(physeq = ., Week == "w0") %>%
    tax_glom(physeq = ., taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    otu_table() %>% t() %>% data.frame()
colnames(selbal_dat) <- ps %>%
    subset_samples(physeq = ., Week == "w0") %>%
    tax_glom(taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    tax_table() %>% data.frame() %>% .[, TAXRANK]
y <- ps %>%
    subset_samples(physeq = ., Week == "w0") %>%
    sample_data() %>% .$Strain
z <- ps %>%
    subset_samples(physeq = ., Week == "w0") %>%
    sample_data() %>% data.frame()
selbal <- selbal.cv(x = as.matrix(selbal_dat), y = y, maxV = 12,
                       n.fold = 5, n.iter = 10,
                       covar = data.frame(Sex = z$Sex),
                       zero.rep = "one",
                       logit.acc = 'AUC')
fig_3a <- selbal$global.plot

TAXRANK <- "Genus"
selbal_dat <- ps %>%
    subset_samples(physeq = ., Week == "w0") %>%
    tax_glom(physeq = ., taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    otu_table() %>% t() %>% data.frame()
colnames(selbal_dat) <- ps %>%
    subset_samples(physeq = ., Week == "w0") %>%
    tax_glom(taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    tax_table() %>% data.frame() %>% .[, TAXRANK]
y <- ps %>%
    subset_samples(physeq = ., Week == "w0") %>%
    sample_data() %>% .$Strain
z <- ps %>%
    subset_samples(physeq = ., Week == "w0") %>%
    sample_data() %>% data.frame()
selbal <- selbal.cv(x = as.matrix(selbal_dat), y = y, maxV = 12,
                       n.fold = 5, n.iter = 10,
                       covar = data.frame(Sex = z$Sex),
                       zero.rep = "one",
                       logit.acc = 'AUC')
fig_3b <- selbal$global.plot

# CD week 8
TAXRANK <- "Phylum"
selbal_dat <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "CD") %>%
    tax_glom(physeq = ., taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    otu_table() %>% t() %>% data.frame()
colnames(selbal_dat) <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "CD") %>%
    tax_glom(taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    tax_table() %>% data.frame() %>% .[, TAXRANK]
y <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "CD") %>%
    sample_data() %>% .$Strain
z <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "CD") %>%
    sample_data() %>% data.frame()
selbal <- selbal.cv(x = as.matrix(selbal_dat), y = y, maxV = 12,
                           n.fold = 5, n.iter = 10,
                           covar = data.frame(Sex = z$Sex),
                           zero.rep = "one",
                           logit.acc = 'AUC')
fig_3c <- selbal$global.plot
grid.draw(fig_3c)

TAXRANK <- "Genus"
selbal_dat <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "CD") %>%
    tax_glom(physeq = ., taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    otu_table() %>% t() %>% data.frame()
colnames(selbal_dat) <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "CD") %>%
    tax_glom(taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    tax_table() %>% data.frame() %>% .[, TAXRANK]
y <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "CD") %>%
    sample_data() %>% .$Strain
z <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "CD") %>%
    sample_data() %>% data.frame()
selbal <- selbal.cv(x = as.matrix(selbal_dat), y = y, maxV = 12,
                    n.fold = 5, n.iter = 10,
                    covar = data.frame(Sex = z$Sex),
                    zero.rep = "one",
                    logit.acc = 'AUC')
fig_3d <- selbal$global.plot
grid.draw(fig_3d)

# HFD week 8
TAXRANK <- "Phylum"
selbal_dat <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "HFD") %>%
    tax_glom(physeq = ., taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    otu_table() %>% t() %>% data.frame()
colnames(selbal_dat) <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "HFD") %>%
    tax_glom(taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    tax_table() %>% data.frame() %>% .[, TAXRANK]
y <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "HFD") %>%
    sample_data() %>% .$Strain
z <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "HFD") %>%
    sample_data() %>% data.frame()
selbal <- selbal.cv(x = as.matrix(selbal_dat), y = y, maxV = 12,
                    n.fold = 5, n.iter = 10,
                    covar = data.frame(Sex = z$Sex),
                    zero.rep = "one",
                    logit.acc = 'AUC')
fig_3e <- selbal$global.plot
grid.draw(fig_3e)

TAXRANK <- "Genus"
selbal_dat <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "HFD") %>%
    tax_glom(physeq = ., taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    otu_table() %>% t() %>% data.frame()
colnames(selbal_dat) <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "HFD") %>%
    tax_glom(taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    tax_table() %>% data.frame() %>% .[, TAXRANK]
y <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "HFD") %>%
    sample_data() %>% .$Strain
z <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == "HFD") %>%
    sample_data() %>% data.frame()
selbal <- selbal.cv(x = as.matrix(selbal_dat), y = y, maxV = 12,
                    n.fold = 5, n.iter = 10,
                    covar = data.frame(Sex = z$Sex),
                    zero.rep = "one",
                    logit.acc = 'AUC')
fig_3f <- selbal$global.plot
grid.draw(fig_3f)

# B6 - ALR week 8
TAXRANK <- "Phylum"
selbal_dat <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "ALR") %>%
    tax_glom(physeq = ., taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    otu_table() %>% t() %>% data.frame()
colnames(selbal_dat) <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "ALR") %>%
    tax_glom(taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    tax_table() %>% data.frame() %>% .[, TAXRANK]
y <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "ALR") %>%
    sample_data() %>% .$Diet
z <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "ALR") %>%
    sample_data() %>% data.frame()
selbal <- selbal.cv(x = as.matrix(selbal_dat), y = y, maxV = 12,
                    n.fold = 5, n.iter = 10,
                    covar = data.frame(Sex = z$Sex),
                    zero.rep = "one",
                    logit.acc = 'AUC')
fig_3g <- selbal$global.plot
grid.draw(fig_3g)

TAXRANK <- "Genus"
selbal_dat <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "ALR") %>%
    tax_glom(physeq = ., taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    otu_table() %>% t() %>% data.frame()
colnames(selbal_dat) <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "ALR") %>%
    tax_glom(taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    tax_table() %>% data.frame() %>% .[, TAXRANK]
y <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "ALR") %>%
    sample_data() %>% .$Diet
z <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "ALR") %>%
    sample_data() %>% data.frame()
selbal <- selbal.cv(x = as.matrix(selbal_dat), y = y, maxV = 12,
                    n.fold = 5, n.iter = 10,
                    covar = data.frame(Sex = z$Sex),
                    zero.rep = "one",
                    logit.acc = 'AUC')
fig_3h <- selbal$global.plot
grid.draw(fig_3h)

# B6 - BPL week 8
TAXRANK <- "Phylum"
selbal_dat <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "BPL") %>%
    tax_glom(physeq = ., taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    otu_table() %>% t() %>% data.frame()
colnames(selbal_dat) <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "BPL") %>%
    tax_glom(taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    tax_table() %>% data.frame() %>% .[, TAXRANK]
y <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "BPL") %>%
    sample_data() %>% .$Diet
z <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "BPL") %>%
    sample_data() %>% data.frame()
selbal <- selbal.cv(x = as.matrix(selbal_dat), y = y, maxV = 12,
                    n.fold = 5, n.iter = 10,
                    covar = data.frame(Sex = z$Sex),
                    zero.rep = "one",
                    logit.acc = 'AUC')
fig_3i <- selbal$global.plot
grid.draw(fig_3i)

TAXRANK <- "Genus"
selbal_dat <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "BPL") %>%
    tax_glom(physeq = ., taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    otu_table() %>% t() %>% data.frame()
colnames(selbal_dat) <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "BPL") %>%
    tax_glom(taxrank = TAXRANK) %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/4) %>%
    tax_table() %>% data.frame() %>% .[, TAXRANK]
y <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "BPL") %>%
    sample_data() %>% .$Diet
z <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == "BPL") %>%
    sample_data() %>% data.frame()
selbal <- selbal.cv(x = as.matrix(selbal_dat), y = y, maxV = 12,
                    n.fold = 5, n.iter = 10,
                    covar = data.frame(Sex = z$Sex),
                    zero.rep = "one",
                    logit.acc = 'AUC')
fig_3j <- selbal$global.plot
grid.draw(fig_3j)

# Figure 4 ----------------------------------------------------------------

#
# Phylum
#

ps_sub <- microbiome::aggregate_taxa(x = ps_clr %>%
                                         subset_samples(., !is.na(Blood_glucose_0_min)),
                                     level = "Phylum")
corr_df <- otu_table(ps_sub) %>%
    data.frame() %>% t()

corr_df <- cbind(
    sample_data(ps_sub) %>%
        data.frame() %>%
        dplyr::select(Strain, Diet2, Sex_num = Sex, Weight,
                      Blood_glucose_0_min,
                      Blood_glucose_45_min) %>%
        dplyr::mutate(Sex_num = as.numeric(Sex_num) ) %>%
        dplyr::mutate(Var = factor(paste(Strain, Diet2))),
    corr_df)

ppcor::pcor.test(x = corr_df$Weight, y = corr_df$Blood_glucose_0_min, z = corr_df[ , c("Sex_num")] )

cor_phyla <- data.frame(Time = NULL, Tax = NULL, Var = NULL, Correlation = NULL, p_val = NULL, Method = NULL)

for( bc in 4:6 ) {
    for( i in 8:ncol(corr_df)) {
        for( ii in levels(corr_df$Var) ) {
            tmp <- corr_df %>%
                dplyr::filter(Var == ii)
            res <- ppcor::pcor.test(x = tmp[, bc],
                                    y = tmp[, i],
                                    z = tmp[ , c("Sex_num")], method = "spearman" )
            cor_phyla <- rbind(cor_phyla,
                               data.frame(Time = colnames(tmp)[bc],
                                          Tax = colnames(tmp)[i],
                                          Var = ii,
                                          Correlation = res$estimate,
                                          p_val = res$p.value,
                                          Method = res$Method))
        }

    }
}

cor_phyla$Correlation[ is.nan(cor_phyla$Correlation) ] <- 0
cor_phyla$p_val[ is.nan(cor_phyla$p_val) ] <- 1
cor_phyla$p_clean <- cor_phyla$p_val
cor_phyla$p_clean <- format( round(cor_phyla$p_clean, 4), scientific = F ) # avoid exponentials
cor_phyla$p_clean[cor_phyla$p_clean >= 0.05] <- ""

fig_4a <- cor_phyla %>%
    dplyr::mutate(Time = gsub(pattern = "Blood_glucose_", replacement = "", Time)) %>%
    dplyr::mutate(Time = gsub(pattern = "_", replacement = "", Time)) %>%
    dplyr::mutate(Time = gsub(pattern = "ipGTTAUC", replacement = "AUC ipGTT", Time)) %>%
    dplyr::mutate(Time = factor(x = Time,
                                levels = c('0min', '45min', 'Weight') ) ) %>%
    ggplot(data = ., aes(x=Var, y=factor(Tax, levels = rev(levels(factor(Tax)))), fill=Correlation)) +
    geom_tile() +
    geom_text(aes(Var, Tax,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    facet_wrap(Time~., scales = "free_x", ncol = 3) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig_4a


#
# Genus
#

ps_sub <- microbiome::aggregate_taxa(x = ps_clr %>%
                                         subset_samples(., !is.na(Blood_glucose_0_min)),
                                     level = "Genus") %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 1/6)
corr_df <- otu_table(ps_sub) %>%
    data.frame() %>% t()

corr_df <- cbind(
    sample_data(ps_sub) %>%
        data.frame() %>%
        dplyr::select(Strain, Diet2, Sex_num = Sex, Weight,
                      Blood_glucose_0_min,
                      Blood_glucose_45_min) %>%
        dplyr::mutate(Sex_num = as.numeric(Sex_num) ) %>%
        dplyr::mutate(Var = factor(paste(Strain, Diet2))),
    corr_df)

cor_genus <- data.frame(Time = NULL, Tax = NULL, Var = NULL, Correlation = NULL, p_val = NULL, Method = NULL)

for( bc in 4:6 ) {
    for( i in 8:ncol(corr_df)) {
        for( ii in levels(corr_df$Var) ) {
            tmp <- corr_df %>%
                dplyr::filter(Var == ii)
            res <- ppcor::pcor.test(x = tmp[, bc],
                                    y = tmp[, i],
                                    z = tmp[ , c("Sex_num")], method = "spearman" )
            cor_genus <- rbind(cor_genus,
                               data.frame(Time = colnames(tmp)[bc],
                                          Tax = colnames(tmp)[i],
                                          Var = ii,
                                          Correlation = res$estimate,
                                          p_val = res$p.value,
                                          Method = res$Method))
        }

    }
}

corr_df <- corr_df[corr_df$Tax != "Unknown", ]
cor_genus$Correlation[ is.nan(cor_genus$Correlation) ] <- 0
cor_genus$p_val[ is.nan(cor_genus$p_val) ] <- 1
cor_genus$p_clean <- cor_genus$p_val
cor_genus$p_clean <- format( round(cor_genus$p_clean, 4), scientific = F ) # avoid exponentials
cor_genus$p_clean[cor_genus$p_clean >= 0.05] <- ""
# keep only genera with p < 0.05
corr_genera <- cor_genus$Tax[cor_genus$p_val < 0.05]

fig_4b <- cor_genus %>%
    dplyr::filter(Tax %in% corr_genera) %>%
    dplyr::mutate(Time = gsub(pattern = "Blood_glucose_", replacement = "", Time)) %>%
    dplyr::mutate(Time = gsub(pattern = "_", replacement = "", Time)) %>%
    dplyr::mutate(Time = gsub(pattern = "ipGTTAUC", replacement = "AUC ipGTT", Time)) %>%
    dplyr::mutate(Time = factor(x = Time,
                                levels = c('0min', '45min', 'Weight') ) ) %>%
    ggplot(data = ., aes(x=Var, y=factor(Tax, levels = rev(levels(factor(Tax)))), fill=Correlation)) +
    geom_tile() +
    geom_text(aes(Var, Tax,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    facet_wrap(Time~., scales = "free_x", ncol = 3, labeller = label_value, ) +

    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig_4b

# Figure 5 ----------------------------------------------------------------

# concentrating on week 8 HFD samples
ps_aldex2 <- subset_samples(physeq = ps_kegg, Week == "w8" & Diet == 'HFD')
y <- sample_data(ps_aldex2)$Strain
cnts_aldex2 <- otu_table(ps_aldex2) %>% data.frame() %>% t()

aldex.clr <- aldex.clr(
    reads = cnts_aldex2,
    conds = y,
    mc.samples = 512, denom = "all", verbose = TRUE, useMC = TRUE)

picr.aldex.ttest <- aldex.ttest(aldex.clr, paired.test=FALSE, verbose=TRUE)
picr.aldex.effect <- aldex.effect(aldex.clr, verbose=TRUE)
picr.aldex.all <- data.frame(picr.aldex.ttest, picr.aldex.effect, row.names = rownames(picr.aldex.ttest))

aldex_tt <- picr.aldex.all %>%
    dplyr::select(Effect = effect, q = we.eBH) %>%
    tibble::rownames_to_column(var = 'KEGG')

fig_5 <- EnhancedVolcano::EnhancedVolcano(
    toptable = aldex_tt,
    lab = aldex_tt$KEGG,
    x = 'Effect', y = 'q', xlim = c(-3,3), ylim = c(0,4),
    pointSize = 0.75, ylab = "-log10 q-values",
    labSize = 4, labCol = "grey25",
    pCutoff = 0.01, FCcutoff = 1, drawConnectors = TRUE, arrowheads = FALSE, colConnectors = "grey45",
    maxoverlapsConnectors = Inf, title = "BPL vs ALR (HFD diet)", subtitle = "ALDEx2 Welch's t-test")
fig_5 <- fig_5 +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
fig_5

rm(aldex.clr)

# Supplementary figure 2 --------------------------------------------------

sigLevel <- 0.05
TAXRANK <- "Phylum"

# Corncob week 0: Compare Strains
corncob_dat <- ps %>%
    subset_samples(physeq = ., Week == "w0") %>%
    tax_glom(physeq = ., taxrank = TAXRANK)
# baseline: BPL
da_analysis <-
    differentialTest(formula = ~ Strain + Sex,
                     phi.formula = ~ Strain + Sex, # model to be fitted to the dispersion
                     formula_null = ~ Sex, # Formula for mean under null, without response
                     phi.formula_null = ~ Strain + Sex, # Formula for overdispersion under null, without response
                     test = "Wald", boot = FALSE, B = 0,
                     data = corncob_dat,
                     fdr_cutoff = sigLevel)
da_plot <- plot(da_analysis, level = c("Phylum")) + ggplot2::xlim(c(-10,10))
da_data <- data.frame(
    Taxa = corncob::otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = corncob_dat, level = TAXRANK),
    Taxa_veri = da_plot$data$taxa,
    p_fdr = da_analysis$p_fdr[da_analysis$p_fdr < sigLevel & !is.na(da_analysis$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax
)
suppl_fig2a <- da_data %>%
    dplyr::mutate(Taxa = gsub(pattern = "_", replacement = " ", Taxa)) %>%
    ggplot(data = ., aes(x = Effect, y = factor(Taxa, levels = rev(levels(factor(Taxa)))))) +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          #axis.ticks.x = element_blank(),
          text = element_text(size = 12),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")
    ) +
    xlim(-4, 4) +
    xlab("Effect size (BPL vs ALR)") + ylab("")
suppl_fig2a

# Corncob week 8: Compare CD
corncob_dat <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == 'CD') %>%
    tax_glom(physeq = ., taxrank = TAXRANK)
da_analysis <-
    differentialTest(formula = ~ Strain + Sex,
                     phi.formula = ~ Strain + Sex, # model to be fitted to the dispersion
                     formula_null = ~ Sex, # Formula for mean under null, without response
                     phi.formula_null = ~ Strain + Sex, # Formula for overdispersion under null, without response
                     test = "Wald", boot = FALSE, B = 0,
                     data = corncob_dat,
                     fdr_cutoff = sigLevel)
da_plot <- plot(da_analysis, level = c("Phylum")) + ggplot2::xlim(c(-10,10))
da_data <- data.frame(
    Taxa = corncob::otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = corncob_dat, level = TAXRANK),
    Taxa_veri = da_plot$data$taxa,
    p_fdr = da_analysis$p_fdr[da_analysis$p_fdr < sigLevel & !is.na(da_analysis$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax
)

suppl_fig2b <- da_data %>%
    dplyr::mutate(Taxa = gsub(pattern = "_", replacement = " ", Taxa)) %>%
    ggplot(data = ., aes(x = Effect, y = factor(Taxa, levels = rev(levels(factor(Taxa)))))) +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          #axis.ticks.x = element_blank(),
          text = element_text(size = 12),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")
    ) +
    xlim(-4, 4) +
    xlab("Effect size (BPL vs ALR)") + ylab("")
suppl_fig2b

# Corncob week 8: Compare HFD
corncob_dat <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Diet == 'HFD') %>%
    tax_glom(physeq = ., taxrank = TAXRANK)
da_analysis <-
    differentialTest(formula = ~ Strain + Sex,
                     phi.formula = ~ Strain + Sex, # model to be fitted to the dispersion
                     formula_null = ~ Sex, # Formula for mean under null, without response
                     phi.formula_null = ~ Strain + Sex, # Formula for overdispersion under null, without response
                     test = "Wald", boot = FALSE, B = 0,
                     data = corncob_dat,
                     fdr_cutoff = sigLevel)
da_plot <- plot(da_analysis, level = c("Phylum")) + ggplot2::xlim(c(-10,10))
da_data <- data.frame(
    Taxa = corncob::otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = corncob_dat, level = TAXRANK),
    Taxa_veri = da_plot$data$taxa,
    p_fdr = da_analysis$p_fdr[da_analysis$p_fdr < sigLevel & !is.na(da_analysis$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax
)

suppl_fig2c <- da_data %>%
    dplyr::mutate(Taxa = gsub(pattern = "_", replacement = " ", Taxa)) %>%
    ggplot(data = ., aes(x = Effect, y = factor(Taxa, levels = rev(levels(factor(Taxa)))))) +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          #axis.ticks.x = element_blank(),
          text = element_text(size = 12),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")
    ) +
    xlim(-4, 4) +
    xlab("Effect size (BPL vs ALR)") + ylab("")
suppl_fig2c

# Corncob week 8: Compare ALR
corncob_dat <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == 'ALR') %>%
    tax_glom(physeq = ., taxrank = TAXRANK)
da_analysis <-
    differentialTest(formula = ~ Diet + Sex,
                     phi.formula = ~ Diet + Sex, # model to be fitted to the dispersion
                     formula_null = ~ Sex, # Formula for mean under null, without response
                     phi.formula_null = ~ Diet + Sex, # Formula for overdispersion under null, without response
                     test = "Wald", boot = FALSE, B = 0,
                     data = corncob_dat,
                     fdr_cutoff = sigLevel)
da_plot <- plot(da_analysis, level = c("Phylum")) + ggplot2::xlim(c(-10,10))
da_data <- data.frame(
    Taxa = corncob::otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = corncob_dat, level = TAXRANK),
    Taxa_veri = da_plot$data$taxa,
    p_fdr = da_analysis$p_fdr[da_analysis$p_fdr < sigLevel & !is.na(da_analysis$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax
)

suppl_fig2d <- da_data %>%
    dplyr::mutate(Taxa = gsub(pattern = "_", replacement = " ", Taxa)) %>%
    ggplot(data = ., aes(x = Effect, y = factor(Taxa, levels = rev(levels(factor(Taxa)))))) +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          #axis.ticks.x = element_blank(),
          text = element_text(size = 12),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")
    ) +
    xlim(-4, 4) +
    xlab("Effect size (HFD vs CD)") + ylab("")
suppl_fig2d

# Corncob week 8: Compare BPL
corncob_dat <- ps %>%
    subset_samples(physeq = ., Week == "w8" & Strain == 'BPL') %>%
    tax_glom(physeq = ., taxrank = TAXRANK)
da_analysis <-
    differentialTest(formula = ~ Diet + Sex,
                     phi.formula = ~ Diet + Sex, # model to be fitted to the dispersion
                     formula_null = ~ Sex, # Formula for mean under null, without response
                     phi.formula_null = ~ Diet + Sex, # Formula for overdispersion under null, without response
                     test = "Wald", boot = FALSE, B = 0,
                     data = corncob_dat,
                     fdr_cutoff = sigLevel)
da_plot <- plot(da_analysis, level = c("Phylum")) + ggplot2::xlim(c(-10,10))
da_data <- data.frame(
    Taxa = corncob::otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = corncob_dat, level = TAXRANK),
    Taxa_veri = da_plot$data$taxa,
    p_fdr = da_analysis$p_fdr[da_analysis$p_fdr < sigLevel & !is.na(da_analysis$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax
)

suppl_fig2e <- da_data %>%
    dplyr::mutate(Taxa = gsub(pattern = "_", replacement = " ", Taxa)) %>%
    ggplot(data = ., aes(x = Effect, y = factor(Taxa, levels = rev(levels(factor(Taxa)))))) +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          #axis.ticks.x = element_blank(),
          text = element_text(size = 12),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")
    ) +
    xlim(-4, 4) +
    xlab("Effect size (HFD vs CD)") + ylab("")
suppl_fig2e

# Save workspace ----------------------------------------------------------

save.image(file='data/99_plots_paper.RData')

#
# Some pretty panel plotting
# Fig1
# Fig2
# Fig3
# Fig4
# Fig5

# Suppl Fig 1
# Suppl Fig 2

# Fig 1 a,b ---------------------------------------------------------------

lay <- rbind(c(1,1,1,2,2))

# use arrangeGrob instead of grid.arrange to not draw the figure
fig1 <- gridExtra::arrangeGrob(fig_1a, fig_1b,
                               layout_matrix = lay) %>%
    as_ggplot(.) +
    cowplot::draw_plot_label(label = c("A", "B"), size = 15,
                             x = c(0, 0.5)
    )
fig1
ggsave(filename = "plots_paper/panelplot_Fig1.pdf", width = 12, height = 6, units = 'in', plot = fig1)

# Fig 2 a,b ---------------------------------------------------------------

lay <- rbind(c(1,1,1,2,2,2))

# use arrangeGrob instead of grid.arrange to not draw the figure
fig2 <- gridExtra::arrangeGrob(fig_2a, fig_2b,
                               layout_matrix = lay) %>%
    as_ggplot(.) +
    cowplot::draw_plot_label(label = c("A", "B"), size = 15,
                             x = c(0, 0.5)
    )
fig2
ggsave(filename = "plots_paper/panelplot_Fig2.pdf", width = 12, height = 6, units = 'in', plot = fig2)

# Fig 3 a..j --------------------------------------------------------------

lay <- rbind(c(1,2),
             c(3,4),
             c(5,6),
             c(7,8),
             c(9,10))

# use arrangeGrob instead of grid.arrange to not draw the figure
fig3 <- gridExtra::arrangeGrob(fig_3a, fig_3b,
                               fig_3c, fig_3d,
                               fig_3e, fig_3f,
                               fig_3g, fig_3h,
                               fig_3i, fig_3j,
                               layout_matrix = lay) %>%
    as_ggplot(.) +
    cowplot::draw_plot_label(label = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
                             size = 15,
                             x = c(0.05, 0.55, 0.05, 0.55, 0.05, 0.55, 0.05, 0.55, 0.05, 0.55),
                             y = c(1, 1, 0.8, 0.8, 0.6, 0.6, 0.4, 0.4, 0.2, 0.2)
    )
# fig3
ggsave(filename = "plots_paper/panelplot_Fig3.pdf", width = 24, height = 30, units = 'in', plot = fig3)

# Fig 4 a,b ---------------------------------------------------------------

lay <- rbind(c(1,1,1),
             c(2,2,2))

# use arrangeGrob instead of grid.arrange to not draw the figure
fig4 <- gridExtra::arrangeGrob(fig_4a, fig_4b,
                               layout_matrix = lay) %>%
    as_ggplot(.) +
    cowplot::draw_plot_label(label = c("A", "B"), size = 15,
                             x = c(0, 0),
                             y = c(1, 0.50)
    )
fig4
ggsave(filename = "plots_paper/panelplot_Fig4.pdf", width = 12, height = 12, units = 'in', plot = fig4)

# Fig 5 -------------------------------------------------------------------

fig_5
ggsave(filename = "plots_paper/panelplot_Fig5.pdf", width = 9, height = 9, units = 'in', plot = fig_5)

# Suppl Fig 1 -------------------------------------------------------------

suppl_fig1
ggsave(filename = "plots_paper/panelplot_Supp_Fig1.pdf", width = 6, height = 6, units = 'in', plot = suppl_fig1)

# Suppl Fig 2 -------------------------------------------------------------

lay <- rbind(c(1,2,3),
             c(4,5,NA))

# use arrangeGrob instead of grid.arrange to not draw the figure
suppl_fig2 <- gridExtra::arrangeGrob(suppl_fig2a, suppl_fig2b, suppl_fig2c,
                                    suppl_fig2d, suppl_fig2e,
                               layout_matrix = lay) %>%
    as_ggplot(.) +
    cowplot::draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 15,
                             x = c(0, 0.33, 0.66, 0, 0.33),
                             y = c(1, 1, 1, 0.5, 0.5)
    )
suppl_fig2
ggsave(filename = "plots_paper/panelplot_Supp_Fig2.pdf", width = 15, height = 9, units = 'in', plot = suppl_fig2)
