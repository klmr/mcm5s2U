#+ echo=FALSE
knitr::opts_chunk$set(dev = c('png', 'pdf'))

#+ echo=TRUE
io = modules::import('ebi-predocs/ebits/io')
dplyr = modules::import_package('dplyr', attach = TRUE)
tidyr = modules::import_package('tidyr')
cds_ = modules::import('./cds')

dir.create('data', showWarnings = FALSE)

load_counts = function (file) {
    data = io$read_table(file, header = TRUE, sep = ',', row.names = 1)
    colnames(data) = sub('_STAR_WS220\\.mis2\\.map10\\.bam', '',
                                sub('.*RNA_', '', colnames(data)))
    colnames(data) = sub('FF', 'Control', colnames(data))
    colnames(data) = sub('FS|37C', 'Treatment', colnames(data))
    data
}

starvation_data = load_counts('raw-data/deseq_WS220_genes_RNA_FF-vs-FS.txt')
heatshock_data = load_counts('raw-data/deseq_WS220_genes_RNA_FF-vs-37C.txt')

load_cds = function () {
    path = 'data/cds.rds'
    cds = try(readRDS(path), silent = TRUE)

    if (inherits(cds, 'try-error')) {
        bm = modules::import_package('biomaRt')
        ensembl = bm$useMart('ensembl', 'celegans_gene_ensembl')
        cds = bm$getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'coding'),
                       mart = ensembl, bmHeader = TRUE) %>% tbl_df()
        saveRDS(cds, path)
    }

    cds
}

cds = load_cds()

coding_sequences = cds %>%
    select(`Ensembl Gene ID`,
           Gene = `Associated Gene Name`,
           Sequence = `Coding sequence`) %>%
    filter(cds_$valid_cds(Sequence)) %>%
    group_by(Gene) %>%
    arrange(-nchar(Sequence)) %>%
    slice(1) %>%
    ungroup()

cds_cu = cds_$cu$cu(coding_sequences)

aggregate_vsd = function (data)
    data %>%
    select(Gene, matches('-vsd$')) %>%
    `colnames<-`(sub('-vsd', '', colnames(.))) %>%
    rowwise() %>%
    mutate(Control = mean(Control.1, Control.2, Control.3),
           Treatment = mean(Treatment.1, Treatment.2, Treatment.3)) %>%
    ungroup() %>%
    select(Gene, Control, Treatment)
# FIXME: Estimate variance, and use in calculation of statistic below.

starvation_vsd = aggregate_vsd(starvation_data)
heatshock_vsd = aggregate_vsd(heatshock_data)

mcm5s2U_codons = io$read_table('raw-data/mcm5s2U-codons.tsv')$V2

modules::import_package('ggplot2', attach = TRUE)
theme_set(theme_bw())

filter_de_genes = function (data)
    data %>%
    filter(padj < 0.01) %>%
    mutate(Which = ifelse(log2FoldChange < 0, 'Control', 'Treatment')) %>%
    inner_join(starvation_vsd, by = 'Gene') %>%
    mutate(Value = ifelse(Which == 'Control', Control, Treatment)) %>%
    select(Gene, Which, Value)

starvation_de = filter_de_genes(starvation_data)
heatshock_de = filter_de_genes(heatshock_data)

codon_usage_summary = function (de_data)
    inner_join(cds_cu, de_data, by = 'Gene') %>%
    mutate(Value = CU * Value) %>%
    group_by(Which, Codon) %>%
    summarize(Value = sum(Value)) %>%
    mutate(Value = Value / sum(Value)) %>%
    ungroup() %>%
    mutate(Interesting = Codon %in% mcm5s2U_codons)

starvation_de_cu = codon_usage_summary(starvation_de)
heatshock_de_cu = codon_usage_summary(heatshock_de)

#+ starvation-de-codon-usage, fig.width=20
ggplot(starvation_de_cu) +
    aes(Codon, Value, fill = Which, alpha = Interesting) +
    scale_alpha_discrete(range = c(0.5, 1)) +
    geom_bar(stat = 'identity', position = 'dodge')

#+ heatshock-de-codon-usage, fig.width=20
ggplot(heatshock_de_cu) +
    aes(Codon, Value, fill = Which, alpha = Interesting) +
    scale_alpha_discrete(range = c(0.5, 1)) +
    geom_bar(stat = 'identity', position = 'dodge')

starvation_diff = starvation_de_cu %>%
    tidyr$spread(Which, Value) %>%
    mutate(Difference = Treatment - Control)

heatshock_diff = heatshock_de_cu %>%
    tidyr$spread(Which, Value) %>%
    mutate(Difference = Treatment - Control)

# H0: codon usage varies randomly between highly expressed genes in FF and FS,
# and this is also true for the codons of interest.

lambda = 3

test_enrichment = function (de_diff_d) {
    # <http://stats.stackexchange.com/a/62653/3512>
    pred_interval_d = de_diff_d %>% filter(! Interesting) %>% .$Difference
    # P(|X - μ| ≥ λσ) ≤ 4/(9λ²) for λ > √(8/3)
    # → 1 - P > 0.95

    # Actually we know the population mean under the null hypothesis (= 0) but we
    # may as well estimate it from the data, and indeed it’s almost 0.
    de_diff_d %>%
    mutate(mean = mean(pred_interval_d),
           sd = sd(pred_interval_d),
           λ = abs(Difference - mean) / sd,
           # FIXME: This isn’t a p-value! It’s the probability of H0.
           p = ifelse(λ > sqrt(8 / 3), 4 / (9 * λ ^ 2), 1),
           Significance = symnum(p, corr = FALSE,
                                 cutpoints = c(0, 0.01, 0.05, 1),
                                 symbols = c('**', '*', ' ')) %>% as.character)
}

starvation_enrichment = test_enrichment(starvation_diff)
heatshock_enrichment = test_enrichment(heatshock_diff)

interval_lines = function (data)
    geom_hline(yintercept = first(data$mean) + lambda * c(1, -1) * first(data$sd))

#+ starvation-enrichment, fig.width=20
ggplot(starvation_enrichment) +
    aes(Codon, Difference, fill = factor(Interesting, labels = c('mcm5s2U', 'rest'))) +
    geom_bar(stat = 'identity', position = 'dodge') +
    interval_lines(starvation_enrichment) +
    annotate('text', label = '~95% prediction interval', x = 25, y = 0.008) +
    geom_text(aes(label = Significance), vjust = -0.5) +
    labs(fill = 'Codon type', y = 'Difference between FF and FS')

#+ heatshock-enrichment, fig.width=20
ggplot(heatshock_enrichment) +
    aes(Codon, Difference, fill = factor(Interesting, labels = c('mcm5s2U', 'rest'))) +
    geom_bar(stat = 'identity', position = 'dodge') +
    interval_lines(heatshock_enrichment) +
    annotate('text', label = '~95% prediction interval', x = 30, y = 0.007) +
    geom_text(aes(label = Significance), vjust = -0.5) +
    labs(fill = 'Codon type', y = 'Difference, between FF and 37°C')

## Summarise the plots for easier digestion.

plot_summary = function (data) {
    data = data %>%
        mutate(Codon = ifelse(Interesting, Codon, 'rest')) %>%
        # Alper wants this order of codons.
        mutate(Codon = factor(Codon,
                              levels = c('AAA', 'GAA', 'CAA', 'rest'),
                              ordered = TRUE))

    distr = filter(data, ! Interesting)$Difference
    mean = mean(distr)
    sd = sd(distr)

    ggplot(data) +
        aes(Codon, Difference) +
        geom_blank() +
        geom_bar(data = filter(data, Codon != 'rest'),
                 stat = 'identity', position = 'dodge', width = 0.5) +
        geom_boxplot(data = filter(data, Codon == 'rest'),
                     size = 1.2, width = 0.75) +
        geom_hline(aes(yintercept = Limit, color = Limits),
                   data = data_frame(Limit = mean + lambda * c(1, -1) * sd,
                                     Limits = '95% prediction interval'),
                   linetype = 'dashed', show.legend = TRUE) +
        coord_cartesian(ylim = c(-0.01, 0.01))
}

#+ starvation-enrichment-summary, fig.width = 6, fig.height = 4
plot_summary(starvation_enrichment)

#+ heatshock-enrichment-summary, fig.width = 6, fig.height = 4
plot_summary(heatshock_enrichment)
