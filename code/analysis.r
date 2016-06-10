io = modules::import('ebi-predocs/ebits/io')
dplyr = modules::import_package('dplyr', attach = TRUE)
tidyr = modules::import_package('tidyr')
cu = modules::import('klmr/codons/codon_usage')

dir.create('data', showWarnings = FALSE)

valid_cds = function (cds) {
    n = nchar(cds)
    n > 0 &
        cds != 'Sequence unavailable' &
        n %% 3 == 0 &
        substr(cds, 1, 3) == 'ATG' &
        substr(cds, n - 2, n) %in% cu$stop_codons
}

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
    filter(valid_cds(Sequence)) %>%
    group_by(Gene) %>%
    arrange(-nchar(Sequence)) %>%
    slice(1) %>%
    ungroup()

cds_cu = cu$cu(coding_sequences)

aggregate_vsd = function (data)
    data %>%
    select(Gene, matches('-vsd$')) %>%
    `colnames<-`(sub('-vsd', '', colnames(.))) %>%
    rowwise() %>%
    mutate(Control = mean(Control.1, Control.2, Control.3),
           Treatment = mean(Treatment.1, Treatment.2, Treatment.3)) %>%
    ungroup() %>%
    select(Gene, Control, Treatment)

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

de_diff_cu = inner_join(cds_cu, starvation_de, by = 'Gene') %>%
    mutate(Value = CU * Value) %>%
    group_by(Which, Codon) %>%
    summarize(Value = sum(Value)) %>%
    mutate(Value = Value / sum(Value)) %>%
    ungroup() %>%
    mutate(Interesting = Codon %in% mcm5s2U_codons)

ggplot(de_diff_cu) +
    aes(Codon, Value, fill = Which, alpha = Interesting) +
    scale_alpha_discrete(range = c(0.5, 1)) +
    geom_bar(stat = 'identity', position = 'dodge')

de_diff_d = de_diff_cu %>%
    tidyr$spread(Which, Value) %>%
    mutate(Difference = Treatment - Control, magnitude = abs(Difference))

# H0: codon usage varies randomly between highly expressed genes in FF and FS,
# and this is also true for the codons of interest.

# <http://stats.stackexchange.com/a/62653/3512>
pred_interval_d = de_diff_d %>% filter(! Interesting) %>% .$Difference
# Actually we know the population mean under the null hypothesis (= 0) but we
# may as well estimate it from the data, and indeed it’s almost 0.
pred_µ = mean(pred_interval_d)
pred_σ = sd(pred_interval_d)
λ_95 = 3
# P(|X - μ| ≥ λσ) ≤ 4/(9λ²) for λ > √(8/3)
# → 1 - P > 0.95

de_diff_d = de_diff_d %>%
    mutate(λ = abs(Difference - pred_µ) / pred_σ,
           p = ifelse(λ > sqrt(8 / 3), 4 / (9 * λ ^ 2), 1),
           Significance = symnum(p, corr = FALSE,
                                 cutpoints = c(0, 0.01, 0.05, 1),
                                 symbols = c('**', '*', ' ')) %>% as.character)

ggplot(de_diff_d) +
    aes(Codon, Difference, fill = factor(Interesting, labels = c('Yes', 'No'))) +
    geom_bar(stat = 'identity', position = 'dodge') +
    geom_hline(yintercept = pred_µ + λ_95 * c(1, -1) * pred_σ) +
    annotate('text', label = '~95% prediction interval', x = 30, y = 0.007) +
    geom_text(aes(label = Significance), vjust = 1.5) +
    labs(fill = 'Codon of interest',
         title = 'Difference in codon usage between Control and Treatment')
