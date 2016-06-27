cu = modules::import('klmr/codons/codon_usage')

valid_cds = function (cds) {
    n = nchar(cds)
    n > 0 &
        cds != 'Sequence unavailable' &
        n %% 3 == 0 &
        substr(cds, 1, 3) == 'ATG' &
        substr(cds, n - 2, n) %in% cu$stop_codons
}
