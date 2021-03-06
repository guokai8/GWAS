qq_plot
```
gg_qqplot = function(xs, ci=0.95) {
    N = length(xs)
    df = data.frame(observed=-log10(sort(xs)),
                    expected=-log10(1:N / N),
                    cupper=-log10(qbeta(ci,     1:N, N - 1:N + 1)),
                    clower=-log10(qbeta(1 - ci, 1:N, N - 1:N + 1)))
    log10Pe = expression(paste("Expected -log"[10], plain(P)))
    log10Po = expression(paste("Observed -log"[10], plain(P)))
    ggplot(df) +
        geom_point(aes(expected, observed), shape=1, size=1) +
        geom_abline(intercept=0, slope=1, alpha=0.5) +
        geom_line(aes(expected, cupper), linetype=2) +
        geom_line(aes(expected, clower), linetype=2) +
        xlab(log10Pe) +
        ylab(log10Po)
}
```
gg_plot(pvalue)
