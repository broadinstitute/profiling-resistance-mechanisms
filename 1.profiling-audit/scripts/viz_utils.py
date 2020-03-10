import plotnine as gg


def save_figure(
    main_figure,
    file_base,
    extensions=[".png", ".pdf", ".svg"],
    dpi=300,
    height=6,
    width=8,
):
    for extension in extensions:
        output_file = "{}{}".format(file_base, extension)
        main_figure.save(
            filename=output_file, dpi=dpi, height=height, width=width, verbose=False
        )


def plot_replicate_correlation(
    df,
    batch,
    plate,
    facet_string=None,
    split_samples=False,
    output_file_base=None,
    output_file_extensions=[".png", ".pdf", ".svg"],
    dpi=300,
    height=3,
    width=3.5,
):
    correlation_gg = (
        gg.ggplot(
            df,
            gg.aes(x="replicate_info", y="pairwise_correlation", fill="replicate_info"),
        )
        + gg.geom_boxplot(
            alpha=0.3, outlier_alpha=0, width=0.8, notchwidth=0.25, fatten=1.5
        )
        + gg.geom_jitter(shape=".", size=0.001, alpha=0.3, width=0.3, height=0)
        + gg.scale_fill_manual(
            name="Replicate",
            labels={"True": "True", "False": "False"},
            values=["#B99638", "#2DB898"],
        )
        + gg.xlab("Replicates")
        + gg.ylab("Pearson Correlation")
        + gg.ggtitle("{}: {}".format(batch, plate))
        + gg.theme_bw()
        + gg.theme(
            subplots_adjust={"wspace": 0.2},
            title=gg.element_text(size=9),
            axis_text=gg.element_text(size=5),
            axis_title=gg.element_text(size=8),
            legend_text=gg.element_text(size=6),
            legend_title=gg.element_text(size=7),
            strip_text=gg.element_text(size=4, color="black"),
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
        )
    )

    if split_samples:
        assert facet_string, "To split samples, specify a facet_string"
        correlation_gg += gg.facet_wrap(facet_string)

    if output_file_base:
        save_figure(
            correlation_gg, output_file_base, output_file_extensions, dpi, height, width
        )

    return correlation_gg


def plot_replicate_density(
    df,
    batch,
    plate,
    output_file_base=None,
    output_file_extensions=[".png", ".pdf", ".svg"],
    dpi=300,
    height=1.5,
    width=2,
):
    density_gg = (
        gg.ggplot(df, gg.aes(x="pairwise_correlation", fill="replicate_info"))
        + gg.geom_density(alpha=0.3)
        + gg.scale_fill_manual(
            name="Replicate",
            labels={"True": "True", "False": "False"},
            values=["#B99638", "#2DB898"],
        )
        + gg.xlab("Pearson Correlation")
        + gg.ylab("Density")
        + gg.ggtitle("{}: {}".format(batch, plate))
        + gg.theme_bw()
        + gg.theme(
            title=gg.element_text(size=9),
            axis_text=gg.element_text(size=5),
            axis_title=gg.element_text(size=8),
            legend_text=gg.element_text(size=6),
            legend_title=gg.element_text(size=7),
            strip_text=gg.element_text(size=4, color="black"),
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
        )
    )

    if output_file_base:
        save_figure(
            density_gg, output_file_base, output_file_extensions, dpi, height, width
        )

    return density_gg
