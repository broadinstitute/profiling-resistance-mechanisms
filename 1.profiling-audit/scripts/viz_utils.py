import numpy as np
import plotnine as gg


def save_figure(
    main_figure,
    file_base,
    extensions=[".png", ".pdf", ".svg"],
    dpi=500,
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
    dpi=500,
    height=4,
    width=5,
    return_plot=False
):
    correlation_gg = (
        gg.ggplot(
            df,
            gg.aes(x="group_replicate", y="similarity_metric", fill="group_replicate"),
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
            title=gg.element_text(size=5),
            axis_text=gg.element_text(size=4),
            axis_title=gg.element_text(size=5),
            legend_text=gg.element_text(size=4),
            legend_title=gg.element_text(size=5),
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
    if return_plot:
        return correlation_gg


def plot_replicate_density(
    df,
    batch,
    plate,
    cutoff,
    percent_strong,
    output_file_base=None,
    output_file_extensions=[".png", ".pdf", ".svg"],
    dpi=300,
    height=1.5,
    width=2,
    return_plot=False
):
    density_gg = (
        gg.ggplot(df, gg.aes(x="similarity_metric", fill="group_replicate"))
        + gg.geom_density(alpha=0.3)
        + gg.scale_fill_manual(
            name="Replicate",
            labels={"True": "True", "False": "False"},
            values=["#B99638", "#2DB898"],
        )
        + gg.xlab("Pearson Correlation")
        + gg.ylab("Density")
        + gg.geom_vline(xintercept=cutoff, color="red", linetype="dashed")
        + gg.ggtitle(f"Batch: {batch}; Plate: {plate}\nPercent Strong: {np.round(percent_strong * 100, 2)}%")
        + gg.theme_bw()
        + gg.theme(
            title=gg.element_text(size=5),
            axis_text=gg.element_text(size=4),
            axis_title=gg.element_text(size=5),
            legend_text=gg.element_text(size=4),
            legend_title=gg.element_text(size=5),
            strip_text=gg.element_text(size=4, color="black"),
            strip_background=gg.element_rect(colour="black", fill="#fdfff4"),
        )
    )

    if output_file_base:
        save_figure(
            density_gg, output_file_base, output_file_extensions, dpi, height, width
        )
    
    if return_plot:
        return density_gg
