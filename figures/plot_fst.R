#Basic plotting
library(ggplot2)
library(cowplot)
#library(tidyverse)
theme_set(theme_cowplot())
suppressMessages(library(viridis))
library(stringr)

data <- read.csv("sim_results_fst.csv")
data$group <- ifelse(grepl("const", data$fname), "msprime", "simple")
#head(data)

# Make a plot. R is trying to be funny again - the way we get color grouping with factors to work
# is by string evaluation of the string "factor(column)"... WTF
make_plot <- function( ax, ay, ag, fact, tx, ty ) {
    tt = str_to_title(str_replace(ag, "_", " "))

    # Holy craps, we need to find the number of levels from a column specified as a string...
    # Who invented this language?
    # numlvl = eval(parse(text=paste0( "levels(factor(data$",ag,"))")))
    # numlvl = eval(paste0( "levels(factor(data$",ag,"))"))

    # Nope neither works, try something else. Double brackets apparently make a column
    # addressable by string. So obvious.
    numlvl = nlevels(factor(data[[ag]]))
    if( fact ){
        # print(levels(factor(data[[ag]])))
        # sorted_labels <- paste(sort(as.double(levels(factor(data[[ag]])))))
        # print(sorted_labels)
        # data[[ag]] <- factor(data[[ag]], levels = sorted_labels)
        colorag=paste0("factor(",ag,")")
    }

    # Plot stuff. We place the legend not at the side, as that seems to be the easiest way
    # to get all plots to have the same actual width of the plot, without having to put them
    # all into the same grid... https://stackoverflow.com/q/16255579/4184258 screw R.
    # Also, there are some points that are marginally outside of the range... let's keep them
    # by transforming the coordinates instead of limiting the axes. Why oh why.
    ggplot(data, aes_string(x=ax, y=ay, color=colorag)) +
        geom_abline(slope=1) +
        geom_point() +
        scale_color_manual(
            values=viridis::viridis(numlvl, alpha = 0.8, begin = 0.1, end = 0.8, direction = -1)
        ) +
        # xlim(0, 0.5) +
        # ylim(-0.1, 0.5) +
        # coord_cartesian(xlim=c(0, 0.5), ylim=c(-0.07, 0.5)) +
        coord_fixed(xlim=c(0, 0.5), ylim=c(-0.07, 0.5), ratio = 1) +
        xlab(tx) +
        ylab(ty) +
        labs(color=tt) +
        theme(legend.position = "top") +
        facet_wrap(vars(group), ncol = 2) +
        theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position = c(0.87, 0.15))

    # ggplot does not have any proper way of setting one image dimension based on the other and
    # the aspect ratio... WTF... see https://stackoverflow.com/q/16422847/4184258
    # So what we do here instead is just eyeballing it to achieve as little white space as possible
    # around the image. This sucks.
    # ggsave(paste0("fst_plot_png/",ax,"-",ay,"-",ag,".png"), bg="white", width=12, height=7.7)
    ggsave(paste0("fst_plot_pdf/",ax,"-",ay,"-",ag,".pdf"), bg="white", width=9, height=5.5)
}

########################################################################################################
#     Group const/stat
########################################################################################################

# Not interesting to have these plots, as we are facet splitting by group now anyway.

# print( "At group" )
#
# # Spence Nei, Spence Hudson
# make_plot(
#     "true_hudson_fst", "est_spence_hudson", "group", FALSE,
#     "True FST (Hudson)", "Spence Estimator FST (Hudson)"
# )
# make_plot(
#     "true_nei_fst", "est_spence_nei", "group", FALSE,
#     "True FST (Nei)", "Spence Estimator FST (Nei)"
# )
#
# # Kofler
# # make_plot(
# #     "true_hudson_fst", "est_kofler", "group", FALSE,
# #     "True FST (Hudson)", "Kofler Estimator FST"
# # )
# make_plot(
#     "true_nei_fst", "est_kofler", "group", FALSE,
#     "True FST (Nei)", "Kofler Estimator FST"
# )
#
# # Karlsson
# make_plot(
#     "true_hudson_fst", "est_karlsson", "group", FALSE,
#     "True FST (Hudson)", "Karlsson Estimator FST"
# )
# # make_plot(
# #     "true_nei_fst", "est_karlsson", "group", FALSE,
# #     "True FST (Nei)", "Karlsson Estimator FST"
# # )

########################################################################################################
#     Factor groups
########################################################################################################

# Plotting the interesting vars. We are not plotting "num_sites" here,
# as this is just a constant num for the const/msprime, and another number for the stat/simple model.
for (ag in c( "pool_size", "read_depth", "seq_error" )) {

    # Even the print statement can't take a list of stuff, we need the super obviously named
    # paste0 here to get this to work. All right, R, keep your secrets.
    print(paste0( "At ", ag ))

    # Spence Nei, Spence Hudson
    make_plot(
        "true_hudson_fst", "est_spence_hudson", ag, TRUE,
        "True FST (Hudson)", "Unbiased Pool-Seq Estimator FST (Hudson)"
    )
    make_plot(
        "true_nei_fst", "est_spence_nei", ag, TRUE,
        "True FST (Nei)", "Unbiased Pool-Seq Estimator FST (Nei)"
    )

    # Kofler
    # make_plot(
    #     "true_hudson_fst", "est_kofler", ag, TRUE,
    #     "True FST (Hudson)", "Kofler Estimator FST"
    # )
    make_plot(
        "true_nei_fst", "est_kofler", ag, TRUE,
        "True FST (Nei)", "Kofler Estimator FST"
    )

    # Karlsson
    make_plot(
        "true_hudson_fst", "est_karlsson", ag, TRUE,
        "True FST (Hudson)", "Karlsson Estimator FST"
    )
    # make_plot(
    #     "true_nei_fst", "est_karlsson", ag, TRUE,
    #     "True FST (Nei)", "Karlsson Estimator FST"
    # )

}

# warnings()
