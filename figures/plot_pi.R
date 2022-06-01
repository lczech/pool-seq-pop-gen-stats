#Basic plotting
library(ggplot2)
library(cowplot)
#library(tidyverse)
theme_set(theme_cowplot())
suppressMessages(library(viridis))
library(stringr)

data <- read.csv("sim_results_pi.csv")
data$group <- ifelse(
    grepl("const", data$fname),
    "const",
    ifelse(
        grepl("growth", data$fname),
        "growth",
        "stat"
    )
)
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
        ag=paste0("factor(",ag,")")
    }

    # Plot stuff. We place the legend not at the side, as that seems to be the easiest way
    # to get all plots to have the same actual width of the plot, without having to put them
    # all into the same grid... https://stackoverflow.com/q/16255579/4184258 screw R.
    # Also, there are some points that are marginally outside of the range... let's keep them
    # by transforming the coordinates instead of limiting the axes. Why oh why.
    ggplot(data, aes_string(x=ax, y=ay, color=ag)) +
        geom_abline(slope=1) +
        geom_point() +
        scale_color_manual(
            values=viridis::viridis(numlvl, alpha = 0.8, begin = 0.1, end = 0.8, direction = -1)
        ) +
        # xlim(0, 0.5) +
        # ylim(-0.1, 0.5) +
        coord_cartesian(xlim=c(0, 0.3), ylim=c(0, 0.3)) +
        xlab(tx) +
        ylab(ty) +
        labs(color=tt) +
        theme(legend.position = "top")
    ggsave(paste0("pi_plot/",ax,"-",ay,"-",ag,".png"), bg="white", width=8, height=8)
}

########################################################################################################
# Group const/stat
########################################################################################################

print( "At group" )

make_plot(
    "true_pairwise_het", "est_pi_within", "group", FALSE,
    "True Pairwise Heterozygosity", "Pi Within"
)

########################################################################################################
# Factor groups
########################################################################################################

for (ag in c( "pool_size", "read_depth", "seq_error", "num_sites" )) {

    # Even the print statement can't take a list of stuff, we need the super obviously named
    # paste0 here to get this to work. All right, R, keep your secrets.
    print(paste0( "At ", ag ))

    make_plot(
        "true_pairwise_het", "est_pi_within", ag, TRUE,
        "True Pairwise Heterozygosity", "Pi Within"
    )

}

warnings()
