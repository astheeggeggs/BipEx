library(ggplot2)
library(dplyr)
library(ggsci)
library(latex2exp)
library(ggrastr)
library(ggrepel)

create_pretty_forest <- function(df, title, save_figure=FALSE, file='file_out',
    width=160, height=90, scaling=1, y_label='', pass_order=NULL, print_p=TRUE,
    horizontal=TRUE, hline_at=0, ggplot_theme=theme_classic)
{
    if (!is.null(pass_order)) {
        df$label <- factor(df$label, levels = df$label[pass_order])
    }

    p <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
        geom_pointrange() + 
        geom_hline(yintercept=hline_at, lty=2)

    if (horizontal==TRUE) p <- p + coord_flip()

    p <- p + labs(title=title, y=y_label, x='') +
        ggplot_theme()

    if (horizontal == TRUE) {
        p <- p + geom_text(aes(y=mean, label=p_vals), vjust=-1)
    } else {
        p <- p + geom_text(aes(y=mean, label=p_vals), hjust=-0.2)
    }

    if (print_p)  print(p)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
    }
    
    return(p)
}

create_pretty_hist <- function(df, aest, x_label, threshold=NULL, threshold_max=NULL,
    file='file_out', title='', binwidth=0.002, width=160, height=90, scaling=1,
    save_figure=FALSE, xlim=NULL, key_label='', print_p=TRUE, title.hjust=0.5, ggplot_theme=theme_classic)
{
    p <- ggplot(df, aest)
    geom_histogram(binwidth=binwidth, fill='#aec7e8', color='#1f77b4')
    ylim <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range
    p <- ggplot(df, aest)

    if (is.null(aest$fill)) {
        p <- p + geom_histogram(binwidth=binwidth, fill='#aec7e8', color='#1f77b4')
    } else {
        p <- p + geom_histogram(binwidth=binwidth, color='grey50')
    }

    if (!is.null(xlim)) {
        p <- p + coord_cartesian(xlim=xlim)
    } 

    p <- p + labs(x=x_label, y='Count', title=title, fill=key_label) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
    scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
    ggplot_theme() +
    theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
    axis.title.y = element_text(margin=ggplot2::margin(r=10)),
    plot.title = element_text(hjust=title.hjust))

    if (!is.null(threshold)) p <- p + geom_vline(xintercept=threshold, linetype='dashed')
    if (!is.null(threshold_max)) p <- p + geom_vline(xintercept=threshold_max, linetype='dashed')

    
    if (print_p) print(p)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
    }
    
    return(p)
}

create_pretty_density <- function(df, aest, x_label, threshold=NULL, threshold_max=NULL,
    file='file_out', title='', binwidth=0.002, width=160, height=90, scaling=1,
    save_figure=FALSE, xlim=NULL, key_label='', print_p=FALSE, title.hjust=0.5,
    ggplot_theme=theme_classic)
{
    p <- ggplot(df, aest)
    ylim <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range
    p <- ggplot(df, aest)

    if (is.null(aest$fill)) {
        p <- p + geom_density(fill='#aec7e8', color='#1f77b4')
    } else {
        p <- p + geom_density(color='grey50')
    }

    if (!is.null(xlim)) {
        p <- p + coord_cartesian(xlim=xlim)
    } 

    p <- p + labs(x=x_label, y='Count', title=title, fill=key_label) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
    scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
    ggplot_theme() +
    theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
    axis.title.y = element_text(margin=ggplot2::margin(r=10)),
    plot.title = element_text(hjust=title.hjust))

    if (!is.null(threshold)) {
        p <- p + geom_vline(xintercept=threshold, linetype='dashed') #+
        # annotate('text', label=paste0(x_label, ' = ', threshold), x=threshold, y=0.8*ylim[2], angle=90, vjust=-1)
    }
    if (!is.null(threshold_max)) {
        p <- p + geom_vline(xintercept=threshold_max, linetype='dashed') #+
        # annotate('text', label=paste0(x_label, ' = ', threshold_max), x=threshold_max, y=0.8*ylim[2], angle=90, vjust=2)
    }

    if (print_p) print(p)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
    }
    
    return(p)
}

create_pretty_boxplots <- function(df, aes, aes_col, threshold=NULL,
    threshold_max=NULL, file='file_out', title='', x_label='', y_label='',
    key_label='', xlim=NULL, legend=FALSE, save_figure=FALSE, 
    width=160, height=90, scaling=1, facet=FALSE, facet_grid=NULL, jitter_size=0.5,
    outlier.shape=NA, n_ticks=10, print_p=FALSE, alpha=0.6, title.hjust=0.5, ggplot_theme=theme_classic)
{
    p = ggplot(df, aes) +
        geom_boxplot(outlier.shape=outlier.shape, coef=0, color='grey50', fill='grey95', show.legend=FALSE) + 
        geom_jitter(width=0.2, height=0, size=jitter_size, aes_col, show.legend=legend, alpha=alpha, stroke=0.05) + 
        coord_flip(ylim=xlim) +
        labs(title=title, x=y_label, y=x_label, color=key_label) + 
        scale_color_d3('category20') +
        scale_y_continuous(breaks = scales::pretty_breaks(n=n_ticks)) +
        guides(color = guide_legend(override.aes = list(size=2))) +
        ggplot_theme() +
        theme(axis.title.x = element_text(margin = ggplot2::margin(t=10)),
              plot.title = element_text(hjust=title.hjust))

    if (!is.null(threshold)) {
        p <- p + geom_hline(yintercept=threshold, linetype='dashed')
    }
    if (!is.null(threshold_max)) {
        p <- p + geom_hline(yintercept=threshold_max, linetype='dashed')
    }
    if (facet){
        p <- p + facet_grid
    }

    if (print_p) print(p)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
    }
    return(p)
}

create_pretty_cumulative <- function(df, aes, x_label, threshold, threshold_max=NULL,
    file='file_out', title='', width=160, height=90, scaling=1, save_figure=FALSE,
    xlim=c(0,1), key_label='', print_p=FALSE, title.hjust=0.5, ggplot_theme=theme_classic)
{
    # These next ones are cdfs.
    p = ggplot(df, aes) + 
        stat_ecdf(geom='line', pad=FALSE) +
        geom_vline(xintercept=threshold, linetype='dashed') +
        coord_cartesian(xlim=xlim) +
        labs(x=x_label, y='Cumulative Proportion of Samples', title=title, color=key_label) +
        scale_color_d3('category10') +
        scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
        ggplot_theme() +
        theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
              axis.title.y = element_text(margin=ggplot2::margin(r=10)),
              plot.title = element_text(hjust=title.hjust))

    # p <- p + annotate('text', label=paste0(x_label, ' = ', threshold), x=threshold, y=mean(ylim), angle=90, vjust=-1)
    # ylim <- ggplot_build(p)$panel$ranges[[1]]$y.range
    # print(ylim)

    if (!is.null(threshold_max)) {
        p <- p + geom_vline(xintercept=threshold_max, linetype='dashed')# +
        # annotate('text', label=paste0(x_label, ' = ', threshold_max), x=threshold_max, y=mean(ylim), angle=90, vjust=2)
    }
    if (print_p) print(p)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
    }
    return(p)
}

create_pretty_scatter <- function(dt, aes, file='file_out', save_figure=FALSE, 
    key_label='', title='', limits=NULL, width=160, height=90, presentation=FALSE,
    add_final_layer=FALSE, final_layer=NULL, n_x_ticks=10, n_y_ticks=10, x_label=NULL,
    y_label=NULL, print_p=FALSE, gradient=FALSE, midpoint=0, gradient_title="", title.hjust=0.5,
    legend_max_min=NULL, n_legend_step=4, xlim=NULL, ylim=NULL, ggplot_theme=theme_classic, alpha=0.6, size=1)
{
    p <- ggplot(dt, aes)
    if ("size" %in% names(aes)) {
        p <- p + geom_point_rast(alpha=alpha)
    } else {
        p <- p + geom_point_rast(alpha=alpha, size=size)
    }
    p <- p +
        scale_color_d3('category20', limits=limits) +
        labs(title=title, color=paste0(key_label)) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=n_x_ticks)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=n_y_ticks)) +
        ggplot_theme() +
        theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
              axis.title.y = element_text(margin=ggplot2::margin(r=10)),
              plot.title = element_text(hjust=title.hjust))

    if (add_final_layer) {
        cat("Adding final layer...\n")
        p <- p + guides(fill=FALSE) + geom_point_rast(mapping=aes, data=final_layer) + scale_color_d3('category20')
    }

    if (gradient) {
        print(gradient_title)
        if(!is.null(legend_max_min)) {
            p <- p + scale_color_gradient2(
                low="blue", high="red", mid='grey50', midpoint=midpoint,
                name=gradient_title,
                breaks=seq(legend_max_min[1], legend_max_min[2], (legend_max_min[2]-legend_max_min[1])/n_legend_step))
        } else {
            p <- p + scale_color_gradient2(
                low="blue", high="red", mid='grey50', midpoint=midpoint,
                name=gradient_title)
        }
    }
    
    cat("Adding axis labels...\n")
    if (!is.null(x_label)) p <- p + labs(x=x_label)
    if (!is.null(y_label)) p <- p + labs(y=y_label)
    if (!is.null(xlim)) p <- p + xlim(xlim[1], xlim[2])
    if (!is.null(ylim)) p <- p + ylim(ylim[1], ylim[2])

    cat("Saving figure...\n")
    if (print_p) print(p)

    if (save_figure) {
        
        if (presentation == TRUE) {
            width <- 160
            height <- 90
            ggsave(paste0(file, '.jpg'), p, width=width, height=height, units='mm')
        } else {
            ggsave(paste0(file, '.jpg'), p, width=width, height=height, units='mm')
            ggsave(paste0(file, '.pdf'), p, width=width, height=height, units='mm')
        }
    
    }

    return(p)
}

create_pretty_qq_plot <- function(dt, aes, file='file_out', save_figure=FALSE,
    plot_title='', limits=NULL, width=110, height=100, n_x_ticks=10, n_y_ticks=10,
    x_label=TeX("$-\\log_{10}(\\mathit{p}_{permutation})$"), 
    y_label=TeX("$-\\log_{10}(\\mathit{p}_{observed})$"),
    n_to_include=NULL, cex_labels=1, print_p=TRUE, gradient=FALSE,
    gradient_title="", pval_col="pval", title.hjust=0.5, legend_max_min=NULL, n_legend_step=4,
    xlim=NULL, ylim=NULL, y_x_col="grey40", ggplot_theme=theme_classic, alpha=0.6)
{
    cat("Creating scatter-plot...\n")
    dt <- data.table(dt)
    setkeyv(dt, cols=pval_col)
    p <- create_pretty_scatter(dt, aes, file=file, save_figure=FALSE,
        title=plot_title, limits=limits,
        width=width, height=height, n_x_ticks=n_x_ticks, n_y_ticks=n_y_ticks,
        x_label=x_label, y_label=y_label, gradient=gradient, gradient_title=gradient_title,
        title.hjust=title.hjust, legend_max_min=legend_max_min, n_legend_step=n_legend_step,
        xlim=xlim, ylim=ylim, ggplot_theme=ggplot_theme, alpha=alpha)

    cat("Adding y=x line...\n")
    p <- p + geom_abline(intercept=0, slope=1, color=y_x_col) #+ coord_fixed()

    if (!is.null(n_to_include)) {
        cat("Adding labels...\n")
        p <- p + geom_label_repel(data=dt[(nrow(dt)-n_to_include+1):nrow(dt), ],
            aes(label=labels), box.padding = 0.5, label.padding=0.1, point.padding = 0.2,
            color = 'grey30', segment.color = 'grey50',
            size=cex_labels, segment.size=0.1, show.legend = FALSE)
    }

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width, height=height, units='mm')
        ggsave(paste0(file, '.pdf'), p, width=width, height=height, units='mm')
    }
    
    cat("Created scatter-plot...\n")
    if (print_p) print(p)

    return(p)
}

create_pretty_pointrange <- function(dt, y, ymin, ymax, colors="phenotype",
    break_on="phenotype", x_labels="phenotype", threshold=NULL, file='file_out',
    title='', width=165, height=110, save_figure=FALSE, xlim=NULL, key_label='',
    print_p=FALSE, title.hjust=0.5, split_on=NULL, spacing=7, x_label="", y_label="",
    title_size=NULL, manual_colours=NULL, colour_levels=NULL, break_on_levels=NULL,
    ggplot_theme=theme_classic, ylim=NULL)
{   
    if(!is.null(colour_levels)) {
        dt[[colors]] <- factor(dt[[colors]], levels=colour_levels)
    } else {
        dt[[colors]] <- factor(dt[[colors]])
    }

    if(!is.null(break_on_levels)) {
        dt[[break_on]] <- factor(dt[[break_on]], levels=break_on_levels)
    } else {
        dt[[break_on]] <- factor(dt[[break_on]])
    }

    dt <- dt %>% arrange(dt[[break_on]], dt[[colors]])

    dt1 <- cbind(dt, aux=rep(1,length(dt[,1]))) 
    dt1 <- within(dt1, {aux = unlist(by(aux, dt1[[break_on]], cumsum))})
    dt1$aux <- dt1$aux + as.numeric(factor(dt1[[break_on]])) * (length(unique(factor(dt1[[break_on]]))) + spacing)

    print(head(dt1))

    define_breaks <- function(dt, split_on) {
        breaks <- c()
        for (split in unique(dt[[split_on]])) {
            breaks <- c(breaks, mean(unlist(dt %>% filter(dt[[split_on]] == split) %>% select(aux))))
        }
        return(breaks)
    }

    p <- ggplot(dt1, aes(x=dt1$aux, y = dt1[[y]], ymin = dt1[[ymin]], ymax = dt1[[ymax]])) +
    geom_pointrange(aes(color = dt1[[colors]]), size = 0.5) + 
    scale_x_continuous("", breaks=define_breaks(dt1, break_on), labels=unique(dt1[[break_on]])) +
    labs(x=x_label, y=y_label, title=title, color=key_label)

    if (!is.null(manual_colours)) {
        p <- p + scale_color_manual(labels=colour_levels, values=manual_colours)
    } else {
        p <- p + scale_color_d3('category10')
    }
    p <- p +
    ggplot_theme() +
    theme(axis.title.x = element_text(margin=ggplot2::margin(t=10)),
          axis.title.y = element_text(margin=ggplot2::margin(r=10)),
          plot.title = element_text(hjust=title.hjust))

    if(!is.null(threshold)) {
        p <- p + geom_hline(yintercept=threshold, linetype='dashed', col="grey40")
    }

    if (!is.null(title_size)) 
        p <- p + theme(plot.title = element_text(size=title_size))

    if (!is.null(ylim)) p <- p + ylim(ylim[1], ylim[2])

    if(print_p) print(p)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width, height=height, units='mm')
        ggsave(paste0(file, '.pdf'), p, width=width, height=height, units='mm')
    }
    return(p)
}
