library(data.table)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)
library(latex2exp)
# install.packages('devtools')
# devtools::install_github('VPetukhov/ggrastr')
library(ggrastr)

chr_lengths_38 <- c(248956422, 242193529, 198295559,190214555,181538259,170805979,159345973,
    145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,
    83257441,80373285,58617616,64444167,46709983,50818468,156040895)

chr_lengths_37 <- c(249250621, 243199373,198022430,191154276,180915260,171115067,159138663,
    146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,
    81195210,78077248,59128983,63025520,48129895,51304566,155270560)

get_start_and_end <- function(chr_lengths) {
    start <- rep(0, length(chr_lengths))
    start[1] <- 1
    end <- rep(0, length(chr_lengths))
    end[1] <- chr_lengths[1]
    for(chr in 2:length(chr_lengths)) {
        start[chr] <- start[chr-1] + chr_lengths[chr-1]
        end[chr] <- end[chr-1] + chr_lengths[chr]
    }
    return(list(start=start, end=end))
}

make_manhattan_plot = function(contigs, positions, pvals, log_OR=NULL, labels=NULL, size_by_p=FALSE, buffer=100000000, title='', threshold=5, chr_lengths=chr_lengths_38,
  colour_aes=NULL, log_p_vals=FALSE, significance_T=5e-8, ggplot_theme=theme_bw, two_tone=TRUE, by_OR=FALSE, colour_1='#2b59a1', colour_2='#5fb756',
  save_figure=FALSE, file='file_out', scaling=1, width=230, height=100, title_size=NULL, minus_log_p_max=NULL)
{  
    contigs_ <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
    start_end <- get_start_and_end(chr_lengths)
    dt_contigs <- data.frame(contig=contigs_, start=start_end$start, end=start_end$end) %>%
    mutate(middle = floor(start + (end-start)/2),
           length = (end-start)) %>%
    mutate(shifted_position=middle + (contig - 1) * buffer)

    dt_plot <- data.frame(contig=contigs, position=as.integer(positions), pval=as.numeric(pvals), labels=labels) %>%
    mutate(x = dt_contigs[gsub('X', '23', contig), 'start'] + position + (as.integer(gsub('X', '23', contig))-1)*buffer)

    # Include two tone chromosome plotting
    if (two_tone) {
        dt_plot <- dt_plot %>% mutate(colour=ifelse((as.integer(gsub('X', '23', contig)) %% 2) == 0, colour_1, colour_2))
    } else {
        dt_plot <- dt_plot %>% mutate(colour=colour_1) 
    }

    if (by_OR & !is.null(log_OR)) {
        dt_plot$log_OR <- log_OR
        dt_plot <- dt_plot %>% mutate(colour = case_when(
            dt_plot$log_OR < -1 ~ "blue3",
            ((dt_plot$log_OR > -1) & (dt_plot$log_OR <= 0)) ~ "cornflowerblue",
            ((dt_plot$log_OR > 0) & (dt_plot$log_OR <= 1)) ~ "indianred3",
            dt_plot$log_OR > 1 ~ "red",
            TRUE ~ "grey40")
        )
        dt_plot$colour <- factor(dt_plot$colour, levels = c("red", "indianred3", "cornflowerblue", "blue3"))
    }

    # Were log p-values passed?
    if(!log_p_vals) {
        dt_plot <- dt_plot %>% mutate(y = -log10(pval)) %>% select(x, y, colour, labels)
    } else {
        dt_plot <- dt_plot %>% mutate(y= pval) %>% select(x, y, colour, labels)
    }

    if (size_by_p) {
        dt_plot <- dt_plot %>% mutate(size = case_when(
            dt_plot$y > 3 ~ 20,
            ((dt_plot$y > -log10(0.05)) & (dt_plot$y <= 3)) ~ 5,
            TRUE ~ 0.5)
        )
    }

    if (size_by_p) {
        p <- ggplot(dt_plot, aes(x=x,y=y,label=labels, col=colour, size=size)) + geom_point_rast()
    } else {
        p <- ggplot(dt_plot, aes(x=x,y=y,label=labels, col=colour)) + geom_point_rast(size=0.5)
    }

    p <- p + geom_hline(yintercept=-log10(significance_T), color='#E15759', linetype='dashed') +
        scale_x_continuous(breaks=dt_contigs$shifted_position, labels=gsub(23, 'X', dt_contigs$contig)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
        labs(x='Chromosome', y=TeX("$-\\log_{10}(\\mathit{p})$"), title=title) + ggplot_theme()

    if ((by_OR == FALSE) & (size_by_p==FALSE)) {
        p <- p + theme(legend.position = "none")
    }

    if (two_tone) { p <- p + scale_color_manual(values=c(colour_1,colour_2)) }
    if (by_OR) { 
        p <- p + scale_colour_manual(name="Odds ratio",
            values = levels(dt_plot$colour), labels=c("> 10", "1 - 10", "0.1 - 1", "< 0.1"))
    }

    if (!is.null(labels))
        p <- p + geom_label_repel(
            data=subset(dt_plot, y > threshold), size = 5,
            aes(label=labels), color='grey30', box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50')

    if (!is.null(title_size)) 
        p <- p + theme(plot.title = element_text(size=title_size))

    if (!is.null(minus_log_p_max))
        p <- p + ylim(0, minus_log_p_max)

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
    }

    return(list(p=p, dt=dt_plot))
}
