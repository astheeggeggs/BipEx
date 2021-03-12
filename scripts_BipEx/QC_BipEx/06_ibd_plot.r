library(ggplot2)
library(ggsci)
library(dplyr)
library(latex2exp)
library(data.table)

# File locations and thresholds defined in r_options_BipEx.r
source("r_functions_and_parameters/r_options_BipEx.r")

# Watch out for O'Donovan!
IBD_FILE <- "gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/06_ibd.tsv"
df <- fread(IBD_FILE, sep='\t', stringsAsFactors=TRUE, header=TRUE, data.table=FALSE)

df$inferred_relationship <- 'Siblings'
df$inferred_relationship[df$ibd.PI_HAT <= IBD_THRESHOLD] <- 'Unrelated'
df$inferred_relationship[df$ibd.Z0 < 0.05 & df$ibd.Z1 < 0.05] <- 'Duplicate/Monozygotic twins'
df$inferred_relationship[df$ibd.Z0 < 0.05 & df$ibd.Z1 > 0.9] <- 'Parent-Offspring'

p <- ggplot(df, aes(x=ibd.Z0)) +
  geom_point(aes(y=ibd.Z1, color=inferred_relationship), size=0.5) +
  geom_abline(intercept=(2-2*IBD_THRESHOLD), slope=-2, linetype='dashed') +
  scale_color_d3('category20', limits=c('Unrelated', 'Siblings', 'Parent-Offspring', 'Duplicate/Monozygotic twins')) +
  labs(x='Proportion of loci with 0 shared alleles',
       y='Proportion of loci with 1 shared allele',
       title="",
       color='Inferred Relationship') +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
  theme_minimal() +
  theme(axis.title.x = element_text(margin=margin(t=10)),
        axis.title.y = element_text(margin=margin(r=10)),
        plot.title = element_text(hjust=0.5))
print(p)
ggsave(paste0(PLOTS, '06_IBD_plot', '.jpg'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '06_IBD_plot', '.pdf'), p, width=160, height=90, units='mm')

# # Create IBD plot similar to that used in Giulio's paper.
# dfpoly <- data.frame(x1=c(IBD_THRESHOLD, 0.5, IBD_THRESHOLD), y1=c(2*IBD_THRESHOLD, 1, 1), x2=c(0.5,1,1), y2=c(1,0,1))
# p <- ggplot(df, aes(x=ibd.PI_HAT)) +
#     geom_polygon(data=dfpoly, mapping=aes(x=x1, y=y1), fill='grey90') +
#     geom_polygon(data=dfpoly, mapping=aes(x=x2, y=y2), fill='grey90') +
#     geom_point(data=df, aes(y=ibd.Z1, color=inferred_relationship)) +
#     xlim(IBD_THRESHOLD,1) +
#     labs(x=TeX('Proportion IBD ($\\hat{\\pi}$)'), 
#          y='Proportion IBD1 (Z1)',
#          color='Inferred Relationship') +
#     theme_minimal() +
#     theme(axis.title.x = element_text(margin=margin(t=10)),
#         axis.title.y = element_text(margin=margin(r=10)),
#         plot.title = element_text(hjust=0.5))
# print(p)

# ggsave(paste0(PLOTS, '06_IBD_alt_plot', '.jpg'), p, width=160, height=90, units='mm')
# ggsave(paste0(PLOTS, '06_IBD_alt_plot', '.pdf'), p, width=160, height=90, units='mm')
