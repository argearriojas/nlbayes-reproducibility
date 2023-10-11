library(ggplot2)
library(ggrepel)
library(cowplot)


make.figure.4 <- function() {

    get.scatter <- function(experiment.name) {
        df <- read.csv(paste0('data/oe_',experiment.name,'_on_net_regulonbrca_nlbayes_and_viper.csv'))
        xmax <- max(- log10(df$viper.pvalue[!is.na(df$viper.pvalue)]))
    
    ggplot(df, aes(
            x=-log10(viper.pvalue),
            y=posterior.p,
            label=ifelse((posterior.p > 0.2 & viper.pvalue < 0.05) | rank(viper.pvalue) <= 3 | rank <= 3, symbol, '')
        )) +
        theme_bw(base_size = 14) +
        geom_vline(xintercept = -log10(c(0.05))) + 
        geom_hline(yintercept = c(0.2)) +
        geom_point(
            color=ifelse((df$posterior.p >= 0.2 & df$viper.pvalue <= 0.05), 'red', ifelse((df$posterior.p < 0.2 & df$viper.pvalue > 0.05), 'gray', 'black')),
            na.rm = TRUE
        ) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10, face="bold")) +
        theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face="bold")) +
        theme(strip.background = element_rect(colour="black", fill="white",
                                            size=0.5, linetype="solid")) +
        theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
        
        coord_cartesian(xlim = c(NA, 0.5 + xmax ), ylim = c(NA, 1.2 )) + 
        geom_text_repel(
            color=ifelse((df$posterior.p >= 0.2 & df$viper.pvalue <= 0.05), 'red', ifelse((df$posterior.p < 0.2 & df$viper.pvalue > 0.05), 'gray', 'black')),
            force=10,
            force_pull=0.1,
            size = 2.5, fontface = "bold",
            box.padding = unit(0.6, "lines"),
            max.overlaps = 300,
            nudge_x = 0.25, 
            nudge_y = 0.25,
            hjust=0.25,
            segment.size = 0.1,
            na.rm = TRUE
        ) +
        labs(title = paste(paste0(toupper(experiment.name))), y="TF posterior mean", x="- log10( viper p-value )") + 
        theme(
            plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
            axis.title.x = element_text(size=12, face="bold", hjust = 0.5),
            axis.title.y = element_text(size=12, face="bold")
        )
    }

    p1 <- get.scatter('e2f3')
    p2 <- get.scatter('myc')
    p3 <- get.scatter('ras')

    plot_title <- ggdraw() + 
        draw_label("Results comparison between NLBayes and VIPER", fontface = 'bold', x = 0, hjust = 0) +
        theme(plot.margin = margin(0, 0, 0, 7))
        
    plot_row <- plot_grid(p1, p2, p3, labels = c('A)', 'B)', 'C)'), nrow = 1)
    pg <- plot_grid(plot_title, plot_row, ncol = 1, rel_heights = c(0.1, 1))

    if (!dir.exists('figures')) {
        dir.create('figures')
    }
    ggsave2('figures/fig4.png', pg, height = 3, width = 10)
    pg
}
