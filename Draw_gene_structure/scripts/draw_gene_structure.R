library(ggplot2)
library(ggtranscript)
library(dplyr)
library(optparse)


option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="enter a processed gff file"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="enter output picture name"),
  make_option(c("-f", "--format"), type = "character", default = "png", 
              help="select picture format, default:png, {png, jpg, pdf, svg, pptx}")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

gtf = read.delim(args$input, header = F, quote = "")

colnames(gtf) = c("Chr","source","type","start","end","score","strand","phase","transcript_name","gene_name","domain")
gtf$transcript_name=as.factor(gtf$transcript_name)

colors_list <- c("#B9C9DE", "#494D6B", "#C2ADAA", "#BDD3EE", "#F2F5F7", "#83AFDA")

domain_count <- length(unique(gtf$domain[gtf$domain != "."]))
domain_colors <- colors_list[1:domain_count]

get_heigth_width <- function(gtf) {
  basis_heigth <- 0.8
  heigth_factor <- 0.5
  basis_width <- 2
  width_factor <- 1
  mrna_data <- filter(gtf, type=="mRNA")
  mrna_max_length <- max(mrna_data$end - mrna_data$start)
  mrna_count <- length(rownames(mrna_data))
  plot_width <- basis_width + (mrna_max_length / 1000)* width_factor
  plot_heigth <- basis_heigth + mrna_count * heigth_factor
  return(list(plot_width, plot_heigth))
}


gtf_draw =
  ggplot(filter(gtf, type=="exon"),aes(xstart=start,xend=end,y=transcript_name)) +
  ggtranscript::geom_range(data = filter(gtf, type == "exon"),
                           fill = "white", , color = "grey30", height = 0.25, linewidth = 0.7) +
  ggtranscript::geom_range(data = filter(gtf, type=="CDS"), 
                           color = "grey30", linewidth = 0.7) +
  ggtranscript::geom_range(data = filter(gtf, type=="domain"), aes(fill = domain),
                           color = "black", linewidth = 0.7) +
  scale_fill_manual(values = domain_colors) +
  geom_intron(data = to_intron(filter(gtf,type=="exon"),"transcript_name"),
              aes(strand=strand), arrow.min.intron.length=10,
              arrow = grid::arrow(ends = "last", length = unit(1.5, "mm"))) +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 12,color = "black"),
        legend.position = "top",
        legend.box.spacing = unit(0.2, "lines"),
        legend.title = element_blank())

size <- get_heigth_width(gtf)
plot_width <- size[[1]]
plot_heigth <- size[[2]]

filename = paste0(args$output, ".", args$format)
if (args$format == "pptx") {
  library(eoffice)
  topptx(gtf_draw, filename, width = plot_width, height = plot_heigth, units = "in")
} else {
  ggsave(filename =  paste0(args$output, ".", args$format), gtf_draw,  width = plot_width, height = plot_heigth, units = "in",
        dpi = 600, device = args$format)
}
