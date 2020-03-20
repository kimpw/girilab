#' Regional Plot
#'
#' @param g genedata
#' @param m map data
#' @param t tissue data
#' @param chr chromosome
#' @param tag_p where tag above p
#' @param y_lab label y axis
#' @param name_gene gene name select
#' @param y_min min number tick
#' @param y_max max number tick
#' @param y_ticks how many tick
#' @param cut_ctr zoom in/out
#' @param sig_line1 significant top number
#' @param sig_line2 significan bottom number
#' @param color_tissue color choice
#'
#' @return plot
#' @export
#'
#' @import dplyr
#' @import ggplot2
#'
#'
spec_chr_plot<-function(g,m,t,chr=8,tag_p=3,y_lab="p",name_gene="PENK",y_min=NULL,
                        y_max=NULL,y_ticks=NULL,cut_ctr=2,sig_line1=3,sig_line2=-3,color_tissue=NULL){
  m<-m%>%
    rename(gene=ENSG)
  m$CHR<-as.integer(m$CHR)
  g<-g%>%
    group_by(CHR)%>%
    select(-SNP)%>%
    mutate(Tis="Single SNP",
           BP=BP/1000000,
           P=log10(P))%>%
    filter(CHR!=23&CHR==chr)
  postart<-t%>%
    left_join(m,by='gene')%>%
    mutate(P=-log10(pvalue),START_POS=START_POS/1000000,END_POS=END_POS/1000000,gene_type=paste(gene_name,tissue,sep="_"))%>%
    filter(CHR==chr)%>%
    select(CHR,gene_name,tissue,START_POS,P,gene_type)%>%
    rename(BP=START_POS)
  poend<-t%>%
    left_join(m,by='gene')%>%
    mutate(P=-log10(pvalue),START_POS=START_POS/1000000,END_POS=END_POS/1000000,gene_type=paste(gene_name,tissue,sep="_"))%>%
    filter(CHR==chr)%>%
    select(CHR,gene_name,tissue,END_POS,P,gene_type)%>%
    rename(BP=END_POS)
  popopo<-bind_rows(postart,poend)
  xyn<-as.numeric(popopo[popopo$gene_name==name_gene,"BP"][1,1])
  popopo<-popopo%>%
    filter(BP>xyn-cut_ctr&BP<xyn+cut_ctr)
  popopo<-bind_rows(postart,poend)%>%
    filter(BP>xyn-cut_ctr&BP<xyn+cut_ctr)
  potest<-g%>%
    filter(BP>xyn-cut_ctr&BP<xyn+cut_ctr)
  for_tag<-popopo%>%
    filter(P>tag_p)
  if(is.null(color_tissue)==TRUE) {
    color_tissue= c("bisque3","peachpuff2","lightgreen","firebrick4","red",
                    "orangered","blue4","darkturquoise","dodgerblue","cyan","lightsteelblue",
                      "lightslateblue","skyblue1", "royalblue","slateblue3",
                    "deepskyblue1", "lightpink1","gray32","grey85","goldenrod4","khaki3",
                      "lightgoldenrod","tan2","yellow","tomato","tomato3", "burlywood4",
                    "plum","chocolate4","paleturquoise4","palevioletred3","seagreen1",
                      "darkgreen","honeydew4", "cornsilk2","wheat2","sandybrown","lightsalmon3",
                    "gold","darkslategrey","forestgreen","lightcoral", "deeppink","gray0")
  }
  else {
    color_tissue = color_tissue
  }
  labels_cat <- c(sort(unique(as.character(t$tissue))))
  if(is.null(y_min) == TRUE) {y_min <- round(min(potest$P, na.rm=TRUE) - 1)} else {y_min <- y_min}
  if(is.null(y_max) == TRUE) {y_max <- round(max(popopo$P, na.rm=TRUE) + 1)} else {y_max <- y_max}
  if(is.null(y_ticks) == TRUE) {
    if((y_max - y_min) > 16) {y_ticks <- 16} else {y_ticks <- round(y_max - y_min)}
  } else {y_ticks <- y_ticks}
  break_length <- round(round(y_max-y_min)/y_ticks)
  genom<-ggplot()+
    theme_bw()+
    theme(axis.text.x = element_text(size=12, color = 'black'),
          axis.text.y = element_text(size=12, color = 'black'),
          axis.title.x = element_text(size = 12, face = "bold", color ="black"),
          axis.title.y = element_text(size = 12, face = "bold", color ="black"),
          axis.ticks.x=element_line())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab(paste("Chromosomes",chr,"(MB)", sep=" ")) +
    ylab(y_lab)+
    ggrepel::geom_label_repel(data = for_tag, aes(x = BP, y = P, label = gene_name))+
    geom_point(data = potest, aes(x = BP, y = P)) +
    geom_line(data = popopo, aes(x = BP, y = P, group=gene_type,size=1.5,color = factor(tissue)))+
    scale_colour_manual(name = "Tissue Type", values = c("black", color_tissue), labels = c("Single SNP", labels_cat))+
    guides(shape = "none", size = "none", colour = guide_legend(reverse = TRUE, override.aes = list(size=6))) +
    geom_hline(aes(yintercept = 0), size = 3) +
    geom_hline(yintercept = sig_line1, size = .5)+
    geom_hline(yintercept = sig_line2, size = .5)+
    scale_y_continuous(breaks=seq(y_min, y_max, break_length)) +
    expand_limits(y=c(y_min, y_max))
  return(genom)
}
