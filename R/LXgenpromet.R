
LXgenpromet <- function(gene_meta,pro_meta){

#------安装必要的R包-----------------------------------------
installed_packs <- installed.packages()[,1]
packs <- c("openxlsx","dplyr","purrr","tidyverse","ggplot2")
not_packs <- packs[!packs %in% installed_packs]

if(length(not_packs)>0){
  packs_fun <- function(i){install.packages(i)}
  sapply(not_packs,packs_fun,simplify = T)
}  

lib_fun <- function(i){library(i,character.only = T)}
sapply(packs,lib_fun,simplify = T)

bio_pack <- c("limma")
not_bio <- bio_pack[!bio_pack %in% installed_packs]

if(length(not_bio)>0){
  bio_fun <- function(i){BiocManager::install(i)}
  sapply(not_bio,bio_fun,simplify = T)
}  

lib_fun <- function(i){library(i,character.only = T)}
sapply(bio_pack,lib_fun,simplify = T)

#------建立文件夹---------------
if(!dir.exists("analysis results"))
  dir.create("analysis results")

#------读取数据---------------------------------------------------

gene_meta_path <- read.csv(gene_meta)
gene_meta_kegg_pathways <- dplyr::filter(gene_meta_path,gene_meta_path$pathway_source=="KEGG")
gene_meta_kegg_pathways$pathway_name <- str_extract(gene_meta_kegg_pathways$pathway_name,".*(?= -)")
gene_meta_kegg_pathways <- dplyr::filter(gene_meta_kegg_pathways,num_overlapping_genes>0 & num_overlapping_metabolites>0)

pro_meta_path <- read.csv(pro_meta)
pro_meta_kegg_pathways <- dplyr::filter(pro_meta_path,pro_meta_path$pathway_source=="KEGG")
pro_meta_kegg_pathways$pathway_name <- str_extract(pro_meta_kegg_pathways$pathway_name,".*(?= -)")
pro_meta_kegg_pathways <- dplyr::filter(pro_meta_kegg_pathways,num_overlapping_genes>0 & num_overlapping_metabolites>0)
names(pro_meta_kegg_pathways) <- c("pathway_name","pathway_source","num_overlapping_proteins","overlapping_proteins",
                                   "num_all_pathway_proteins","P_proteins","Q_proteins","num_overlapping_metabolites",
                                   "overlapping_metabolites","num_all_pathway_metabolites","P_metabolites","Q_metabolites",
                                   "P_joint","Q_joint"  )

write.xlsx(gene_meta_kegg_pathways,"analysis results/gene_meta_kegg_pathways.xlsx")
write.xlsx(pro_meta_kegg_pathways,"analysis results/pro_meta_kegg_pathways.xlsx")

#---------kegg通路取交集------------------------------------------

gene_meta_kegg_pathways$pathway_name <- trimws(gene_meta_kegg_pathways$pathway_name)
pro_meta_kegg_pathways$pathway_name <- trimws(pro_meta_kegg_pathways$pathway_name)

venn <- inner_join(gene_meta_kegg_pathways,pro_meta_kegg_pathways,"pathway_name") 
names(venn)
venn_pathways <- filter(venn,P_joint.x<0.05 & P_joint.y<0.05)

write.xlsx(venn,"analysis results/gene_pro_meta_joint_pathways_all.xlsx")
write.xlsx(venn_pathways,"analysis results/gene_pro_meta_joint_pathways_p005.xlsx")

venn_path <- venn_pathways[,c(1,13,4,9,26,17,22)]
names(venn_path) <- c("pathways","P_gene_meta","overlapping_genes","overlapping_metabolites_x",
                      "P_protein_meta","overlapping_proteins","overlapping_metabolites_y")

gene_pro_meta_path <- data.frame(pathways=rep(venn_path$pathways,2),
                                 log2_p=c(-log2(venn_path$P_gene_meta),-log2(venn_path$P_protein_meta)),
                                 type=c(rep("genes-metabolites",nrow(venn_path)),rep("proteins-metabolites",nrow(venn_path)))
)

y_p <- max(gene_pro_meta_path$log2_p)

height_y <- y_p*1.2

nrow_path <- nrow(gene_pro_meta_path)/2

joint_title_size <- case_when(nrow_path>=30 ~12,
                              nrow_path>=20 ~12,
                              TRUE ~14)

joint_x_size <- case_when(nrow_path>=20 ~9,
                          nrow_path>=10 ~10,
                          TRUE ~12)

joint_y_size <- case_when(nrow_path>=20 ~12,
                          nrow_path>=10 ~12,
                          TRUE ~12)

joint_legend_size <- case_when(nrow_path>=30 ~12,
                               nrow_path>=20 ~12,
                               TRUE ~12)


bar_width <- case_when(nrow_path>=30 ~0.9,
                       nrow_path>=20 ~0.8,
                       TRUE ~0.7)

f1 <- ggplot(gene_pro_meta_path, aes(x = pathways, y = log2_p,fill=type))+
  geom_bar(position = "dodge",stat = "identity",width = bar_width)+
  scale_fill_manual(values=c("#008b8b","#f08080","#87ceeb"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(colour = "black", face="bold",size=12)) +
  labs(x="",y = "-log2(Pvalue)",title = 'Gene-protein-metabolite Joint Pathways')+
  scale_y_continuous(expand = c(0, 0),limits = c(0, height_y))

f1

log05 <- -log2(0.05)
log01 <- -log2(0.01)

line1 <- geom_hline(yintercept = c(log05),
                    linewidth = 0.6,
                    color = "blue",
                    lty = "dashed")
line2 <- geom_hline(yintercept = c(log01),
                    linewidth = 0.6,
                    color = "red",
                    lty = "dashed")

y1 <- geom_text(x=nrow(gene_pro_meta_path)/3,y=log05+0.8,label = c("p<0.05"),
                size=4,color="blue",fontface="italic")
y2 <- geom_text(x=nrow(gene_pro_meta_path)/3,y=log01+0.8,label = c("p<0.01"),
                size=4,color="blue",fontface="italic")

f2 <- f1+line1+line2+y1+y2

mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="black",face="bold",size =joint_title_size),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"),
        plot.margin = unit(c(t=0.5, r=0.5, b=0.5, l=2), "cm")
  )+
  theme(plot.title = element_text(hjust = 0.5))


xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=joint_x_size,angle =45,hjust=1))+
  theme(axis.text.y = element_text(face="bold",color="black",size=joint_y_size))

legend_theme <- theme(
  legend.title = element_blank(),
  legend.text = element_text(size = joint_legend_size, face = "bold"),
  legend.direction = "vertical",
  #legend.position = c(0.5,0.9),
  legend.background = element_blank()
)

f3 <- f2+mytheme+xytheme

ggsave("analysis results/Gene_protein_metabolite_Joint_pathways 01.png",
       f3,width=1200, height =1000, dpi=150,units = "px")

f4 <- f3+
  labs(fill="")+
  theme(legend.direction = "horizontal",
        legend.position = c(0.5,0.92),
        legend.text = element_text(size=14,face = "bold") )

ggsave("analysis results/Gene_protein_metabolite_Joint_pathways 02.png",
       f4,width=1200, height =1000, dpi=150,units = "px")

print("--------------------------------------------------------------")

print("The results can be found in the folder of <analysis results>")

f4

}









