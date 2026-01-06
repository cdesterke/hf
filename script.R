list.files()

data<-read.csv("GSE121893_human_heart_sc_umi.csv",h=T,row.names=1)

library(data.table)

library(dplyr)
meta%>%inner_join(
meta<-fread("GSE121893_human_heart_sc_info.txt")





clusters<-fread("GSE121893_all_heart_cell_cluster_info.txt")



meta<-as.data.frame(meta)

row.names(meta)<-meta$ID
all(colnames(data)==row.names(meta))
library(Seurat)

set  <- CreateSeuratObject(counts = data, min.cells = 50, min.features=200, project = "hf")

set <- AddMetaData(object = set,  metadata = meta)

save(set,file="set.rda")
Idents(set)<-"Type"



set <- subset(set, subset = nFeature_RNA < 5000 & nCount_RNA < 20000)
VlnPlot(set, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

table(set$X384.Well.Plate.Location)

set[["RNA"]]<-split(set[["RNA"]],set$Individual)

all<-set



all <- NormalizeData(all)
all <- FindVariableFeatures(all)
all <- ScaleData(all)
all <- RunPCA(all)
ElbowPlot(all)



all <- IntegrateLayers(
  object = all, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

# re-join layers after integration
all[["RNA"]] <- JoinLayers(all[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)
Idents(all)<-"Individual"
hf<-subset(all,idents=c("C1","C2","D1","D2","D4","D5"))




##hf<-subset(hf,idents=c("C1","C2","D1","D2","D4","D5"))
hf <- FindNeighbors(hf, dims = 1:10, reduction = "harmony")
hf <- FindClusters(hf, resolution = 0.5, cluster.name = "clusters")
hf <- RunUMAP(hf, dims = 1:10, reduction = "harmony", reduction.name = "UMAP")
hf <- RunTSNE(hf, dims = 1:10, reduction = "harmony", reduction.name = "TSNE")

split.by="Type",
library(pals)
DimPlot(hf, reduction = "harmony", group.by ="CellType",split.by="Individual",pt.size=1,cols=glasbey())
DimPlot(hf, reduction = "TSNE",group.by="Type",cols=glasbey(),pt.size=1)

save(hf,file="harmonyHF.rda")




library(ggplot2)
colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=12)
FeaturePlot(hf,reduction="TSNE",pt.size=0.1,features=c("LUM"),
	min.cutoff = "q9",col=c("grey90","darkred"))

Idents(hf)<-"clusters"

markers <- FindAllMarkers(hf, only.pos = TRUE)

library(dplyr)
markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)->markers2


write.table(markers2,file="markers.tsv",row.names=F,sep="\t")


library(RColorBrewer)
library(ggplot2)
colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=16)
FeaturePlot(data,reduction="umap",pt.size=0.1,features="ferroptosis.down1",min.cutoff = "q9",split.by="sample_group",col=c("grey90","darkred"))+
	plotTheme+coord_fixed()



library(dplyr)

meta<-hf[[]]
meta%>%mutate(annot=case_when(clusters=="0"~"Cardiomyocytes",
					clusters=="1"~"Endothelial-cells",
					clusters=="2"~"Endothelial-cells",
					clusters=="3"~"Pericytes",
					clusters=="4"~"Cardiomyocytes",
					clusters=="5"~"Fibroblasts",
					clusters=="6"~"Endothelial-cells",
					clusters=="7"~"Macrophages",
					clusters=="8"~"NKT-cells",
					clusters=="9"~"Epicardial-cells"))->meta	

hf$annot<-meta$annot

vector<-c("MYH7","TTN","HSPB7","TNNI3","VWF","FABP4","AQP1","MCF2L","PRG4","IL18","MSLN","KRT18","DCN","LUM","COL1A2","ACKR1","POSTN","LIFR",
"MRC1","CSF1R","CD163","CD8A","NKG7","IL7R","CCL5","PDGFRB","RGS5","CSPG4")
DotPlot(hf,group.by="annot",features=vector)+
coord_flip() + scale_color_gradientn(colors = colgex) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))


library(pals)
DimPlot(hf, reduction = "TSNE", group.by ="annot",pt.size=1,cols=glasbey())




library(chi2residuals)

sub <- preProcess(meta, v1 = "Type", v2 = "annot")
head(sub,n=10)

residuals <- computeResiduals(sub, col1 = "Type", col2 = "annot")


plotResiduals(residuals,
              col1 = "annot",
              col2 = "Type",
              themeSize = 20,
              labelSize = 4,
              colorLow = "turquoise",
              colorHigh = "purple",
              colorLabels = "white",
              title = "Significant residuals p<0.05")


plotNet(residuals,
        var1 = "annot",
        var2 = "Type",
        node_colors = c(annot = "green", Type = "skyblue"),
        edge_colors = c(positive = "red", negative = "blue", nonsignificant = "lightgrey"))




hub<-c("ACE","BCL2","DDR1","LUM","MFAP4","MPV17","NOX4","NPPB","SERPINH1","SIRT6","TGFB3","TNNC1","AGT","LOXL2",
"MYH7","NPPA","TNNI3","ANKRD1","COL1A2","EDNRA","MMP2","SLC9A1","COL14A1","COL1A1","COL3A1",
"MYOCD","POSTN","LOX","TGFB1","TGFB2")


DotPlot(hf,group.by="annot",features=hub)+
scale_color_gradientn(colors = colgex) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))



FeaturePlot(hf,reduction="TSNE",pt.size=0.5,features="AGT",min.cutoff = "q9",col=c("grey90","darkred"))









cardiomyocyte <- list(c("NPPB","TNNC1","MYH7","NPPA","TNNI3","ANKRD1"))
hf<- AddModuleScore(
  object = hf,
  features = cardiomyocyte, name = 'cardiomyocyte',
ctrl=6
)
Idents(hf)<-"annot"
library(pals)
VlnPlot(hf, features = c("cardiomyocyte1"), slot = "data", log = TRUE,pt.size=1,split.by="annot",col=cols25())


Idents(hf)<-"Type"
VlnPlot(
  hf,
  features = c("cardiomyocyte1"),
  slot = "data",
  log = TRUE,
  pt.size = 1,
  split.by = "Type",
  col = cols25()
) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # labels verticaux
    legend.position = "none"                                        # retirer la légende
  )

vector_cardio<-c("NPPB","TNNC1","MYH7","NPPA","TNNI3","ANKRD1")
DotPlot(hf,group.by="annot",features=vector_cardio)+
scale_color_gradientn(colors = colgex) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))



colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=16)
FeaturePlot(hf,reduction="TSNE",pt.size=0.1,features="cardiomyocyte1",min.cutoff = "q9",col=c("grey90","darkred"))+
	plotTheme+coord_fixed()

DotPlot(data,group.by="Cell.class",features=vector)+
coord_flip() + scale_color_gradientn(colors = colgex) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))


## fibro

fibro <- list(c("LUM","COL1A1","COL14A1","MMP2","COL3A1","MFAP4"))
hf<- AddModuleScore(
  object = hf,
  features = fibro, name = 'fibroblast',
ctrl=6
)
Idents(hf)<-"annot"
library(pals)
VlnPlot(hf, features = c("fibroblast1"), slot = "data", log = TRUE,pt.size=1,split.by="annot",col=cols25())

vector_fibro<-c("LUM","COL1A1","COL14A1","MMP2","COL3A1","MFAP4")
Idents(hf)<-"Type"
VlnPlot(
  hf,
  features = c("fibroblast1"),
  slot = "data",
  log = TRUE,
  pt.size = 1,
  split.by = "Type",
  col = cols25()
) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # labels verticaux
    legend.position = "none"                                        # retirer la légende
  )


DotPlot(hf,group.by="Type",features=vector_fibro)+
scale_color_gradientn(colors = colgex) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))


colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=16)
FeaturePlot(hf,reduction="TSNE",pt.size=0.5,features="fibroblast1",min.cutoff = "q9",col=c("grey90","darkred"))+
	plotTheme+coord_fixed()

meta<-hf[[]]

write.table(meta,file="metaHF.tsv",row.names=T,sep="\t")

library(pheatmap)


table(set$CellType,set$Individual)

pheatmap(table(set$CellType,set$Individual),scale="none",cluster_cols=F)
