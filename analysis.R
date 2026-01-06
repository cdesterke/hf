list.files()

load("HeartFailure_Voom_Combat_ESET.rda")

library(Biobase)

data<-as.data.frame(exprs(eset))

pheno<-pData(eset)

all(colnames(data)==row.names(pheno))


df<-read.csv("fdrtx.csv",h=T)

vector<-df$Association.Name

inter<-intersect(vector,row.names(data))

small<-data[row.names(data)%in%inter,]

save(small,file="selection263.rda")


library(transpipe15)
pcatrans(small,pheno,group="group",names=F,alpha=1,pal="Set1")

annot<-pheno
annot$tissue<-NULL
annot$disease_detail<-NULL
annot$age<-NULL


bestheat(small,annot,scale="row",rownames=F,font=10)


### kmeans


trans<-as.data.frame(t(small))

tot.withinss <- vector(mode="character", length=10)
for (i in 1:10){
  Cluster <- kmeans(trans, center=i, nstart=20)
  tot.withinss[i] <- Cluster$tot.withinss
}
Let’s visualize it.

plot(1:10, tot.withinss, type="b", pch=19)

set.seed(101)
Cluster <- kmeans(trans, center=2, nstart=20)
Cluster

prop.table(table(pheno$group,Cluster$cluster))*100
table(pheno$group,Cluster$cluster)
library(cluster)
clusplot(trans, Cluster$cluster, color=T, shade=T, labels=0, lines=0)

pheno$clusters<-Cluster$cluster


pcatrans(small,Cluster,group="cluster",pal="Set1",alpha=0.5,names=F)


pheno<-cbind(pheno,trans)

pheno%>%mutate(level=ifelse(group=="HF","1","2"))->pheno
pheno$level<-as.factor(pheno$level)
pheno$clusters<-as.factor(pheno$clusters)
library(caret)
cm <- confusionMatrix(data = pheno$clusters, reference = pheno$level)

pheno$consensus<-as.factor(pheno$consensus)

table(pheno$clusters,pheno$level)

draw_confusion_matrix(cm)


##pamr

library(pamr)
complete <- small[complete.cases(small), ]
complete<-as.matrix(complete)
d <- list(x=complete,y=pheno$group, geneid=as.character(row.names(complete)), genenames=paste(row.names(complete)))
pamr.menu(d)






# Entraîner le modèle
pamr_model <- pamr.train(d)

# Effectuer une validation croisée
pamr_cv <- pamr.cv(pamr_model, d)

# Afficher les résultats
print(pamr_cv)
pamr.plotcv(pamr_cv)



# Extraire la liste des prédicteurs avec leurs valeurs de FDR
pamr_thresholds <- pamr.listgenes(pamr_model, d, threshold=1.04)


# Calcul du False Discovery Rate (FDR)
myfdr <- pamr.fdr(pamr_model, d)



# Construire le dataframe des résultats FDR
fdr_df <- data.frame(
  Threshold = myfdr$results[, "Threshold"],
  Median_FDR = myfdr$results[, "Median FDR"],
  Percentile_90_FDR = myfdr$results[, "90th percentile of FDR"],
  Number_of_Genes = myfdr$results[, "Number of significant genes"]
)

# Définir le seuil optimal basé sur FDR < 0.05
selected_threshold <- fdr_df$Threshold[min(which(fdr_df$Median_FDR < 0.05))]
selected_genes <- fdr_df$Number_of_Genes[min(which(fdr_df$Median_FDR < 0.05))]

# Graphique avec seuil FDR sélectionné et annotation du nombre de gènes
ggplot(fdr_df, aes(x = Threshold)) +
  geom_line(aes(y = Median_FDR, color = "Median FDR"), size = 1.2, linetype = "solid") +
  geom_line(aes(y = Percentile_90_FDR, color = "90e percentile FDR"), size = .5 ) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = selected_threshold, linetype = "dashed", color = "purple", size = 1) +
  annotate("text", x = selected_threshold, y = 0.1, 
           label = paste0("Selected threshold = ", round(selected_threshold, 2)), color = "purple", vjust = -0.5,hjust=-0.5) +
  annotate("text", x = selected_threshold, y = 0.15, 
           label = paste0("Significant genes = ", selected_genes), color = "purple", vjust = -0.5, hjust=-0.5) +
  scale_color_manual(values = c("Median FDR" = "blue", "90e percentile FDR" = "red")) +
  theme_minimal() +
  labs(title = "False Discovery Rate (FDR) according threshold",
       x = "Threshold",
       y = "Estimate FDR",
       color = "Type of FDR")


## pamr matriceconfusion a partir du seuil définit par FDR


library(caret)



# Définir le seuil optimal basé sur FDR < 0.05
selected_threshold <- myfdr$results[min(which(myfdr$results[, "Median FDR"] < 0.05)), "Threshold"]

# Prédiction des classes avec le seuil optimal
predictions <- pamr.predict(pamr_model, d$x, threshold = selected_threshold)

# Création de la matrice de confusion
conf_matrix <- confusionMatrix(factor(predictions), factor(d$y))

# Extraction des données pour la visualisation
conf_df <- as.data.frame(conf_matrix$table)

# error rate
overall_error_rate <- 1 - sum(diag(conf_matrix$table)) / sum(conf_matrix$table)



# Graphique de la matrice de confusion
ggplot(conf_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = paste0("Confusion matrix FDR<=0.05 | Overall Error Rate: ", round(overall_error_rate, 3)),
       x = "Real classes",
       y = "Class prediction") +
  geom_text(aes(label = Freq), color = "black", size = 5)


# pamr centroids barplot groups


library(ggplot2)
library(dplyr)
library(tidyr)
# Création du dataframe
data <- read.csv("pamr.csv",h=T)
colnames(data)<-c("rank","id","C1","C2","C3")
# Transformation en format long pour ggplot2
data_long <- pivot_longer(data, cols = c("C1", "C2", "C3"), 
                          names_to = "Clusters", values_to = "Score")



# Transformation des données en format long et tri par 'rank' (ordre décroissant)
data_long <- data_long %>%
  arrange(desc(rank))  # Trie les Hallmarks en ordre décroissant selon 'rank'

# Création du barplot
ggplot(data_long, aes(x = reorder(id, desc(rank)), y = Score, fill = Clusters)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +  
  labs(x = "Hallmarks", y = "Pamr predictive scores", fill = "Clusters") +
  coord_flip() +
  theme_minimal()




draw_confusion_matrix <- function(cm) {

  total <- sum(cm$table)
  res <- as.numeric(cm$table)

  # Generate color gradients. Palettes come from RColorBrewer.
  greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
  redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
  getColor <- function (greenOrRed = "green", amount = 0) {
    if (amount == 0)
      return("#FFFFFF")
    palette <- greenPalette
    if (greenOrRed == "red")
      palette <- redPalette
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
  }

  # set the basic layout
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)

  # create the matrix 
  classes = colnames(cm$table)
  rect(150, 430, 240, 370, col=getColor("green", res[1]))
  text(195, 435, classes[1], cex=1.2)
  rect(250, 430, 340, 370, col=getColor("red", res[3]))
  text(295, 435, classes[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=getColor("red", res[2]))
  rect(250, 305, 340, 365, col=getColor("green", res[4]))
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)

  # add in the cm results
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')

  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)

  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}


