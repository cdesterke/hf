x<-read.table("metaHF.tsv",h=T,sep="\t")

# Liste des groupes uniques
groupes <- unique(x$annot)

# Créer un data frame pour stocker les résultats
resultats_df <- data.frame(
  Comparaison = character(),
  Mediane_Groupe1 = numeric(),
  Mediane_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Boucle pour faire les tests de Wilcoxon pour chaque paire de groupes
for (i in 1:(length(groupes) - 1)) {
  for (j in (i + 1):length(groupes)) {
    # Sélection des deux groupes à comparer
    groupe1 <- x$cardiomyocyte1[x$annot == groupes[i]]
    groupe2 <- x$cardiomyocyte1[x$annot == groupes[j]]
    
    # Test de Wilcoxon
    test <- wilcox.test(groupe1, groupe2)
    
    # Ajouter les résultats dans le data frame
    resultats_df <- rbind(resultats_df, data.frame(
      Comparaison = paste(groupes[i], "vs", groupes[j]),
      Mediane_Groupe1 = median(groupe1, na.rm = TRUE),
      Mediane_Groupe2 = median(groupe2, na.rm = TRUE),
      p_value = test$p.value
    ))
  }
}

# Afficher les résultats
print(resultats_df)

# Afficher les résultats sous forme de data frame
print(resultats_df)
##write.csv(resultats_df,file="comparisons_cellsgroups_ferroptosisUP.csv",row.names=F)

cardio<-resultats_df[grepl("Cardiomyocytes",resultats_df$Comparaison),]

write.table(cardio,file="Table_cardiomyocytes_celltypes.tsv",row.names=F,sep="\t")





#### local cardio

# Liste des groupes uniques
groupes <- unique(x$Type)

# Créer un data frame pour stocker les résultats
resultats_df <- data.frame(
  Comparaison = character(),
  Mediane_Groupe1 = numeric(),
  Mediane_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Boucle pour faire les tests de Wilcoxon pour chaque paire de groupes
for (i in 1:(length(groupes) - 1)) {
  for (j in (i + 1):length(groupes)) {
    # Sélection des deux groupes à comparer
    groupe1 <- x$cardiomyocyte1[x$Type == groupes[i]]
    groupe2 <- x$cardiomyocyte1[x$Type == groupes[j]]
    
    # Test de Wilcoxon
    test <- wilcox.test(groupe1, groupe2)
    
    # Ajouter les résultats dans le data frame
    resultats_df <- rbind(resultats_df, data.frame(
      Comparaison = paste(groupes[i], "vs", groupes[j]),
      Mediane_Groupe1 = median(groupe1, na.rm = TRUE),
      Mediane_Groupe2 = median(groupe2, na.rm = TRUE),
      p_value = test$p.value
    ))
  }
}

# Afficher les résultats
print(resultats_df)

# Afficher les résultats sous forme de data frame
print(resultats_df)
write.csv(resultats_df,file="Table_cardiomyocyte_types.csv",row.names=F)

cardio<-resultats_df[grepl("Cardiomyocytes",resultats_df$Comparaison),]

write.table(cardio,file="Table_cardiomyocytes_celltypes.tsv",row.names=F,sep="\t")


### fibroblast cell type


# Liste des groupes uniques
groupes <- unique(x$annot)

# Créer un data frame pour stocker les résultats
resultats_df <- data.frame(
  Comparaison = character(),
  Mediane_Groupe1 = numeric(),
  Mediane_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Boucle pour faire les tests de Wilcoxon pour chaque paire de groupes
for (i in 1:(length(groupes) - 1)) {
  for (j in (i + 1):length(groupes)) {
    # Sélection des deux groupes à comparer
    groupe1 <- x$fibroblast1[x$annot == groupes[i]]
    groupe2 <- x$fibroblast1[x$annot == groupes[j]]
    
    # Test de Wilcoxon
    test <- wilcox.test(groupe1, groupe2)
    
    # Ajouter les résultats dans le data frame
    resultats_df <- rbind(resultats_df, data.frame(
      Comparaison = paste(groupes[i], "vs", groupes[j]),
      Mediane_Groupe1 = median(groupe1, na.rm = TRUE),
      Mediane_Groupe2 = median(groupe2, na.rm = TRUE),
      p_value = test$p.value
    ))
  }
}

# Afficher les résultats
print(resultats_df)

# Afficher les résultats sous forme de data frame
print(resultats_df)
##write.csv(resultats_df,file="comparisons_cellsgroups_ferroptosisUP.csv",row.names=F)

fibro<-resultats_df[grepl("Fibroblasts",resultats_df$Comparaison),]

write.table(fibro,file="Table_fibroblast_celltypes.tsv",row.names=F,sep="\t")


## local fibro

# Liste des groupes uniques
groupes <- unique(x$Type)

# Créer un data frame pour stocker les résultats
resultats_df <- data.frame(
  Comparaison = character(),
  Mediane_Groupe1 = numeric(),
  Mediane_Groupe2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Boucle pour faire les tests de Wilcoxon pour chaque paire de groupes
for (i in 1:(length(groupes) - 1)) {
  for (j in (i + 1):length(groupes)) {
    # Sélection des deux groupes à comparer
    groupe1 <- x$fibroblast1[x$Type == groupes[i]]
    groupe2 <- x$fibroblast1[x$Type == groupes[j]]
    
    # Test de Wilcoxon
    test <- wilcox.test(groupe1, groupe2)
    
    # Ajouter les résultats dans le data frame
    resultats_df <- rbind(resultats_df, data.frame(
      Comparaison = paste(groupes[i], "vs", groupes[j]),
      Mediane_Groupe1 = median(groupe1, na.rm = TRUE),
      Mediane_Groupe2 = median(groupe2, na.rm = TRUE),
      p_value = test$p.value
    ))
  }
}

# Afficher les résultats
print(resultats_df)

# Afficher les résultats sous forme de data frame
print(resultats_df)
write.csv(resultats_df,file="Table_fibroblast_types.csv",row.names=F)


