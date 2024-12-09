# Normalisation des matrices de transition
transition_matrix_2019 <- sweep(matrix2019, 1, rowSums(matrix2019), "/")
transition_matrix_2019[is.na(transition_matrix_2019)] <- 0

transition_matrix_2021 <- sweep(matrix2021, 1, rowSums(matrix2021), "/")
transition_matrix_2021[is.na(transition_matrix_2021)] <- 0

# Calcul des distributions stationnaires en forçant les valeurs réelles
stationary_2019 <- Re(eigen(t(transition_matrix_2019))$vectors[, 1])
stationary_2021 <- Re(eigen(t(transition_matrix_2021))$vectors[, 1])

# Normalisation des distributions stationnaires
stationary_2019 <- stationary_2019 / sum(stationary_2019)
stationary_2021 <- stationary_2021 / sum(stationary_2021)

# Affichage des distributions stationnaires
print("Distribution stationnaire 2019:")
print(stationary_2019)

print("Distribution stationnaire 2021:")
print(stationary_2021)

barplot(cbind(stationary_2019, stationary_2021), beside = TRUE, 
        col = c("blue", "red"), 
        names.arg = state_names,  # Remplacer par les noms des états
        legend.text = c("2019", "2021"),
        main = "Comparaison des distributions stationnaires (2019 vs 2021)",
        xlab = "État", ylab = "Probabilité")

dist_matrix_2019 <- dist(matrix2019)  # Calcul des distances pour 2019
dist_matrix_2021 <- dist(matrix2021)  # Calcul des distances pour 2021

# Appliquer le MDS classique
mds_2019 <- cmdscale(dist_matrix_2019, k = 2)  # k=2 pour une représentation en 2D
mds_2021 <- cmdscale(dist_matrix_2021, k = 2)  # k=2 pour une représentation en 2D

# Exemple de noms des états
state_names <- c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", "Delaware", "District of Columbia", 
                 "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", 
                 "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", 
                 "New Hampshire", "New Jersey", "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", 
                 "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", 
                 "Washington", "West Virginia", "Wisconsin", "Wyoming", "Puerto Rico")

# Visualisation des résultats avec noms des paires d'états
plot(mds_2019, main = "MDS - Distribution des états en 2019", xlab = "Dimension 1", ylab = "Dimension 2", col = "blue", pch = 16)
points(mds_2021, col = "red", pch = 16)  # Superposer les résultats de 2021 en rouge

# Ajouter les noms des paires d'états
text(mds_2019, labels = state_names, pos = 4, col = "blue")  # Noms pour 2019
text(mds_2021, labels = state_names, pos = 4, col = "red")  # Noms pour 2021

# Légende pour le graphique
legend("topright", legend = c("2019", "2021"), col = c("blue", "red"), pch = 16)
