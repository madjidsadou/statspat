matrix2019 <- as.matrix(FLOW2019 )
state_names <- c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", "Delaware", "District of Columbia", 
                 "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", 
                 "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", 
                 "New Hampshire", "New Jersey", "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", 
                 "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", 
                 "Washington", "West Virginia", "Wisconsin", "Wyoming", "Puerto Rico")
dstate <- as.data.frame(state_names)


# Assign state names to the matrix rows and columns
rownames(matrix2019) <- state_names
colnames(matrix2019) <- state_names

# View the matrix with row and column names

rownames(matrix2019) <- state_names
head(matrix2019)

df2019 <- as.data.frame(matrix2019)

matrix2021 <- as.matrix(FLOW2021 )
print(matrix2021)

state_names <- c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", "Delaware", "District of Columbia", 
                 "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", 
                 "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", 
                 "New Hampshire", "New Jersey", "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", 
                 "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", 
                 "Washington", "West Virginia", "Wisconsin", "Wyoming", "Puerto Rico")
# Assign state names to the matrix rows and columns
rownames(matrix2021) <- state_names
colnames(matrix2021) <- state_names

# View the matrix with row and column names

rownames(matrix2021) <- state_names
head(matrix2021)

df2021 <- as.data.frame(matrix2021)


transition_matrix2021 <- sweep(df2021, 1, rowSums(df2021), "/")
transition_matrix2019 <- sweep(df2019, 1, rowSums(df2019), "/")

row_sums <- rowSums(transition_matrix2019)
print(all.equal(unname(row_sums), rep(1, nrow(transition_matrix2019))))

row_sums <- rowSums(transition_matrix2021)
print(all.equal(unname(row_sums), rep(1, nrow(transition_matrix2021))))


stationary_distribution2019 <- eigen(t(transition_matrix2019))$vectors[,1]
stationary_distribution2019 <- Re(stationary_distribution2019)  # Prendre la partie réelle
stationary_distribution2019 <- stationary_distribution2019 / sum(stationary_distribution2019)  # Normalisation



# Calcul de la distribution stationnaire pour la matrice de transition de 2021
stationary_distribution2021 <- eigen(t(transition_matrix2021))$vectors[,1]
stationary_distribution2021 <- Re(stationary_distribution2021)  # Prendre la partie réelle
stationary_distribution2021 <- stationary_distribution2021 / sum(stationary_distribution2021)  # Normalisation

sum(stationary_distribution2021)
sum(stationary_distribution2019)


# Affichage des distributions stationnaires
print("Stationary Distribution 2019:")
print(stationary_distribution2019)

print("Stationary Distribution 2021:")
print(stationary_distribution2021)

rownames(stationary_distribution2019) <- state_names
rownames(stationary_distribution2021) <- state_names

ST2019 <- as.data.frame(stationary_distribution2019)
ST2021 <- as.data.frame(stationary_distribution2021)



cat("Taille de la distribution stationnaire 2019: ", length(stationary_distribution2019), "\n")
cat("Taille de la distribution stationnaire 2021: ", length(stationary_distribution2021), "\n")

# Vérifier la longueur de state_names
cat("Taille de state_names: ", length(state_names), "\n")

# Vérification pour s'assurer que les tailles sont égales
if (length(stationary_distribution2019) == length(state_names) && length(stationary_distribution2021) == length(state_names)) {
  
  # Préparer les distributions stationnaires avec les noms des états
  stationary_df <- data.frame(State = state_names,
                              `Stationary 2019` = stationary_distribution2019,
                              `Stationary 2021` = stationary_distribution2021)
  
} else {
  print("Les tailles des distributions stationnaires ne correspondent pas au nombre d'états.")
}

#State 10: Probability = 0.0702. This suggests that in the long run, the system is most likely to be in state 10 about 7.02% of the time.
#State 52: Probability = 0.0015. This suggests that state 52 has the least chance of being visited in the long run, with only about 0.15% chance of being in that state.

rownames(ST2019) <- state_names
rownames(ST2021) <- state_names


# Function to compute the steps to reach stationary state
steps_to_stationary <- function(transition_matrix, tolerance = 1e-6, max_steps = 1000) {
  # Initialize with uniform distribution
  num_states <- nrow(transition_matrix)
  current_distribution <- rep(1 / num_states, num_states)
  steps <- 0
  
  repeat {
    # Compute the next distribution
    next_distribution <- current_distribution %*% transition_matrix
    
    # Check convergence
    if (max(abs(next_distribution - current_distribution)) < tolerance) {
      break
    }
    
    # Update and count steps
    current_distribution <- next_distribution
    steps <- steps + 1
    
    if (steps >= max_steps) {
      cat("Did not converge within the maximum number of steps.\n")
      break
    }
  }
  
  return(list(steps = steps, stationary_distribution = as.vector(current_distribution)))
}

transition_matrix2019 <- as.matrix(transition_matrix2019 )
transition_matrix2021 <- as.matrix(transition_matrix2021 )

# Calculate steps to stationary for 2019
result2019 <- steps_to_stationary(transition_matrix2019)
steps_2019 <- result2019$steps
stationary_distribution2019 <- result2019$stationary_distribution

# Calculate steps to stationary for 2021
result2021 <- steps_to_stationary(transition_matrix2021)
steps_2021 <- result2021$steps
stationary_distribution2021 <- result2021$stationary_distribution

# Print results
cat("Steps to reach stationary distribution for 2019:", steps_2019, "\n")
cat("Steps to reach stationary distribution for 2021:", steps_2021, "\n")

# Assign names to distributions
names(stationary_distribution2019) <- state_names
names(stationary_distribution2021) <- state_names

# Convert to data frames for better visualization
ST2019 <- as.data.frame(stationary_distribution2019)
ST2021 <- as.data.frame(stationary_distribution2021)

# Visualization (same as in your original code)
if (length(stationary_distribution2019) == length(state_names) && length(stationary_distribution2021) == length(state_names)) {
  
  # Prepare data frame for visualization
  stationary_df <- data.frame(State = state_names,
                              `Stationary 2019` = stationary_distribution2019,
                              `Stationary 2021` = stationary_distribution2021)
  
  # Barplot visualization
  barplot(cbind(stationary_distribution2019, stationary_distribution2021), beside = TRUE, 
          col = c("blue", "red"), 
          names.arg = state_names, 
          legend.text = c("2019", "2021"),
          main = "Comparaison des distributions stationnaires (2019 vs 2021)",
          xlab = "État", ylab = "Probabilité", las = 2, cex.names = 0.6) # Adjust text size for better readability
} else {
  print("Les tailles des distributions stationnaires ne correspondent pas au nombre d'états.")
}

flux_profiles <- t(transition_matrix2019)  # Profils migratoires sortants
dissimilarity_matrix <- as.dist(1 - cor(flux_profiles))
dissimilarity_matrix <- as.matrix(dissimilarity_matrix)

mds_result <- cmdscale(dissimilarity_matrix, k = 2)  # k = dimensions (ici 2D)

# Visualiser les résultats
plot(mds_result, type = "n", main = "MDS des flux migratoires", asp = -2)  # Prépare le graphique vide
text(mds_result, labels = state_names, col = "blue", cex = 0.05)  # Affiche les noms des états


matrice_sans_flux_internes <- df2019
diag(matrice_sans_flux_internes) <- NA

# Calculer la moyenne des lignes en excluant les flux internes
moyennes_lignes_sans_internes <- apply(matrice_sans_flux_internes, 1, function(x) mean(x, na.rm = TRUE))

# Afficher les moyennes des lignes sans les flux internes
print("Moyennes des lignes sans les flux internes :")
print(moyennes_lignes_sans_internes)
rowmean <- as.matrix(moyennes_lignes_sans_internes)

library(ggplot2)
library(sf)
library(dplyr)
library(tigris)

us_states <- states(cb = TRUE)  # "cb = TRUE" permet de récupérer les données simplifiées

# Créer une copie de la matrice pour éviter de modifier l'originale
FLOW2019_quasi_sym <- FLOW2019

# Transformer la matrice en quasi-symétrique
for (i in 1:(nrow(FLOW2019)-1)) {
  for (j in (i+1):ncol(FLOW2019)) {
    # Calculer la moyenne des flux entre les deux états
    avg_flux <- (FLOW2019[i,j] + FLOW2019[j,i]) / 2
    
    # Mettre à jour la matrice avec les valeurs moyennes
    FLOW2019_quasi_sym[i,j] <- avg_flux
    FLOW2019_quasi_sym[j,i] <- avg_flux
  }
}

# Afficher la matrice quasi-symétrique
print(FLOW2019_quasi_sym)
FLOW2019_quasi_sym <- as.matrix(FLOW2019_quasi_sym)

rownames(FLOW2019_quasi_sym) <- colnames(FLOW2019_quasi_sym) <- state_names

# Convertir la matrice en tableau (format long)
flux_data <- as.data.frame(as.table(FLOW2019_quasi_sym))

# Renommer les colonnes pour plus de clarté
colnames(flux_data) <- c("from_state", "to_state", "flow")

# Filtrer les flux internes (même état pour l'origine et la destination)
flux_data <- flux_data[flux_data$from_state != flux_data$to_state, ]

# Afficher le tableau résultant
print(head(flux_data))


library(sf)
library(dplyr)
library(purrr)


# Extraire les centroïdes des états
us_states <- us_states %>%
  mutate(
    CENTROID = st_centroid(geometry),
    INTPTLAT = st_coordinates(CENTROID)[, 2],  # Latitude
    INTPTLONG = st_coordinates(CENTROID)[, 1] # Longitude
  )

# Fusionner les données de flux avec les géométries
flux_sf <- flux_data %>%
  left_join(us_states, by = c("from_state" = "NAME")) %>%
  rename(LONG_from = INTPTLONG, LAT_from = INTPTLAT) %>%
  left_join(us_states, by = c("to_state" = "NAME")) %>%
  rename(LONG_to = INTPTLONG, LAT_to = INTPTLAT)

# Vérification
if (any(is.na(c(flux_sf$LONG_from, flux_sf$LAT_from, flux_sf$LONG_to, flux_sf$LAT_to)))) {
  stop("Des coordonnées manquent après la fusion.")
}

# Créer les géométries des flux
flux_sf <- flux_sf %>%
  mutate(
    geometry = pmap(list(LONG_from, LAT_from, LONG_to, LAT_to), function(lon1, lat1, lon2, lat2) {
      st_sfc(st_linestring(matrix(c(lon1, lon1, lat1, lat2), ncol = 2, byrow = TRUE)), crs = st_crs(us_states))
    })
  )

flux_sf <- flux_sf %>%
  rowwise() %>%
  mutate(
    geometry = st_sfc(
      st_linestring(
        matrix(c(LONG_from, LAT_from, LONG_to, LAT_to), ncol = 2, byrow = TRUE)
      ),
      crs = st_crs(us_states)
    )
  ) %>%
  ungroup()
# Convertir en objet sf
flux_sf <- st_as_sf(flux_sf, crs = st_crs(us_states))

# Visualisation
library(ggplot2)


# Define custom limits for zooming in
x_limits <- c(-160, -60)  # Longitude range (zoom into US, cut west of 60°W)
y_limits <- c(10, 70)     # Latitude range (focus on continental US)

# Plot with zoomed-in view
ggplot() +
  geom_sf(data = us_states, fill = "lightgray", color = "black") +  # Base map
  geom_sf(
    data = flux_sf, 
    aes(geometry = geometry, size = flow, color = flow), 
    alpha = 0.7
  ) +
  scale_color_gradientn(
    colors = c("blue", "yellow", "red"), 
    name = "Flow Density"
  ) +
  coord_sf(
    xlim = x_limits, 
    ylim = y_limits, 
    expand = FALSE  # Prevent extra padding around the map
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Migration Flows Between US States",
    x = "Longitude", 
    y = "Latitude"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    # Adjust legend size
    legend.key.size = unit(0.5, "cm"),  # Smaller key size
    legend.text = element_text(size = 10),  # Smaller text size for legend labels
    legend.title = element_text(size = 12)  # Smaller text size for legend title
  )

transition_matrix2019quasisym <- sweep(FLOW2019_quasi_sym, 1, rowSums(FLOW2019_quasi_sym), "/")

