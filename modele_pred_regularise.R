# ============================================================================
# Modèle de Poisson RÉGULARISÉ avec LOOCV
# Adaptation pour éviter l'overfitting avec peu de données
# ============================================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(corrplot)
library(glmnet)  # Pour la régularisation LASSO/Ridge

# 1. CHARGEMENT ET PRÉPARATION DES DONNÉES
# ============================================================================
cat("=== CHARGEMENT DES DONNÉES ===\n")

data <- read_excel("indicateurs_pour_R.xlsx")

# Normaliser les noms de colonnes
names(data) <- tolower(names(data))
names(data) <- gsub(" ", "_", names(data))

# Calculer la durée de chaque période
data$duree_annees <- sapply(strsplit(data$periode_code, "-"), function(x) {
  debut <- as.numeric(x[1])
  fin <- as.numeric(x[2])
  return(fin - debut + 1)
})

# Calculer le taux annuel
data$taux_annuel <- data$nb_eboulements / data$duree_annees

# Classes de risque basées sur terciles
terciles <- quantile(data$taux_annuel, probs = c(1/3, 2/3))
data$classe_observee <- cut(data$taux_annuel, 
                            breaks = c(-Inf, terciles[1], terciles[2], Inf),
                            labels = c("faible", "moyen", "fort"))

cat("Données chargées :", nrow(data), "périodes\n\n")

# 2. NORMALISATION DES VARIABLES PAR LA DURÉE
# ============================================================================
cat("=== NORMALISATION DES VARIABLES ===\n")

# Variables à normaliser (celles qui sont des comptages cumulés)
vars_a_normaliser <- c(
  "nb_seq_gel_3j", "nb_seq_seche_10j", "jours_gel", "jours_vent_fort_60",
  "jours_gel_degel", "jours_basse_pression", "jours_humide_90",
  "jours_forte_pluie", "nb_seq_humide_3j", "pluie_apres_gel",
  "nb_jours_ampli_sup_10", "energie_vent_cumulee"
)

# Créer des versions normalisées (par an)
for (var in vars_a_normaliser) {
  if (var %in% names(data)) {
    nouveau_nom <- paste0(var, "_par_an")
    data[[nouveau_nom]] <- data[[var]] / data$duree_annees
    cat("Variable normalisée :", nouveau_nom, "\n")
  }
}

cat("\n")

# 3. ANALYSE DE MULTICOLINÉARITÉ
# ============================================================================
cat("=== ANALYSE DE MULTICOLINÉARITÉ ===\n")

# Top 15 meilleures variables identifiées par l'analyse de sélection
# (combinaison score corrélation + pseudo-R²)
top15_variables <- c(
  "nb_jours_pmermin",
  "rafale_max_kmh",
  "nb_seq_depression_3j",
  "jours_vent_fort_60",
  "plus_longue_serie_tempete_consecutive",
  "t02_max",
  "marnage_std_m",
  "energie_vent_cumulee",
  "nb_seq_seche_10j",
  "pression_min_hpa",
  "nb_jours_vent_dir_ouest",
  "t02_moy",
  "ifm",
  "jours_pluie",
  "nb_combinaisons_critiques"
)

# Créer les versions normalisées pour celles qui en ont besoin
vars_normalisees <- c()
for (var in top15_variables) {
  if (grepl("^(nb_|jours_|energie_.*cumulee)", var, ignore.case = TRUE)) {
    nouveau_nom <- paste0(var, "_par_an")
    if (!nouveau_nom %in% names(data)) {
      data[[nouveau_nom]] <- data[[var]] / data$duree_annees
    }
    vars_normalisees <- c(vars_normalisees, nouveau_nom)
  } else {
    vars_normalisees <- c(vars_normalisees, var)
  }
}

vars_normalisees <- vars_normalisees[vars_normalisees %in% names(data)]
cat("Variables sélectionnées pour LOOCV :", length(vars_normalisees), "\n")
print(vars_normalisees)

# Matrice de corrélation entre variables
matrice_corr <- cor(data[, vars_normalisees], use = "complete.obs")

# Visualiser la matrice de corrélation
corrplot(matrice_corr, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, tl.cex = 0.7,
         title = "Matrice de corrélation entre variables (normalisées)",
         mar = c(0,0,2,0))

# Identifier les paires de variables très corrélées (> 0.8)
corr_fortes <- which(abs(matrice_corr) > 0.8 & abs(matrice_corr) < 1, arr.ind = TRUE)
if (nrow(corr_fortes) > 0) {
  cat("\nVariables fortement corrélées entre elles (|r| > 0.8) :\n")
  for (i in 1:nrow(corr_fortes)) {
    var1 <- rownames(matrice_corr)[corr_fortes[i, 1]]
    var2 <- colnames(matrice_corr)[corr_fortes[i, 2]]
    if (var1 < var2) {  # Éviter les doublons
      cat(sprintf("  %s <-> %s : r = %.3f\n", 
                  var1, var2, matrice_corr[corr_fortes[i, 1], corr_fortes[i, 2]]))
    }
  }
  cat("\nATTENTION : Ces variables sont redondantes, garder seulement l'une d'elles.\n")
}

# 4. SÉLECTION RÉDUITE DE VARIABLES (RÈGLE : n/3 maximum)
# ============================================================================
cat("\n=== SÉLECTION RÉDUITE DE VARIABLES ===\n")
cat("Avec 5 périodes, maximum recommandé : 1-2 variables\n\n")

# Calculer les corrélations avec nb_eboulements pour chaque variable normalisée
correlations_target <- sapply(vars_normalisees, function(v) {
  cor(data[[v]], data$nb_eboulements, use = "complete.obs")
})

# Trier par corrélation absolue
correlations_triees <- sort(abs(correlations_target), decreasing = TRUE)

cat("Top 5 variables normalisées par corrélation :\n")
for (i in 1:min(5, length(correlations_triees))) {
  var <- names(correlations_triees)[i]
  cat(sprintf("%d. %s : r = %.3f\n", i, var, correlations_target[var]))
}

# Stratégie 1 : Utiliser uniquement les 2 meilleures variables NON corrélées
cat("\n--- STRATÉGIE 1 : 2 VARIABLES NON CORRÉLÉES ---\n")

variables_selectionnees_1 <- character()
for (var in names(correlations_triees)) {
  if (length(variables_selectionnees_1) == 0) {
    variables_selectionnees_1 <- var
  } else {
    # Vérifier si la variable n'est pas trop corrélée avec celles déjà sélectionnées
    corr_avec_selectionnees <- sapply(variables_selectionnees_1, function(v) {
      abs(cor(data[[var]], data[[v]], use = "complete.obs"))
    })
    if (all(corr_avec_selectionnees < 0.7)) {
      variables_selectionnees_1 <- c(variables_selectionnees_1, var)
    }
  }
  if (length(variables_selectionnees_1) >= 2) break
}

cat("Variables sélectionnées (non corrélées entre elles) :\n")
print(variables_selectionnees_1)

# Stratégie 2 : Score composite (moyenne de plusieurs indicateurs)
cat("\n--- STRATÉGIE 2 : SCORE COMPOSITE ---\n")

# Standardiser les variables (moyenne=0, écart-type=1)
data_std <- data
for (var in vars_normalisees) {
  data_std[[var]] <- scale(data[[var]])
}

# Créer un score composite (moyenne des top 5 variables standardisées)
top5_vars <- names(correlations_triees)[1:min(5, length(correlations_triees))]
data$score_composite <- rowMeans(data_std[, top5_vars], na.rm = TRUE)

cat("Score composite créé à partir de :\n")
print(top5_vars)
cat("\nCorrélation score composite vs nb_eboulements :", 
    round(cor(data$score_composite, data$nb_eboulements), 3), "\n")

# 5. LOOCV AVEC RÉGULARISATION (LASSO)
# ============================================================================
cat("\n=== LOOCV AVEC RÉGULARISATION LASSO ===\n")

n <- nrow(data)
resultats_lasso <- data.frame(
  periode = data$periode_code,
  taux_observe = data$taux_annuel,
  classe_observee = data$classe_observee,
  taux_predit = NA,
  classe_predite = NA
)

# Préparer la matrice de prédicteurs (toutes les variables normalisées)
X <- as.matrix(data[, vars_normalisees])
y <- data$nb_eboulements
offset_log <- log(data$duree_annees)

for (i in 1:n) {
  cat(paste("\nItération", i, "- Test sur:", data$periode_code[i], "\n"))
  
  # Train/Test split
  X_train <- X[-i, , drop = FALSE]
  y_train <- y[-i]
  offset_train <- offset_log[-i]
  
  X_test <- X[i, , drop = FALSE]
  offset_test <- offset_log[i]
  
  # LASSO avec validation croisée pour choisir lambda optimal
  # alpha=1 pour LASSO (pénalité L1 qui fait de la sélection de variables)
  cv_lasso <- cv.glmnet(X_train, y_train, 
                        family = "poisson",
                        offset = offset_train,
                        alpha = 1,
                        nfolds = 3)  # 3-fold car peu de données
  
  # Prédiction avec lambda optimal
  nb_predit <- predict(cv_lasso, newx = X_test, 
                       newoffset = offset_test,
                       s = "lambda.min", 
                       type = "response")
  
  taux_predit <- nb_predit / data$duree_annees[i]
  
  classe_predite <- cut(taux_predit, 
                        breaks = c(-Inf, terciles[1], terciles[2], Inf),
                        labels = c("faible", "moyen", "fort"))
  
  resultats_lasso$taux_predit[i] <- taux_predit
  resultats_lasso$classe_predite[i] <- as.character(classe_predite)
  
  cat(paste("  Taux observé:", round(data$taux_annuel[i], 2), 
            "| Taux prédit:", round(taux_predit, 2), "\n"))
}

# 6. LOOCV AVEC MODÈLE SIMPLE (2 VARIABLES)
# ============================================================================
cat("\n=== LOOCV AVEC MODÈLE SIMPLE (2 VARIABLES) ===\n")

resultats_simple <- data.frame(
  periode = data$periode_code,
  taux_observe = data$taux_annuel,
  classe_observee = data$classe_observee,
  taux_predit = NA,
  classe_predite = NA
)

for (i in 1:n) {
  cat(paste("\nItération", i, "- Test sur:", data$periode_code[i], "\n"))
  
  train_data <- data[-i, ]
  test_data <- data[i, ]
  
  # Modèle avec seulement 2 variables non corrélées
  formule <- as.formula(paste("nb_eboulements ~", 
                              paste(variables_selectionnees_1, collapse = " + "),
                              "+ offset(log(duree_annees))"))
  
  modele <- glm(formule, data = train_data, family = poisson(link = "log"))
  
  nb_predit <- predict(modele, newdata = test_data, type = "response")
  taux_predit <- nb_predit / test_data$duree_annees
  
  classe_predite <- cut(taux_predit, 
                        breaks = c(-Inf, terciles[1], terciles[2], Inf),
                        labels = c("faible", "moyen", "fort"))
  
  resultats_simple$taux_predit[i] <- taux_predit
  resultats_simple$classe_predite[i] <- as.character(classe_predite)
  
  cat(paste("  Taux observé:", round(test_data$taux_annuel, 2), 
            "| Taux prédit:", round(taux_predit, 2), "\n"))
}

# 7. LOOCV AVEC SCORE COMPOSITE
# ============================================================================
cat("\n=== LOOCV AVEC SCORE COMPOSITE ===\n")

resultats_composite <- data.frame(
  periode = data$periode_code,
  taux_observe = data$taux_annuel,
  classe_observee = data$classe_observee,
  taux_predit = NA,
  classe_predite = NA
)

for (i in 1:n) {
  cat(paste("\nItération", i, "- Test sur:", data$periode_code[i], "\n"))
  
  train_data <- data[-i, ]
  test_data <- data[i, ]
  
  # Modèle avec score composite unique
  formule <- as.formula("nb_eboulements ~ score_composite + offset(log(duree_annees))")
  modele <- glm(formule, data = train_data, family = poisson(link = "log"))
  
  nb_predit <- predict(modele, newdata = test_data, type = "response")
  taux_predit <- nb_predit / test_data$duree_annees
  
  classe_predite <- cut(taux_predit, 
                        breaks = c(-Inf, terciles[1], terciles[2], Inf),
                        labels = c("faible", "moyen", "fort"))
  
  resultats_composite$taux_predit[i] <- taux_predit
  resultats_composite$classe_predite[i] <- as.character(classe_predite)
  
  cat(paste("  Taux observé:", round(test_data$taux_annuel, 2), 
            "| Taux prédit:", round(taux_predit, 2), "\n"))
}

# 8. COMPARAISON DES 3 APPROCHES
# ============================================================================
cat("\n=== COMPARAISON DES 3 APPROCHES ===\n\n")

evaluer_modele <- function(resultats, nom_modele) {
  mae <- mean(abs(resultats$taux_observe - resultats$taux_predit))
  erreur_relative <- mean(abs(resultats$taux_observe - resultats$taux_predit) / 
                          resultats$taux_observe) * 100
  r2 <- cor(resultats$taux_observe, resultats$taux_predit)^2
  
  matrice_conf <- table(resultats$classe_observee, resultats$classe_predite)
  accuracy <- sum(diag(matrice_conf)) / sum(matrice_conf)
  
  cat(paste("--- ", nom_modele, " ---\n"))
  cat(paste("MAE                :", round(mae, 2), "éboulements/an\n"))
  cat(paste("Erreur relative    :", round(erreur_relative, 1), "%\n"))
  cat(paste("R²                 :", round(r2, 3), "\n"))
  cat(paste("Accuracy           :", round(accuracy * 100, 1), "%\n\n"))
  
  return(list(mae = mae, erreur_rel = erreur_relative, r2 = r2, accuracy = accuracy))
}

perf_lasso <- evaluer_modele(resultats_lasso, "LASSO (régularisé)")
perf_simple <- evaluer_modele(resultats_simple, "SIMPLE (2 variables)")
perf_composite <- evaluer_modele(resultats_composite, "COMPOSITE (score unique)")

# 9. VISUALISATION COMPARATIVE
# ============================================================================

cat("\n=== GÉNÉRATION DES GRAPHIQUES ===\n")

# Combiner les résultats
tous_resultats <- rbind(
  data.frame(resultats_lasso, modele = "LASSO"),
  data.frame(resultats_simple, modele = "Simple (2 var)"),
  data.frame(resultats_composite, modele = "Composite")
)

# Graphique 1 : Comparatif des 3 approches
png("loocv_01_comparatif_3approches.png", width = 1200, height = 600)
p_comparatif <- ggplot(tous_resultats, aes(x = taux_observe, y = taux_predit, color = modele)) +
  geom_point(size = 4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  facet_wrap(~modele) +
  labs(title = "Comparaison des 3 Approches : Taux Observé vs Prédit",
       subtitle = "La ligne diagonale = prédiction parfaite",
       x = "Taux observé (éboulements/an)",
       y = "Taux prédit (éboulements/an)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 12))
print(p_comparatif)
dev.off()
cat("✓ Graphique 1 sauvegardé : loocv_01_comparatif_3approches.png\n")

# Graphique 2 : Barplot comparatif par période
png("loocv_02_comparaison_par_periode.png", width = 1200, height = 700)
tous_resultats_long <- tous_resultats %>%
  select(periode, taux_observe, taux_predit, modele) %>%
  pivot_longer(cols = c(taux_observe, taux_predit), 
               names_to = "type", values_to = "taux") %>%
  mutate(type = recode(type, 
                       taux_observe = "Observé", 
                       taux_predit = "Prédit"))

p_barplot <- ggplot(tous_resultats_long, aes(x = periode, y = taux, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  facet_wrap(~modele, ncol = 1) +
  scale_fill_manual(values = c("Observé" = "#3498db", "Prédit" = "#e74c3c")) +
  labs(title = "Comparaison Taux Observé vs Prédit par Période et Modèle",
       x = "Période",
       y = "Taux annuel d'éboulements",
       fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "top",
        strip.text = element_text(face = "bold", size = 11))
print(p_barplot)
dev.off()
cat("✓ Graphique 2 sauvegardé : loocv_02_comparaison_par_periode.png\n")

# Graphique 3 : Erreurs par modèle
png("loocv_03_erreurs_par_modele.png", width = 1000, height = 700)
tous_resultats_erreurs <- tous_resultats %>%
  mutate(erreur = taux_predit - taux_observe,
         erreur_abs = abs(erreur))

p_erreurs <- ggplot(tous_resultats_erreurs, aes(x = modele, y = erreur_abs, fill = modele)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  scale_fill_manual(values = c("LASSO" = "#9b59b6", 
                               "Simple (2 var)" = "#e67e22", 
                               "Composite" = "#2ecc71")) +
  labs(title = "Distribution des Erreurs Absolues par Modèle",
       subtitle = "Plus les valeurs sont basses, meilleur est le modèle",
       x = "Modèle",
       y = "Erreur absolue (éboulements/an)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "none")
print(p_erreurs)
dev.off()
cat("✓ Graphique 3 sauvegardé : loocv_03_erreurs_par_modele.png\n")

# Graphique 4 : Métriques de performance
png("loocv_04_metriques_performance.png", width = 1000, height = 700)
performances_long <- performances %>%
  pivot_longer(cols = c(MAE, Erreur_pct, R2, Accuracy), 
               names_to = "metrique", 
               values_to = "valeur") %>%
  mutate(metrique = factor(metrique, 
                          levels = c("MAE", "Erreur_pct", "R2", "Accuracy"),
                          labels = c("MAE\n(éboul./an)", "Erreur\nrelative (%)", 
                                    "R²", "Accuracy")))

p_metriques <- ggplot(performances_long, aes(x = Modele, y = valeur, fill = Modele)) +
  geom_col(width = 0.7) +
  facet_wrap(~metrique, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("LASSO" = "#9b59b6", 
                               "Simple (2 var)" = "#e67e22", 
                               "Composite" = "#2ecc71")) +
  labs(title = "Comparaison des Métriques de Performance",
       subtitle = "Pour MAE et Erreur % : plus bas = meilleur | Pour R² et Accuracy : plus haut = meilleur",
       x = "",
       y = "Valeur") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 11))
print(p_metriques)
dev.off()
cat("✓ Graphique 4 sauvegardé : loocv_04_metriques_performance.png\n")

# Graphique 5 : Classification observée vs prédite (meilleur modèle)
png("loocv_05_classification_meilleur_modele.png", width = 1000, height = 700)
meilleur_modele <- performances$Modele[which.min(performances$MAE)]
resultats_meilleur <- switch(meilleur_modele,
                             "LASSO" = resultats_lasso,
                             "Simple (2 var)" = resultats_simple,
                             "Composite" = resultats_composite)

resultats_meilleur$prediction_correcte <- 
  ifelse(resultats_meilleur$classe_observee == resultats_meilleur$classe_predite, 
         "Correct", "Incorrect")

p_classif <- ggplot(resultats_meilleur, aes(x = periode, y = 1)) +
  geom_segment(aes(x = periode, xend = periode, y = 0.7, yend = 1.3), 
               color = "gray80", linewidth = 2) +
  geom_point(aes(y = 1.3, fill = classe_observee), 
             shape = 21, size = 10, color = "black", stroke = 1.5) +
  geom_text(aes(y = 1.5, label = "OBSERVÉ"), size = 3, fontface = "bold") +
  geom_point(aes(y = 0.7, fill = classe_predite), 
             shape = 24, size = 10, color = "black", stroke = 1.5) +
  geom_text(aes(y = 0.5, label = "PRÉDIT"), size = 3, fontface = "bold") +
  geom_text(aes(y = 1, label = prediction_correcte, 
                color = prediction_correcte), 
            size = 4.5, fontface = "bold") +
  scale_fill_manual(name = "Classe de risque",
                    values = c("faible" = "#2ecc71", 
                              "moyen" = "#f39c12", 
                              "fort" = "#e74c3c")) +
  scale_color_manual(values = c("Correct" = "#27ae60", 
                                "Incorrect" = "#c0392b"),
                     guide = "none") +
  labs(title = paste("Classification par Période - Modèle :", meilleur_modele),
       subtitle = paste("Observé (●) vs Prédit (▲) | Accuracy :", 
                       round(performances$Accuracy[performances$Modele == meilleur_modele] * 100, 1), "%"),
       x = "Période",
       y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        panel.grid = element_blank())
print(p_classif)
dev.off()
cat("✓ Graphique 5 sauvegardé : loocv_05_classification_meilleur_modele.png\n")

# Graphique 6 : Classification pour les 3 modèles (style triangles/ronds)
png("loocv_06_classification_3modeles.png", width = 1400, height = 1000)

# Préparer les données pour les 3 modèles
plot_data <- rbind(
  resultats_lasso %>% mutate(modele = "LASSO"),
  resultats_simple %>% mutate(modele = "Simple (2 var)"),
  resultats_composite %>% mutate(modele = "Composite")
) %>%
  mutate(prediction_correcte = ifelse(classe_observee == classe_predite, "✓ Correct", "✗ Incorrect"))

p_classif_all <- ggplot(plot_data, aes(x = periode, y = modele)) +
  # Ligne de connexion
  geom_segment(aes(x = periode, xend = periode, y = as.numeric(factor(modele)) - 0.15, 
                   yend = as.numeric(factor(modele)) + 0.15), 
               color = "gray80", linewidth = 1.5) +
  # Classe observée (ronds en haut)
  geom_point(aes(y = as.numeric(factor(modele)) + 0.15, fill = classe_observee), 
             shape = 21, size = 8, color = "black", stroke = 1.5) +
  # Classe prédite (triangles en bas)
  geom_point(aes(y = as.numeric(factor(modele)) - 0.15, fill = classe_predite), 
             shape = 24, size = 8, color = "black", stroke = 1.5) +
  # Indicateur correct/incorrect
  geom_text(aes(y = as.numeric(factor(modele)), label = prediction_correcte, 
                color = prediction_correcte), 
            size = 3, fontface = "bold") +
  scale_fill_manual(name = "Classe de risque",
                    values = c("faible" = "#2ecc71", 
                              "moyen" = "#f39c12", 
                              "fort" = "#e74c3c")) +
  scale_color_manual(values = c("✓ Correct" = "#27ae60", 
                                "✗ Incorrect" = "#c0392b"),
                     guide = "none") +
  scale_y_continuous(breaks = 1:3, labels = c("LASSO", "Simple (2 var)", "Composite")) +
  labs(title = "Classification par Période pour les 3 Modèles",
       subtitle = "Rond (●) = Classe Observée | Triangle (▲) = Classe Prédite",
       x = "Période",
       y = "Modèle") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
        axis.text.y = element_text(size = 11, face = "bold"),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 12, color = "gray40"),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5))
print(p_classif_all)
dev.off()
cat("✓ Graphique 6 sauvegardé : loocv_06_classification_3modeles.png\n")

# Graphique 7 : Matrice de confusion pour chaque modèle
png("loocv_07_matrices_confusion.png", width = 1200, height = 400)

# Créer les matrices de confusion
matrices <- list()
for (mod in c("LASSO", "Simple (2 var)", "Composite")) {
  data_mod <- switch(mod,
                     "LASSO" = resultats_lasso,
                     "Simple (2 var)" = resultats_simple,
                     "Composite" = resultats_composite)
  
  mat <- as.data.frame(table(Observé = data_mod$classe_observee, 
                            Prédit = data_mod$classe_predite))
  mat$modele <- mod
  matrices[[mod]] <- mat
}

matrices_df <- do.call(rbind, matrices)

p_matrices <- ggplot(matrices_df, aes(x = Prédit, y = Observé, fill = Freq)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = Freq), size = 8, fontface = "bold") +
  facet_wrap(~modele, ncol = 3) +
  scale_fill_gradient(low = "#ecf0f1", high = "#3498db", name = "Nombre") +
  labs(title = "Matrices de Confusion par Modèle",
       subtitle = "La diagonale = prédictions correctes",
       x = "Classe Prédite",
       y = "Classe Observée") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        axis.text = element_text(size = 11, face = "bold"),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "right")
print(p_matrices)
dev.off()
cat("✓ Graphique 7 sauvegardé : loocv_07_matrices_confusion.png\n")

# Graphique 8 : Taux de bonne classification par période
png("loocv_08_accuracy_par_periode.png", width = 1000, height = 700)

accuracy_periode <- plot_data %>%
  group_by(periode) %>%
  summarise(
    nb_correct = sum(classe_observee == classe_predite),
    nb_total = n(),
    taux_correct = nb_correct / nb_total * 100
  )

p_accuracy_periode <- ggplot(accuracy_periode, aes(x = reorder(periode, taux_correct), 
                                                     y = taux_correct)) +
  geom_col(aes(fill = taux_correct), width = 0.7) +
  geom_text(aes(label = paste0(round(taux_correct, 0), "%")), 
            hjust = -0.2, size = 5, fontface = "bold") +
  scale_fill_gradient(low = "#e74c3c", high = "#2ecc71", guide = "none") +
  coord_flip() +
  ylim(0, 110) +
  labs(title = "Taux de Bonne Classification par Période",
       subtitle = "Sur les 3 modèles testés (3 prédictions par période)",
       x = "Période",
       y = "Taux de classification correcte (%)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        axis.text = element_text(size = 11, face = "bold"))
print(p_accuracy_periode)
dev.off()
cat("✓ Graphique 8 sauvegardé : loocv_08_accuracy_par_periode.png\n")

cat("\n✓ Tous les graphiques LOOCV ont été sauvegardés !\n")
cat("  TOTAL : 8 graphiques générés\n")
cat("  - loocv_01 : Comparatif 3 approches (scatter)\n")
cat("  - loocv_02 : Comparaison par période (barplot)\n")
cat("  - loocv_03 : Erreurs par modèle (boxplot)\n")
cat("  - loocv_04 : Métriques de performance\n")
cat("  - loocv_05 : Classification meilleur modèle\n")
cat("  - loocv_06 : Classification 3 modèles (triangles/ronds)\n")
cat("  - loocv_07 : Matrices de confusion\n")
cat("  - loocv_08 : Accuracy par période\n\n")

# 10. RECOMMANDATION FINALE
# ============================================================================
cat("\n=== RECOMMANDATION FINALE ===\n\n")

performances <- data.frame(
  Modele = c("LASSO", "Simple (2 var)", "Composite"),
  MAE = c(perf_lasso$mae, perf_simple$mae, perf_composite$mae),
  Erreur_pct = c(perf_lasso$erreur_rel, perf_simple$erreur_rel, perf_composite$erreur_rel),
  R2 = c(perf_lasso$r2, perf_simple$r2, perf_composite$r2),
  Accuracy = c(perf_lasso$accuracy, perf_simple$accuracy, perf_composite$accuracy)
)

print(performances)

meilleur <- which.min(performances$MAE)
cat(paste("\nMEILLEURE APPROCHE :", performances$Modele[meilleur], "\n"))
cat("\nAvec seulement 5 périodes, privilégier le modèle le plus simple (2 variables ou composite)\n")
cat("pour éviter l'overfitting, même si les métriques semblent moins bonnes.\n")

# Exporter les résultats
write.csv(resultats_simple, "resultats_loocv_simple.csv", row.names = FALSE)
write.csv(performances, "comparaison_modeles.csv", row.names = FALSE)

cat("\nFichiers exportés : resultats_loocv_simple.csv, comparaison_modeles.csv\n")