# ============================================================================
# SÉLECTION DES VARIABLES EXPLICATIVES - VERSION CORRIGÉE
# Analyse univariée avec normalisation par la durée
# ============================================================================

# Packages nécessaires
library(readxl)
library(ggplot2)
library(dplyr)

# 1. CHARGEMENT ET PRÉPARATION DES DONNÉES
# ============================================================================
cat("=== CHARGEMENT DES DONNÉES ===\n")

# Charger le fichier Excel
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

cat("Nombre de variables disponibles :", ncol(data), "\n")
cat("Nombre de périodes :", nrow(data), "\n")
cat("\nDurée des périodes :\n")
print(data[, c("periode_code", "duree_annees", "nb_eboulements")])

# 2. IDENTIFICATION DES VARIABLES À NORMALISER
# ============================================================================
cat("\n=== NORMALISATION DES VARIABLES PAR LA DURÉE ===\n")

# Variable cible : taux annuel (pas le nombre brut)
data$taux_eboulements_annuel <- data$nb_eboulements / data$duree_annees

# Variables d'identification à exclure
exclusions <- c("periode_code", "nb_eboulements", "nb_points_x", "nb_points_y", 
                "jours_disponibles", "nb_jours", "duree_annees", "taux_eboulements_annuel")

# Identifier toutes les variables numériques
vars_numeriques <- names(data)[sapply(data, is.numeric)]
vars_candidats <- setdiff(vars_numeriques, exclusions)

cat("Variables candidates pour analyse :", length(vars_candidats), "\n\n")

# Créer des versions normalisées pour les variables de comptage
# (celles qui augmentent mécaniquement avec la durée)
vars_a_normaliser <- c()
vars_deja_normalisees <- c()

for (var in vars_candidats) {
  # Vérifier si c'est un comptage ou une valeur déjà intensive
  # Les variables avec "moy", "min", "max", "std", "mediane" sont déjà normalisées
  # Les variables avec "nb_", "jours_", "energie_cumulee" doivent être normalisées
  
  if (grepl("^(nb_|jours_|energie_.*cumulee|pluie_cum|tempete)", var, ignore.case = TRUE)) {
    vars_a_normaliser <- c(vars_a_normaliser, var)
    nouveau_nom <- paste0(var, "_par_an")
    data[[nouveau_nom]] <- data[[var]] / data$duree_annees
    cat("Variable normalisée :", var, "→", nouveau_nom, "\n")
  } else {
    vars_deja_normalisees <- c(vars_deja_normalisees, var)
  }
}

cat("\nVariables déjà normalisées (moyennes, max, etc.) :", length(vars_deja_normalisees), "\n")

# Liste complète des variables à analyser
vars_normalisees <- paste0(vars_a_normaliser, "_par_an")
vars_a_analyser <- c(vars_normalisees, vars_deja_normalisees)
vars_a_analyser <- vars_a_analyser[vars_a_analyser %in% names(data)]

cat("Total de variables à analyser :", length(vars_a_analyser), "\n\n")

# 3. ANALYSE UNIVARIÉE : CORRÉLATIONS ET POISSON
# ============================================================================
cat("=== ANALYSE UNIVARIÉE (avec variables normalisées) ===\n\n")

# Variable cible
target <- "taux_eboulements_annuel"

# Initialiser le tableau de résultats
resultats <- data.frame(
  variable = character(),
  type_variable = character(),
  corr_pearson = numeric(),
  p_value_pearson = numeric(),
  corr_spearman = numeric(),
  p_value_spearman = numeric(),
  beta_poisson = numeric(),
  se_poisson = numeric(),
  p_value_poisson = numeric(),
  irr = numeric(),
  pseudo_r2 = numeric(),
  aic = numeric(),
  stringsAsFactors = FALSE
)

# Modèle nul (pour calculer le pseudo-R²)
# On modélise le TAUX directement (pas le nombre avec offset)
modele_nul <- glm(as.formula(paste(target, "~ 1")), 
                  family = gaussian(link = "identity"), 
                  data = data)
deviance_nulle <- deviance(modele_nul)

# Boucle sur chaque variable explicative
for (var in vars_a_analyser) {
  
  x <- data[[var]]
  y <- data[[target]]
  
  # Vérifier validité
  if (all(is.na(x)) || length(unique(x[!is.na(x)])) <= 1) {
    next
  }
  
  # Données complètes
  donnees_completes <- complete.cases(x, y)
  x_clean <- x[donnees_completes]
  y_clean <- y[donnees_completes]
  
  if (length(x_clean) < 3) next
  
  # Type de variable
  type_var <- ifelse(grepl("_par_an$", var), "Normalisée", "Intensive")
  
  ## CORRÉLATIONS
  test_pearson <- cor.test(y_clean, x_clean, method = "pearson")
  corr_pearson <- test_pearson$estimate
  p_pearson <- test_pearson$p.value
  
  test_spearman <- cor.test(y_clean, x_clean, method = "spearman")
  corr_spearman <- test_spearman$estimate
  p_spearman <- test_spearman$p.value
  
  ## MODÈLE LINÉAIRE (pour le taux, on utilise gaussian)
  tryCatch({
    data_temp <- data.frame(y = y_clean, x = x_clean)
    
    # Modèle linéaire gaussien (car on prédit un taux continu)
    formule <- as.formula("y ~ x")
    modele <- glm(formule, family = gaussian(link = "identity"), data = data_temp)
    
    coefs <- summary(modele)$coefficients
    beta <- coefs[2, 1]
    se <- coefs[2, 2]
    p_value <- coefs[2, 4]
    
    # Pour un modèle linéaire, IRR n'a pas de sens, on garde juste beta
    irr <- NA
    
    # Pseudo-R² 
    pseudo_r2 <- 1 - (deviance(modele) / deviance_nulle)
    
    # AIC
    aic <- AIC(modele)
    
    resultats <- rbind(resultats, data.frame(
      variable = var,
      type_variable = type_var,
      corr_pearson = corr_pearson,
      p_value_pearson = p_pearson,
      corr_spearman = corr_spearman,
      p_value_spearman = p_spearman,
      beta_poisson = beta,
      se_poisson = se,
      p_value_poisson = p_value,
      irr = irr,
      pseudo_r2 = pseudo_r2,
      aic = aic,
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    cat("Erreur avec", var, ":", e$message, "\n")
  })
}

# 4. CLASSEMENT DES VARIABLES
# ============================================================================

resultats_corr <- resultats %>%
  arrange(desc(abs(corr_pearson))) %>%
  mutate(rang_corr = row_number())

resultats_r2 <- resultats %>%
  arrange(desc(pseudo_r2)) %>%
  mutate(rang_r2 = row_number())

# 5. DIAGNOSTIC : POURQUOI LES CORRÉLATIONS SONT-elles ÉLEVÉES ?
# ============================================================================
cat("\n=== DIAGNOSTIC DES CORRÉLATIONS ÉLEVÉES ===\n\n")

# Compter les variables avec corrélation > 0.9
vars_tres_correlees <- resultats %>% filter(abs(corr_pearson) > 0.9)

cat("Nombre de variables avec |r| > 0.9 :", nrow(vars_tres_correlees), "\n\n")

if (nrow(vars_tres_correlees) > 0) {
  cat("ATTENTION : Corrélations anormalement élevées détectées !\n\n")
  cat("EXPLICATIONS POSSIBLES :\n")
  cat("1. FAIBLE NOMBRE D'OBSERVATIONS (n=5)\n")
  cat("   Avec seulement 5 périodes, les corrélations sont très instables.\n")
  cat("   Une seule période extrême peut créer une corrélation artificielle.\n\n")
  
  cat("2. MULTICOLINÉARITÉ ENTRE VARIABLES\n")
  cat("   Beaucoup de variables météo sont inter-corrélées (gel ↔ vent ↔ pression)\n\n")
  
  cat("3. BIAIS TEMPOREL (même après normalisation)\n")
  cat("   Les périodes longues peuvent avoir des caractéristiques communes.\n\n")
  
  # Analyser la distribution des valeurs
  cat("DISTRIBUTION DES TAUX D'ÉBOULEMENTS :\n")
  print(summary(data$taux_eboulements_annuel))
  cat("\n")
  
  # Vérifier s'il y a des outliers
  q1 <- quantile(data$taux_eboulements_annuel, 0.25)
  q3 <- quantile(data$taux_eboulements_annuel, 0.75)
  iqr <- q3 - q1
  outliers <- data$taux_eboulements_annuel < (q1 - 1.5*iqr) | 
              data$taux_eboulements_annuel > (q3 + 1.5*iqr)
  
  if (any(outliers)) {
    cat("PÉRIODES ATYPIQUES DÉTECTÉES :\n")
    print(data[outliers, c("periode_code", "taux_eboulements_annuel")])
    cat("\n→ Ces périodes peuvent artificiellement gonfler les corrélations.\n\n")
  }
}

# 6. AFFICHAGE DES RÉSULTATS
# ============================================================================

cat("=== TOP 15 VARIABLES PAR CORRÉLATION (après normalisation) ===\n")
print(resultats_corr[1:min(15, nrow(resultats_corr)), 
                     c("variable", "type_variable", "corr_pearson", 
                       "p_value_pearson", "pseudo_r2")])

cat("\n=== TOP 15 VARIABLES PAR PSEUDO-R² ===\n")
print(resultats_r2[1:min(15, nrow(resultats_r2)), 
                   c("variable", "type_variable", "pseudo_r2", 
                     "corr_pearson", "aic")])

# 7. RECOMMANDATIONS POUR SÉLECTION
# ============================================================================
cat("\n=== RECOMMANDATIONS POUR SÉLECTION DES VARIABLES ===\n\n")

# Score combiné
resultats_combines <- resultats %>%
  mutate(
    rang_corr = rank(-abs(corr_pearson)),
    rang_r2 = rank(-pseudo_r2),
    score_combine = rang_corr + rang_r2
  ) %>%
  arrange(score_combine)

top12 <- head(resultats_combines, 12)

cat("Top 12 variables (critère combiné) :\n\n")
for (i in 1:nrow(top12)) {
  cat(sprintf("%2d. %-40s | Type: %-10s | r=%6.3f | R²=%5.3f\n", 
              i, 
              top12$variable[i],
              top12$type_variable[i],
              top12$corr_pearson[i], 
              top12$pseudo_r2[i]))
}

cat("\n⚠️ MISE EN GARDE IMPORTANTE :\n")
cat("Avec n=5 périodes, TOUTE corrélation doit être interprétée avec prudence.\n")
cat("Privilégier des modèles simples (1-2 variables max) ou un score composite.\n\n")

# 8. VISUALISATIONS
# ============================================================================

cat("\n=== GÉNÉRATION DES GRAPHIQUES ===\n")

# Sauvegarder tous les graphiques dans des fichiers PNG
png("01_distribution_correlations.png", width = 800, height = 600)
p1 <- ggplot(resultats, aes(x = abs(corr_pearson))) +
  geom_histogram(binwidth = 0.05, fill = "#3498db", color = "white") +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 0.92, y = Inf, label = "Seuil suspect (0.9)", 
           vjust = 1.5, color = "red", fontface = "bold") +
  labs(title = "Distribution des Corrélations (valeur absolue)",
       subtitle = paste("n =", nrow(data), "périodes"),
       x = "|Corrélation de Pearson|",
       y = "Nombre de variables") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))
print(p1)
dev.off()
cat("✓ Graphique 1 sauvegardé : 01_distribution_correlations.png\n")

# Graphique 2 : Top 15 avec indication du type
png("02_top15_variables.png", width = 1000, height = 800)
top15 <- head(resultats_corr, 15)
p2 <- ggplot(top15, aes(x = reorder(variable, abs(corr_pearson)), 
                         y = corr_pearson,
                         fill = type_variable)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Normalisée" = "#2ecc71", "Intensive" = "#e67e22"),
                    name = "Type") +
  coord_flip() +
  labs(title = "Top 15 Variables - Corrélation avec Taux d'Éboulements",
       subtitle = "Variables normalisées par la durée",
       x = "",
       y = "Corrélation de Pearson") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))
print(p2)
dev.off()
cat("✓ Graphique 2 sauvegardé : 02_top15_variables.png\n")

# Graphiques 3-5 : Relation entre taux et top 3 variables
top3_vars <- head(resultats_corr$variable, 3)

for (i in 1:length(top3_vars)) {
  var <- top3_vars[i]
  filename <- paste0("0", i+2, "_relation_", gsub("_par_an", "", var), ".png")
  
  png(filename, width = 800, height = 600)
  p <- ggplot(data, aes_string(x = var, y = "taux_eboulements_annuel")) +
    geom_point(size = 4, color = "#3498db") +
    geom_smooth(method = "lm", se = TRUE, color = "#e74c3c", fill = "gray80") +
    geom_text(aes(label = periode_code), vjust = -1, size = 3.5) +
    labs(title = paste("Relation :", var, "vs Taux d'Éboulements"),
         subtitle = paste("r =", round(resultats_corr$corr_pearson[resultats_corr$variable == var], 3)),
         x = var,
         y = "Taux d'éboulements (par an)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 12))
  print(p)
  dev.off()
  
  cat(paste0("✓ Graphique ", i+2, " sauvegardé : ", filename, "\n"))
}

# Graphique 6 : Top 15 par pseudo-R²
png("06_top15_pseudoR2.png", width = 1000, height = 800)
top15_r2 <- head(resultats_r2, 15)
p6 <- ggplot(top15_r2, aes(x = reorder(variable, pseudo_r2), 
                            y = pseudo_r2,
                            fill = type_variable)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Normalisée" = "#9b59b6", "Intensive" = "#e67e22"),
                    name = "Type") +
  coord_flip() +
  labs(title = "Top 15 Variables - Pouvoir Explicatif (Pseudo-R²)",
       subtitle = "Variables avec le meilleur ajustement en régression",
       x = "",
       y = "Pseudo-R² (McFadden)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))
print(p6)
dev.off()
cat("✓ Graphique 6 sauvegardé : 06_top15_pseudoR2.png\n")

# Graphique 7 : Comparaison Corrélation vs Pseudo-R²
png("07_correlation_vs_pseudoR2.png", width = 1000, height = 800)
p7 <- ggplot(resultats, aes(x = abs(corr_pearson), y = pseudo_r2)) +
  geom_point(aes(color = type_variable, size = abs(corr_pearson)), alpha = 0.6) +
  geom_text(data = head(resultats_corr, 10), 
            aes(label = gsub("_par_an", "", variable)), 
            vjust = -0.8, size = 2.5, check_overlap = TRUE) +
  scale_color_manual(values = c("Normalisée" = "#2ecc71", "Intensive" = "#e67e22"),
                     name = "Type") +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  labs(title = "Relation entre Corrélation et Pouvoir Explicatif",
       subtitle = "Les meilleures variables sont en haut à droite",
       x = "|Corrélation de Pearson|",
       y = "Pseudo-R²") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))
print(p7)
dev.off()
cat("✓ Graphique 7 sauvegardé : 07_correlation_vs_pseudoR2.png\n")

# Graphique 8 : Top 15 par AIC (plus petit = meilleur)
png("08_top15_AIC.png", width = 1000, height = 800)
top15_aic <- resultats %>% arrange(aic) %>% head(15)
p8 <- ggplot(top15_aic, aes(x = reorder(variable, -aic), 
                             y = aic,
                             fill = type_variable)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Normalisée" = "#3498db", "Intensive" = "#e67e22"),
                    name = "Type") +
  coord_flip() +
  labs(title = "Top 15 Variables - Meilleur AIC",
       subtitle = "Plus l'AIC est faible, meilleur est le modèle",
       x = "",
       y = "AIC (Akaike Information Criterion)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))
print(p8)
dev.off()
cat("✓ Graphique 8 sauvegardé : 08_top15_AIC.png\n")

# Graphique 9 : P-values des corrélations (significativité)
png("09_significativite_correlations.png", width = 1000, height = 800)
resultats_sign <- resultats %>%
  arrange(p_value_pearson) %>%
  head(20) %>%
  mutate(significatif = ifelse(p_value_pearson < 0.05, "Significatif (p<0.05)", "Non significatif"))

p9 <- ggplot(resultats_sign, aes(x = reorder(variable, -p_value_pearson), 
                                  y = -log10(p_value_pearson),
                                  fill = significatif)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 15, y = -log10(0.05) + 0.2, 
           label = "Seuil p = 0.05", color = "red", fontface = "bold") +
  scale_fill_manual(values = c("Significatif (p<0.05)" = "#27ae60", 
                               "Non significatif" = "#95a5a6"),
                    name = "") +
  coord_flip() +
  labs(title = "Significativité Statistique des Corrélations (Top 20)",
       subtitle = "-log10(p-value) : plus c'est haut, plus c'est significatif",
       x = "",
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "top")
print(p9)
dev.off()
cat("✓ Graphique 9 sauvegardé : 09_significativite_correlations.png\n")

# Graphique 10 : Comparaison Pearson vs Spearman
png("10_pearson_vs_spearman.png", width = 1000, height = 800)
p10 <- ggplot(resultats, aes(x = corr_pearson, y = corr_spearman)) +
  geom_point(aes(color = type_variable, size = abs(corr_pearson)), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_text(data = head(resultats_corr, 8), 
            aes(label = gsub("_par_an", "", variable)), 
            vjust = -0.8, size = 2.5, check_overlap = TRUE) +
  scale_color_manual(values = c("Normalisée" = "#2ecc71", "Intensive" = "#e67e22"),
                     name = "Type") +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  labs(title = "Comparaison Corrélations : Pearson vs Spearman",
       subtitle = "Si les points sont sur la diagonale, les deux corrélations sont équivalentes",
       x = "Corrélation de Pearson (linéaire)",
       y = "Corrélation de Spearman (monotone)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))
print(p10)
dev.off()
cat("✓ Graphique 10 sauvegardé : 10_pearson_vs_spearman.png\n")

# Graphique 11 : Tableau récapitulatif des Top 12
png("11_tableau_recapitulatif_top12.png", width = 1200, height = 800)
top12_recap <- head(resultats_combines, 12) %>%
  select(variable, corr_pearson, pseudo_r2, p_value_pearson, aic) %>%
  mutate(variable = gsub("_par_an", "", variable),
         variable = substr(variable, 1, 25))  # Tronquer les noms longs

library(tidyr)
top12_long <- top12_recap %>%
  mutate(rang = row_number()) %>%
  pivot_longer(cols = c(corr_pearson, pseudo_r2), 
               names_to = "metrique", 
               values_to = "valeur") %>%
  mutate(metrique = recode(metrique, 
                          corr_pearson = "Corrélation",
                          pseudo_r2 = "Pseudo-R²"))

p11 <- ggplot(top12_long, aes(x = reorder(variable, -rang), 
                               y = valeur, 
                               fill = metrique)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Corrélation" = "#3498db", "Pseudo-R²" = "#9b59b6"),
                    name = "Métrique") +
  coord_flip() +
  labs(title = "Top 12 Variables - Comparaison Corrélation vs Pseudo-R²",
       subtitle = "Score combiné (meilleur compromis entre les deux critères)",
       x = "",
       y = "Valeur") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "top")
print(p11)
dev.off()
cat("✓ Graphique 11 sauvegardé : 11_tableau_recapitulatif_top12.png\n")

cat("\n✓ Tous les graphiques ont été sauvegardés dans le répertoire de travail.\n")
cat("  TOTAL : 11 graphiques générés\n")
cat("  - 01 : Distribution des corrélations\n")
cat("  - 02 : Top 15 par corrélation\n")
cat("  - 03-05 : Relations individuelles (top 3 variables)\n")
cat("  - 06 : Top 15 par pseudo-R²\n")
cat("  - 07 : Corrélation vs Pseudo-R²\n")
cat("  - 08 : Top 15 par AIC\n")
cat("  - 09 : Significativité statistique (p-values)\n")
cat("  - 10 : Pearson vs Spearman\n")
cat("  - 11 : Tableau récapitulatif Top 12\n\n")

# 9. EXPORT
# ============================================================================

write.csv(resultats_corr, "selection_variables_normalises_corr.csv", row.names = FALSE)
write.csv(resultats_r2, "selection_variables_normalises_r2.csv", row.names = FALSE)

cat("\n=== FICHIERS EXPORTÉS ===\n")
cat("- selection_variables_normalises_corr.csv\n")
cat("- selection_variables_normalises_r2.csv\n")