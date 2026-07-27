# ==============================================================================
# Sturgeon cranial bone shape analysis - geometric morphometrics
# Species: Huso huso, Acipenser gueldenstaedtii, A. ruthenus, A. stellatus,
#          A. baerii + 2 hybrids (Bester, Sevbel)
# 11 cranial bones: 6 digitised in 2D, 5 already registered in 3D
#
# Pipeline:
#   1. Generalized Procrustes Analysis (2D bones) / read pre-aligned coords (3D)
#   2. Average left/right sides per individual, per-bone PCA
#   3. Taxonomic shape signal (procD.lm) - all species & Acipenser only
#   4. Allometry test (shape ~ log(centroid size) * species)
#   5. Bilateral (directional/fluctuating) asymmetry
#   6. Pairwise PERMANOVA between species, per bone
#   7. Hierarchical clustering (Ward.D2) on per-bone PC scores
#   8. Figures: taxonomic-signal barplot, R2 heatmaps, dendrograms,
#      per-bone PCA scatterplots
# ==============================================================================

library(geomorph)
library(vegan)
library(ggplot2)
library(ggdendro)


# ---- paths ------------------------------------------------------------------

# >>> SET YOUR OWN DATA FOLDER HERE <<<
# Folder containing all_landmarks_2D.csv and all_landmarks_3D.csv.
# Can be an absolute path, or a path relative to your working directory
# (e.g. "data" if you run the script from the repository root).
DATA_DIR <- "path/to/your/data"

f2d <- file.path(DATA_DIR, "all_landmarks_2D.csv")
f3d <- file.path(DATA_DIR, "all_landmarks_3D.csv")

# output folder for this run (tables/ + figures/); created automatically
out <- file.path(DATA_DIR, "output")
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures/pca"), showWarnings = FALSE)


# ---- params -------------------------------------------------------------------

n_perm <- 9999   # permutations for procD.lm / adonis2

pure_species <- c("Huso huso",
                  "Acipenser gueldenstaedtii",
                  "Acipenser ruthenus",
                  "Acipenser stellatus")


# ---- read data ------------------------------------------------------------

d2 <- read.csv(f2d, stringsAsFactors = FALSE)   # 2D landmarks, long format
d3 <- read.csv(f3d, stringsAsFactors = FALSE)   # 3D landmarks, already Procrustes-aligned

cat("2D:", nrow(d2), "rows;  3D:", nrow(d3), "rows\n")


# convert a wide landmark table (one row per specimen, LMi_X/_Y[/_Z] columns)
# into the p x k x n array format required by geomorph
to_array <- function(df, p, k) {
  if (k == 2) {
    cols <- paste0("LM", rep(1:p, each = 2), c("_X", "_Y"))
  } else {
    cols <- paste0("LM", rep(1:p, each = 3), c("_X", "_Y", "_Z"))
  }
  M <- as.matrix(df[, cols])
  arr <- array(t(M), dim = c(k, p, nrow(df)))
  aperm(arr, c(2, 1, 3))
}


# ---- per-bone shape data --------------------------------------------------
# GPA for the 2D-digitised bones; 3D bones are read in already aligned

bones <- list()

for (b in unique(d2$bone)) {
  s <- d2[d2$bone == b, ]
  p <- s$n_lm[1]
  gpa <- gpagen(to_array(s, p, 2), print.progress = FALSE)  # Generalized Procrustes Analysis
  bones[[b]] <- list(
    coords = gpa$coords,
    meta = data.frame(species = s$species, side = s$side,
                      fish = s$fish, Csize = gpa$Csize),
    dims = 2)
  cat(sprintf("  GPA %-22s p=%d n=%d\n", b, p, nrow(s)))
}

for (b in unique(d3$bone)) {
  s <- d3[d3$bone == b, ]
  p <- s$n_lm[1]
  bones[[b]] <- list(
    coords = to_array(s, p, 3),
    meta = data.frame(species = s$species, side = s$side,
                      fish = s$fish, Csize = s$centroid),
    dims = 3)
  cat(sprintf("  3D  %-22s p=%d n=%d\n", b, p, nrow(s)))
}

bone_list <- sort(names(bones))


# ---- average L/R sides + per-bone PCA -------------------------------------
# collapse paired bones to a single shape per individual (mean of L/R),
# then run an ordinary PCA on the averaged Procrustes coordinates

for (b in bone_list) {
  m <- bones[[b]]$meta
  C <- bones[[b]]$coords
  m$key <- paste(m$species, m$fish)   # unique individual ID
  
  uk <- unique(m$key)
  p <- dim(C)[1]; k <- dim(C)[2]
  Cavg <- array(NA, dim = c(p, k, length(uk)))
  meta_avg <- data.frame(species = NA, fish = NA, Csize = NA)[0, ]
  
  for (i in seq_along(uk)) {
    ix <- which(m$key == uk[i])
    if (length(ix) == 1) {
      Cavg[, , i] <- C[, , ix]
    } else {
      Cavg[, , i] <- apply(C[, , ix], c(1, 2), mean)   # average L/R landmarks
    }
    meta_avg <- rbind(meta_avg,
                      data.frame(species = m$species[ix[1]],
                                 fish = m$fish[ix[1]],
                                 Csize = mean(m$Csize[ix])))
  }
  
  Y <- two.d.array(Cavg)
  pc <- prcomp(Y)
  
  bones[[b]]$coords_avg <- Cavg
  bones[[b]]$meta_avg <- meta_avg
  bones[[b]]$pca <- pc$x[, 1:min(8, ncol(pc$x)), drop = FALSE]   # keep up to 8 PCs
  bones[[b]]$pca_var <- pc$sdev^2 / sum(pc$sdev^2)
}


# sample sizes per bone x species
ss <- data.frame()
for (b in bone_list) {
  tab <- table(bones[[b]]$meta_avg$species)
  ss <- rbind(ss, data.frame(bone = b,
                             species = names(tab),
                             n = as.integer(tab)))
}
ss_wide <- reshape(ss, idvar = "bone", timevar = "species", direction = "wide")
ss_wide[is.na(ss_wide)] <- 0
write.csv(ss_wide, file.path(out, "tables/sample_size.csv"), row.names = FALSE)


# PCA variance summary (how much shape variance the retained PCs capture)
pca_var <- data.frame()
for (b in bone_list) {
  v <- bones[[b]]$pca_var
  pca_var <- rbind(pca_var,
                   data.frame(bone = b,
                              dims = bones[[b]]$dims,
                              total_pcs = length(v),
                              cum_8 = sum(v[1:min(8, length(v))]),
                              cum_95_n = which(cumsum(v) >= 0.95)[1]))
}
write.csv(pca_var, file.path(out, "tables/pca_coverage.csv"), row.names = FALSE)


# ---- taxonomic shape signal -----------------------------------------------
# procD.lm (RRPP) on raw Procrustes coordinates, species as predictor

tax_all <- data.frame()
tax_acip <- data.frame()
tax_allom <- data.frame()      # Type I SS, shape ~ logCS * sp
tax_allom2 <- data.frame()     # Type II SS, sanity check

for (b in bone_list) {
  
  C <- bones[[b]]$coords_avg
  m <- bones[[b]]$meta_avg
  
  # drop species with a single specimen for this bone
  keep <- m$species %in% names(table(m$species))[table(m$species) >= 2]
  C <- C[, , keep, drop = FALSE]
  m <- m[keep, ]
  
  if (length(unique(m$species)) < 2) next
  
  # all four species
  gdf <- geomorph.data.frame(coords = C, sp = factor(m$species))
  fit <- procD.lm(coords ~ sp, data = gdf, iter = n_perm,
                  RRPP = TRUE, print.progress = FALSE)
  a <- anova(fit)$table
  tax_all <- rbind(tax_all,
                   data.frame(bone = b, dims = bones[[b]]$dims, n = dim(C)[3],
                              n_sp = length(unique(m$species)),
                              R2 = a$Rsq[1], F = a$F[1], Z = a$Z[1], p = a[1, "Pr(>F)"]))
  
  # Acipenser only (excludes Huso huso)
  ix <- m$species != "Huso huso"
  if (sum(ix) >= 4 && length(unique(m$species[ix])) >= 2) {
    gdf2 <- geomorph.data.frame(coords = C[, , ix, drop = FALSE],
                                sp = factor(m$species[ix]))
    fit2 <- procD.lm(coords ~ sp, data = gdf2, iter = n_perm,
                     RRPP = TRUE, print.progress = FALSE)
    a2 <- anova(fit2)$table
    tax_acip <- rbind(tax_acip,
                      data.frame(bone = b, dims = bones[[b]]$dims, n = sum(ix),
                                 n_sp = length(unique(m$species[ix])),
                                 R2 = a2$Rsq[1], F = a2$F[1], Z = a2$Z[1],
                                 p = a2[1, "Pr(>F)"]))
  }
  
  # allometry: shape ~ logCS * sp (do species differ once size is accounted for?)
  if (dim(C)[3] >= 6 && all(m$Csize > 0, na.rm = TRUE)) {
    
    gdf3 <- geomorph.data.frame(coords = C, sp = factor(m$species),
                                logCS = log(m$Csize))
    
    fit3 <- procD.lm(coords ~ logCS * sp, data = gdf3, iter = n_perm,
                     RRPP = TRUE, print.progress = FALSE)
    a3 <- anova(fit3)$table
    tax_allom <- rbind(tax_allom,
                       data.frame(bone = b, dims = bones[[b]]$dims, n = dim(C)[3],
                                  R2_size = a3$Rsq[1], R2_sp = a3$Rsq[2], R2_int = a3$Rsq[3],
                                  p_size = a3[1, "Pr(>F)"], p_sp = a3[2, "Pr(>F)"],
                                  p_int = a3[3, "Pr(>F)"]))
    
    # Type II SS as a sanity check on the Type I (sequential) result above
    fit4 <- procD.lm(coords ~ logCS * sp, data = gdf3, iter = n_perm,
                     RRPP = TRUE, SS.type = "II", print.progress = FALSE)
    a4 <- anova(fit4)$table
    tax_allom2 <- rbind(tax_allom2,
                        data.frame(bone = b, R2_size = a4$Rsq[1], R2_sp = a4$Rsq[2],
                                   p_sp = a4[2, "Pr(>F)"]))
  }
}

tax_all$p_FDR <- p.adjust(tax_all$p, method = "BH")
tax_acip$p_FDR <- p.adjust(tax_acip$p, method = "BH")
tax_allom$p_sp_FDR <- p.adjust(tax_allom$p_sp, method = "BH")

write.csv(tax_all,    file.path(out, "tables/tax_signal_all.csv"),     row.names = FALSE)
write.csv(tax_acip,   file.path(out, "tables/tax_signal_acip.csv"),    row.names = FALSE)
write.csv(tax_allom,  file.path(out, "tables/tax_signal_allom.csv"),   row.names = FALSE)
write.csv(tax_allom2, file.path(out, "tables/tax_signal_allom_t2.csv"),row.names = FALSE)


# ---- bilateral symmetry ---------------------------------------------------
# uses raw (un-averaged) coords, paired (L/R) bones only
# object.sym = FALSE: these are discrete paired bones, not a midline structure

bilat <- data.frame()

for (b in bone_list) {
  
  m <- bones[[b]]$meta
  C <- bones[[b]]$coords
  
  if (!any(m$side %in% c("L", "R"))) next   # skip unpaired bones
  
  m$key <- paste(m$species, m$fish)
  ct <- table(m$key, m$side)
  if (!all(c("L", "R") %in% colnames(ct))) next
  
  paired <- rownames(ct)[ct[, "L"] > 0 & ct[, "R"] > 0]
  if (length(paired) < 4) next
  
  ix <- m$key %in% paired & m$side %in% c("L", "R")
  if (sum(ix) < 8) next
  
  res <- try(bilat.symmetry(A = C[, , ix],
                            ind = factor(m$key[ix]),
                            side = factor(m$side[ix]),
                            object.sym = FALSE,
                            iter = 999, print.progress = FALSE),
             silent = TRUE)
  if (inherits(res, "try-error")) next
  
  tab <- res$shape.anova
  bilat <- rbind(bilat,
                 data.frame(bone = b, n_paired = length(paired),
                            R2_ind = tab$Rsq[grep("^ind$", rownames(tab))],
                            R2_DA  = tab$Rsq[grep("^side$", rownames(tab))],
                            R2_FA  = tab$Rsq[grep("ind:side", rownames(tab))],
                            p_DA   = tab[grep("^side$", rownames(tab)), "Pr(>F)"],
                            p_FA   = tab[grep("ind:side", rownames(tab)), "Pr(>F)"]))
}

if (nrow(bilat) > 0) {
  bilat$p_DA_FDR <- p.adjust(bilat$p_DA, method = "BH")
  write.csv(bilat, file.path(out, "tables/bilat_symmetry.csv"), row.names = FALSE)
}


# ---- pairwise PERMANOVA ---------------------------------------------------
# only test species pairs with n>=3 in each group

pw <- data.frame()

for (b in bone_list) {
  
  m <- bones[[b]]$meta_avg
  C <- bones[[b]]$coords_avg
  spp <- sort(unique(m$species))
  pairs <- combn(spp, 2)
  
  for (i in 1:ncol(pairs)) {
    
    s1 <- pairs[1, i]; s2 <- pairs[2, i]
    n1 <- sum(m$species == s1); n2 <- sum(m$species == s2)
    if (n1 < 3 || n2 < 3) next
    
    keep <- m$species %in% c(s1, s2)
    Y <- two.d.array(C[, , keep, drop = FALSE])
    sp <- m$species[keep]
    
    # euclidean distance on Procrustes coordinates == Procrustes distance
    ad <- adonis2(Y ~ sp, method = "euclidean", permutations = n_perm)
    
    # smallest p-value the permutation design could ever produce for this
    # pair (limited by how many distinct L/R relabellings exist at this n);
    # flags pairs where "n.s." may just mean "underpowered", not "no signal"
    mp <- 1 / choose(n1 + n2, min(n1, n2))
    tier <- if (mp > 0.05) "low" else if (mp > 0.014) "marginal" else "ok"
    
    pw <- rbind(pw,
                data.frame(bone = b, dims = bones[[b]]$dims,
                           sp1 = s1, sp2 = s2, n1 = n1, n2 = n2,
                           R2 = ad$R2[1], F = ad$F[1], p = ad[1, "Pr(>F)"],
                           min_p = mp, tier = tier))
  }
}

pw$p_FDR <- p.adjust(pw$p, method = "BH")
write.csv(pw, file.path(out, "tables/pairwise_PERMANOVA.csv"), row.names = FALSE)

# ordered version for the supplement
pw_supp <- pw[order(pw$bone, pw$sp1, pw$sp2), ]
write.csv(pw_supp, file.path(out, "tables/pairwise_supplement.csv"),
          row.names = FALSE)


# ---- dendrograms ------------------------------------------------------------
# core: 4 native species, n>=3 per cell
# full: all 7 taxa, n>=1 (illustrative - some hybrids have n=1 for some bones)

build_means <- function(taxa, min_n) {
  
  parts <- list()
  for (b in bone_list) {
    m <- bones[[b]]$meta_avg
    P <- bones[[b]]$pca
    
    if (!is.null(taxa)) {
      keep <- m$species %in% taxa
      m <- m[keep, ]; P <- P[keep, , drop = FALSE]
    }
    
    cnt <- table(m$species)
    keep_sp <- names(cnt)[cnt >= min_n]
    if (!length(keep_sp)) next
    
    df <- data.frame(species = m$species, P)
    df <- df[df$species %in% keep_sp, ]
    mn <- aggregate(df[, -1], list(species = df$species), mean)
    colnames(mn)[-1] <- paste(b, colnames(mn)[-1], sep = "_")
    parts[[b]] <- mn
  }
  
  M <- parts[[1]]
  for (k in 2:length(parts)) M <- merge(M, parts[[k]], by = "species", all = TRUE)
  
  X <- as.matrix(M[, -1])
  rownames(X) <- M$species
  X[, colSums(is.na(X)) == 0]
}

mat_core <- build_means(pure_species, min_n = 3)
hc_core <- hclust(dist(mat_core), method = "ward.D2")
saveRDS(hc_core, file.path(out, "tables/dendrogram_core.rds"))

mat_full <- build_means(NULL, min_n = 1)
hc_full <- hclust(dist(mat_full), method = "ward.D2")
saveRDS(hc_full, file.path(out, "tables/dendrogram_full.rds"))

cat("core dendrogram:", ncol(mat_core), "features\n")
cat("full dendrogram:", ncol(mat_full), "features\n")


# ---- plotting -----------------------------------------------------------------

sp_col <- c("Huso huso"                 = "#d62728",
            "Acipenser gueldenstaedtii" = "#1f77b4",
            "Acipenser ruthenus"        = "#2ca02c",
            "Acipenser stellatus"       = "#ff7f0e",
            "Acipenser baerii"          = "#9467bd",
            "Bester"                    = "#8c564b",
            "Sevbel"                    = "#e377c2")

sp_pch <- c("Huso huso" = 16, "Acipenser gueldenstaedtii" = 16,
            "Acipenser ruthenus" = 16, "Acipenser stellatus" = 16,
            "Acipenser baerii" = 17, "Bester" = 4, "Sevbel" = 8)


# manual 68% confidence ellipse (covariance eigen-decomposition), because
# ggplot2::stat_ellipse requires n>=4 and several groups here are smaller
ellipse_pts <- function(xy, grp, level = 0.68) {
  res <- data.frame()
  for (g in unique(grp)) {
    pts <- xy[grp == g, , drop = FALSE]
    if (nrow(pts) < 3) next
    S <- cov(pts)
    if (any(is.na(S)) || det(S) <= 0) next
    eg <- eigen(S, symmetric = TRUE)
    chi <- qchisq(level, df = 2)
    th <- seq(0, 2 * pi, length.out = 80)
    circ <- cbind(sqrt(chi * eg$values[1]) * cos(th),
                  sqrt(chi * eg$values[2]) * sin(th))
    rot <- circ %*% t(eg$vectors)
    co <- sweep(rot, 2, colMeans(pts), "+")
    res <- rbind(res, data.frame(PC1 = co[, 1], PC2 = co[, 2], sp = g))
  }
  res
}

pair_lab <- function(a, b) paste(sort(c(a, b)), collapse = " vs ")


# ---- Fig 1: taxonomic shape signal per bone --------------------------------

tax_all$level  <- "All four species"
tax_acip$level <- "Acipenser only"
common <- intersect(colnames(tax_all), colnames(tax_acip))
tx <- rbind(tax_all[, common], tax_acip[, common])
tx$level <- factor(tx$level, levels = c("All four species", "Acipenser only"))
tx$signif <- ifelse(tx$p_FDR < 0.05, "p_FDR < 0.05", "n.s.")
tx$bone <- factor(tx$bone, levels = tax_all$bone[order(tax_all$R2)])

g1 <- ggplot(tx, aes(R2, bone, fill = signif)) +
  geom_col(color = "black", size = 0.3, width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", R2)), hjust = -0.1, size = 3) +
  scale_fill_manual(values = c("p_FDR < 0.05" = "#2c7fb8", "n.s." = "grey75"),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0.01)) +
  facet_wrap(~ level) +
  labs(x = expression(R^2), y = NULL,
       title = "Taxonomic shape signal per bone") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(out, "figures/Fig1_tax_signal.png"), g1,
       width = 11, height = 5.5, dpi = 300)


# ---- Fig 2: pairwise R2 heat map (4 native species) ------------------------

pw4 <- pw[pw$sp1 %in% pure_species & pw$sp2 %in% pure_species, ]
pw4$pair <- mapply(pair_lab, pw4$sp1, pw4$sp2)

po <- aggregate(pw4$R2, list(pair = pw4$pair), mean)
po <- po$pair[order(po$x)]
bo <- aggregate(pw4$R2, list(bone = pw4$bone), mean)
bo <- bo$bone[order(-bo$x)]
pw4$pair <- factor(pw4$pair, levels = po)
pw4$bone <- factor(pw4$bone, levels = bo)

pw4$lab <- sprintf("%.2f", pw4$R2)
pw4$lab[pw4$p_FDR < 0.05] <- paste0(pw4$lab[pw4$p_FDR < 0.05], "*")
pw4$lab[pw4$tier == "low"] <- paste0(pw4$lab[pw4$tier == "low"], "\u2020")

g2 <- ggplot(pw4, aes(bone, pair, fill = R2)) +
  geom_tile(color = "white", size = 0.6) +
  geom_text(aes(label = lab, color = R2 > 0.5), size = 3.1, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20")) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b",
                      limits = c(0, 1), name = expression(R^2)) +
  labs(x = NULL, y = NULL,
       title = "Pairwise PERMANOVA R-squared",
       caption = "* p_FDR < 0.05;  \u2020 low-power tier") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        panel.grid = element_blank())

ggsave(file.path(out, "figures/Fig2_pairwise.png"), g2,
       width = 10, height = 5, dpi = 300)


# ---- Fig 3: pairwise R2 within Acipenser only -------------------------------

acip_sp <- pure_species[grep("^Acipenser", pure_species)]
pwa <- pw[pw$sp1 %in% acip_sp & pw$sp2 %in% acip_sp, ]
pwa$pair <- mapply(pair_lab, pwa$sp1, pwa$sp2)

# reuse the global FDR (calibrated on the full set of pairs, not just this subset)
pwa$lab <- sprintf("%.2f", pwa$R2)
pwa$lab[pwa$p_FDR < 0.05] <- paste0(pwa$lab[pwa$p_FDR < 0.05], "*")
pwa$lab[pwa$tier == "low"] <- paste0(pwa$lab[pwa$tier == "low"], "\u2020")

po2 <- aggregate(pwa$R2, list(pair = pwa$pair), mean)
po2 <- po2$pair[order(po2$x)]
bo2 <- aggregate(pwa$R2, list(bone = pwa$bone), mean)
bo2 <- bo2$bone[order(-bo2$x)]
pwa$pair <- factor(pwa$pair, levels = po2)
pwa$bone <- factor(pwa$bone, levels = bo2)

g3 <- ggplot(pwa, aes(bone, pair, fill = R2)) +
  geom_tile(color = "white", size = 0.6) +
  geom_text(aes(label = lab, color = R2 > 0.5), size = 3.3, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20")) +
  scale_fill_gradient(low = "#f7fbff", high = "#006d2c",
                      limits = c(0, 1), name = expression(R^2)) +
  labs(x = NULL, y = NULL,
       title = "Pairwise R-squared within Acipenser",
       caption = "* p_FDR < 0.05;  \u2020 low-power tier") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        panel.grid = element_blank())

ggsave(file.path(out, "figures/Fig3_pairwise_acip.png"), g3,
       width = 10, height = 3.5, dpi = 300)


# ---- Fig 4: dendrograms ------------------------------------------------------

dendro_plot <- function(hc, ttl, cap = NULL) {
  dd <- dendro_data(as.dendrogram(hc))
  ggplot() +
    geom_segment(data = dd$segments,
                 aes(x, y, xend = xend, yend = yend), color = "grey30") +
    geom_text(data = dd$labels,
              aes(x, y, label = label), hjust = 1, size = 4, fontface = "italic") +
    coord_flip() +
    scale_y_reverse(expand = c(0.3, 0)) +
    labs(title = ttl, caption = cap) +
    theme_void(base_size = 11)
}

ggsave(file.path(out, "figures/Fig4a_dendro_core.png"),
       dendro_plot(hc_core,
                   "Ward.D2, 4 native species (n>=3 per cell)",
                   "Topological inferences from this version"),
       width = 8, height = 4, dpi = 300)

ggsave(file.path(out, "figures/Fig4b_dendro_full.png"),
       dendro_plot(hc_full,
                   "Ward.D2, all 7 taxa",
                   "Illustrative - A. baerii / Sevbel n=1 for some bones"),
       width = 8, height = 5, dpi = 300)


# ---- per-bone PCA scatterplots ----------------------------------------------

for (b in bone_list) {
  
  P <- bones[[b]]$pca
  if (nrow(P) < 5) next
  
  df <- data.frame(PC1 = P[, 1], PC2 = P[, 2],
                   sp = bones[[b]]$meta_avg$species)
  el <- ellipse_pts(df[, c("PC1", "PC2")], df$sp)
  
  p <- ggplot(df, aes(PC1, PC2)) +
    geom_point(aes(color = sp, shape = sp), size = 3) +
    scale_color_manual(values = sp_col, name = NULL) +
    scale_fill_manual(values = sp_col, guide = "none") +
    scale_shape_manual(values = sp_pch, name = NULL) +
    labs(title = b,
         subtitle = sprintf("n=%d (%dD; 68%% ellipses, n>=3)",
                            nrow(df), bones[[b]]$dims)) +
    theme_classic(base_size = 11) +
    theme(legend.text = element_text(face = "italic"))
  
  if (nrow(el) > 0) {
    p <- p + geom_polygon(data = el, aes(PC1, PC2, fill = sp, group = sp),
                          alpha = 0.15, color = NA, inherit.aes = FALSE)
  }
  
  ggsave(file.path(out, "figures/pca", paste0(gsub(" ", "_", b), ".png")),
         p, width = 7.5, height = 5, dpi = 300)
}


# ---- session info (package/R versions used for this run, for reproducibility) --

sink(file.path(out, "tables/sessionInfo.txt"))
print(sessionInfo())
sink()
