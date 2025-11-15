# 0)  Set working directory 
#     (Assumes this script lives in a /code/ subfolder inside the project.)
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::hasFun("getActiveDocumentContext")) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  message("Working directory set to script folder: ", getwd())
}
source("diagnostic_fcns2.R")
source("glmm_stability.R")
source("dfbeta.r")

library(tidyverse)   # readr, dplyr, ggplot2, etc.
library(lme4)        # mixed‑effects models
library(lmerTest)    # p‑values for lme4 models
library(emmeans)     # estimated marginal means
library(broom.mixed) # tidy model output

#PREP DATA--------------------------------------------------
    
      
      root    <- normalizePath("..")          # parent of /code/ = project root
      csv_dir <- file.path(root, "glm data")  # location of CSVs and output PNGs
      
      # 1)  User setting 
      distance_col <- "distance_01"            # change if you prefer another column
      
      # 3)  Helper to load and prep one CSV 
      prep_df <- function(fname) {
        df <- read_csv(file.path(csv_dir, fname), show_col_types = FALSE)
        if (!distance_col %in% names(df)) {
          stop(sprintf("Column '%s' missing in %s", distance_col, fname))
        }
        df %>%
          mutate(stage  = factor(stage,  levels = c("before", "after")),
                 metric = factor(metric),
                 pair   = factor(pair))
      }
      
      # 4)  Load data 
      cat("\nLoading data …\n")
      
      #spec_whole      <- prep_df("glm_data_spec_whole.csv")
      spec_partner    <- prep_df("glm_data_spec_paired.csv")
      spec_nonpartner <- prep_df("glm_data_spec_non_paired.csv")
      #seq_whole       <- prep_df("glm_data_seq_whole.csv")
      seq_partner     <- prep_df("glm_data_seq_partner.csv")
      seq_nonpartner  <- prep_df("glm_data_seq_non_partner.csv")
      
      

#seq_partner----------------------------------------------------------------------

mod_seq_partner <- lmer(distance_01~ stage  + (1|metric) + (1|pair) +(1|session_focal)+ (1|session_partner), 
                        data = seq_partner,    REML = FALSE)

#check what random slopes to include
  rs1 <- fe.re.tab(fe.model = "distance_01~ stage", re='(1|metric) + (1|pair) +(1|session_focal)+ (1|session_partner)', 
                   data=seq_partner)
  rs1$summary
    # stage within metric - yes
    # stage within pair - yes
    # stage within session focal - yes
    # stage within session partner - yes

#dummy code data to center it
  t.data <- rs1$data
  t.data$stage.after <- t.data$stage.after-mean(t.data$stage.after) 
  


#include random slopes
  mod_seq_partner_rs <- lmer(distance_01~ stage  + (1+stage.after|metric) + (1+stage.after|pair) +(1+stage.after|session_focal)+ (1+stage.after|session_partner), data = t.data, REML = FALSE)
  summary(mod_seq_partner_rs)

#basic diagnostics
  diagnostics.plot(mod_seq_partner) #looks ok

#BLUPs
 ranef.diagn.plot(mod_seq_partner) #good enough


#full-null comparison
 mod_seq_partner_rs_null <- lmer(distance_01~ (1+stage.after|metric) + (1+stage.after|pair) +(1+stage.after|session_focal)+ (1+stage.after|session_partner), data = t.data, REML = FALSE)
 as.data.frame(anova(mod_seq_partner_rs_null, mod_seq_partner_rs,  test="Chisq")) 
  #not-sig

#seq_nonpartner----------------------------------------------------------------------
 
 mod_seq_nonpartner <- lmer(distance_01~ stage  + (1|metric) + (1|pair) +(1|session_focal)+ (1|session_partner), 
                         data = seq_nonpartner,    REML = FALSE)
 
#check which random slopes to include 
  rs1 <- fe.re.tab(fe.model = "distance_01~ stage", re='(1|metric) + (1|pair) +(1|session_focal)+ (1|session_partner)', 
                  data=seq_nonpartner)
   rs1$summary
     # stage within metric - yes, include
     # stage within pair - yes, include
     # stage within session focal - yes, include
     # stage within session partner - yes, include
 
 #dummy code data to center
   t.data <- rs1$data
   t.data$stage.after <- t.data$stage.after-mean(t.data$stage.after) 
 
 #include random slopes
   mod_seq_nonpartner_rs <- lmer(distance_01~ stage  + (1+stage.after|metric) + (1+stage.after|pair) +(1+stage.after|session_focal)+ (1|session_partner), data = t.data, REML = FALSE)
   summary(mod_seq_nonpartner_rs)
        #converges without the random slope within session partner - not sure why? please look into it
   
 #basic diagnostics
  diagnostics.plot(mod_seq_nonpartner_rs) 
      #looks ok
 
 #BLUPs
  ranef.diagn.plot(mod_seq_nonpartner_rs) 
      #something wrong  - session partner
 
 #full-null comparison
   mod_seq_nonpartner_rs_null <- lmer(distance_01~ (1+stage.after|metric) + (1+stage.after|pair) +(1+stage.after|session_focal)+ (1|session_partner), data = t.data, REML = FALSE)
   as.data.frame(anova(mod_seq_nonpartner_rs_null, mod_seq_nonpartner_rs,  test="Chisq")) 
        #full-null is significant - good
   
     
     #spec_partner----------------------------------------------------------------------
     
     mod_spec_partner <- lmer(distance_01~ stage  + (1|metric) + (1|pair) +(1|session_focal)+ (1|session_partner), 
                            data = spec_partner,    REML = FALSE)
     #run model
     rs1 <- fe.re.tab(fe.model = "distance_01~ stage", re='(1|metric) + (1|pair) +(1|session_focal)+ (1|session_partner)', 
                      data=spec_partner)
     rs1$summary
     # stage within metric - yes, include
     # stage within pair - yes, include
     # stage within session focal - yes, include
     # stage within session partner - yes, include
     
     #dummy code data to center
     t.data <- rs1$data
     t.data$stage.after <- t.data$stage.after-mean(t.data$stage.after) 
     
     #include random slopes
     mod_spec_partner_rs <- lmer(distance_01~ stage  + (1+stage.after|metric) + (1+stage.after|pair) +(1+stage.after|session_focal)+ (1+stage.after|session_partner), data = t.data, REML = FALSE)
     summary(mod_spec_partner_rs)
     
     #basic diagnostics
     diagnostics.plot(mod_spec_partner_rs) 
     #looks ok
     
     #BLUPs
     ranef.diagn.plot(mod_spec_partner_rs) 
     
     #full-null comparison
     mod_spec_partner_rs_null <- lmer(distance_01~ 1 + (1+stage.after|metric) + (1+stage.after|pair) +(1+stage.after|session_focal)+ (1+stage.after|session_partner), data = t.data, REML = FALSE)
     as.data.frame(anova(mod_spec_partner_rs_null, mod_spec_partner_rs,  test="Chisq")) 
        #not significant
     
     
     #spec_nonpartner----------------------------------------------------------------------
     
     mod_spec_nonpartner <- lmer(distance_01~ stage  + (1|metric) + (1|pair) +(1|session_focal)+ (1|session_partner), 
                              data = spec_nonpartner,    REML = FALSE)
     #run model
     rs1 <- fe.re.tab(fe.model = "distance_01~ stage", re='(1|metric) + (1|pair) +(1|session_focal)+ (1|session_partner)', 
                      data=spec_nonpartner)
     rs1$summary
     # stage within metric - yes, include
     # stage within pair - yes, include
     # stage within session focal - yes, include
     # stage within session partner - yes, include
     
     #dummy code data to center
     t.data <- rs1$data
     t.data$stage.after <- t.data$stage.after-mean(t.data$stage.after) 
     
     #include random slopes
     mod_spec_nonpartner_rs <- lmer(distance_01~ stage  + (1+stage.after|metric) + (1+stage.after|pair) +(1+stage.after|session_focal)+ (1|session_partner), data = t.data, REML = FALSE)
     summary(mod_spec_nonpartner_rs)
          #similar problem to earlier - not converging if i include all random slopes. But removing it from within 'session_partner' solves the issue
     
     #basic diagnostics
     diagnostics.plot(mod_spec_nonpartner_rs) 
     #looks ok
     
     #BLUPs
     ranef.diagn.plot(mod_spec_nonpartner_rs) 
     
     #full-null comparison
     mod_spec_nonpartner_rs_null <- lmer(distance_01~ (1+stage.after|metric) + (1+stage.after|pair) +(1+stage.after|session_focal)+ (1|session_partner), data = t.data, REML = FALSE)
     as.data.frame(anova(mod_spec_nonpartner_rs_null, mod_spec_nonpartner_rs,  test="Chisq")) 
     #nonsignificant

##########
     # Helper pulls the "stage = after" fixed effect + CI and tags pairing/modality
     get_stage_coef2 <- function(model, pairing, modality) {
       broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) |>
         dplyr::filter(term %in% c("stageafter", "stage.after")) |>
         dplyr::transmute(
           pairing  = pairing,                   # "Partner" / "Non-partner"
           modality = modality,                  # "Sequence structure" / "Call structure"
           estimate, conf.low, conf.high
         )
     }
     
     coef_stage <- dplyr::bind_rows(
       get_stage_coef2(mod_seq_nonpartner_rs,  "Non-partner", "Sequence structure"),
       get_stage_coef2(mod_spec_nonpartner_rs, "Non-partner", "Call structure"),
       get_stage_coef2(mod_seq_partner_rs,     "Partner",     "Sequence structure"),
       get_stage_coef2(mod_spec_partner_rs,    "Partner",     "Call structure")
     ) |>
       dplyr::mutate(
         ylab = factor(
           paste0(pairing, "\n", modality),
           levels = c(
             "Non-partner\nSequence structure",
             "Non-partner\nCall structure",
             "Partner\nSequence structure",
             "Partner\nCall structure"
           )
         ),
         pairing  = factor(pairing,  levels = c("Partner","Non-partner")),
         modality = factor(modality, levels = c("Sequence structure","Call structure"))
       )
     
     # --- New colours (viridis family)
     PARTNER_COL    <- "#2c728e"  # teal
     NONPARTNER_COL <- "#b8de29"  # lime
     
     coef_plot <- ggplot(coef_stage, aes(x = estimate, y = ylab)) +
       geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
       geom_pointrange(
         aes(xmin = conf.low, xmax = conf.high,
             colour = pairing, shape = modality),
         size = 0.8
       ) +
       scale_color_manual(values = c("Partner" = PARTNER_COL,
                                     "Non-partner" = NONPARTNER_COL)) +
       scale_shape_manual(values = c("Sequence structure" = 16,   # circle
                                     "Call structure"     = 17))  + # triangle
       labs(
         x = "Coefficient for stage",
         y = "Fitted models",
         colour = NULL,
         shape  = NULL
       ) +
       theme_classic(base_size = 14) +
       theme(
         plot.title = element_blank(),
         axis.text.y = element_text(hjust = 0.5)   # ← center-aligned y labels
       )
     
     print(coef_plot)
     
     ggsave(file.path(csv_dir, "stage_coef_partner_nonpartner.png"),
            coef_plot, dpi = 300, width = 8, height = 4.5)
     
     ggsave(file.path(csv_dir, "stage_coef_partner_nonpartner.svg"),
            coef_plot, width = 8, height = 4.5)