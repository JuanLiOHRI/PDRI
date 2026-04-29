# Create AUC/ROC results
ROC_summary <- function(df_temp, group_vec, test, type1, type2,
                        sex = "All", cut_all = NULL, direction = NULL) {
  
  if (is.null(direction)) {
    roc_1 <- roc(df_temp$group, df_temp$value, levels=group_vec, na.rm=TRUE, ci=TRUE) 
  } else {
    roc_1 <- roc(df_temp$group, df_temp$value, levels=group_vec, na.rm=TRUE, ci=TRUE, direction = direction) 
  }
  
  # ------- dfPlt ------- 
  # PPV, NPV
  ppv <- rep(0, length(roc_1$thresholds))
  npv <- rep(0, length(roc_1$thresholds))
  
  for (i_cut in seq_len(length(roc_1$thresholds))) {
    cut_i <-  roc_1$thresholds[i_cut]
    
    if (roc_1$direction == "<") { # controls < cases
      ppv[i_cut] <- nrow(df_temp %>% filter(value >= cut_i, group == group_vec[2]))/
        nrow(df_temp %>% filter(value >= cut_i))
      npv[i_cut] <- nrow(df_temp %>% filter(value < cut_i, group == group_vec[1]))/
        nrow(df_temp %>% filter(value < cut_i))
    } else { # controls > cases
      ppv[i_cut] <- nrow(df_temp %>% filter(value <= cut_i, group == group_vec[2]))/
        nrow(df_temp %>% filter(value <= cut_i))
      npv[i_cut] <- nrow(df_temp %>% filter(value > cut_i, group == group_vec[1]))/
        nrow(df_temp %>% filter(value > cut_i))
    }
  }
  
  # LR+, LR-
  LR_p <- roc_1$sensitivities/(1-roc_1$specificities)
  LR_n <- (1-roc_1$sensitivities)/roc_1$specificities
  
  dfPlt <- data.frame(
    thresholds = roc_1$thresholds,
    sensitivity = roc_1$sensitivities,
    specificity = roc_1$specificities,
    ppv = ppv,
    npv = npv,
    LR_p = LR_p,
    LR_n = LR_n,
    tests = rep(test, length(roc_1$sensitivities)),
    group = rep(paste0(group_vec[2], " vs ", group_vec[1]),
                length(roc_1$sensitivities)),
    sex = rep(sex, length(roc_1$sensitivities)),
    type1 = rep(type1, length(roc_1$sensitivities)),
    type2 = rep(type2, length(roc_1$sensitivities))
  )
  
  # ------- res ------- 
  temp <- (roc_1$sensitivities + roc_1$specificities)
  ind <- which(temp == max(temp, na.rm = T))[1]
  # threshold when combine both sexes
  if (is.null(cut_all)) {
    cut_all <- roc_1$thresholds[ind]
    sens_all <- roc_1$sensitivities[ind]
    spec_all <- roc_1$specificities[ind]
    ppv_all <- ppv[ind]
    npv_all <- npv[ind]
  } else {
    if (roc_1$direction == "<") { # controls < cases
      sens_all <- nrow(df_temp %>% filter(value >= cut_all, group == group_vec[2]))/
        nrow(df_temp %>% filter(group == group_vec[2]))
      spec_all <- nrow(df_temp %>% filter(value < cut_all, group == group_vec[1]))/
        nrow(df_temp %>% filter(group == group_vec[1]))
      ppv_all <- nrow(df_temp %>% filter(value >= cut_all, group == group_vec[2]))/
        nrow(df_temp %>% filter(value >= cut_all))
      npv_all <- nrow(df_temp %>% filter(value < cut_all, group == group_vec[1]))/
        nrow(df_temp %>% filter(value < cut_all))
    } else { # controls > cases
      sens_all <- nrow(df_temp %>% filter(value <= cut_all, group == group_vec[2]))/
        nrow(df_temp %>% filter(group == group_vec[2]))
      spec_all <- nrow(df_temp %>% filter(value > cut_all, group == group_vec[1]))/
        nrow(df_temp %>% filter(group == group_vec[1]))
      ppv_all <- nrow(df_temp %>% filter(value <= cut_all, group == group_vec[2]))/
        nrow(df_temp %>% filter(value <= cut_all))
      npv_all <- nrow(df_temp %>% filter(value > cut_all, group == group_vec[1]))/
        nrow(df_temp %>% filter(value > cut_all))
    }
  }
  
  sign_PD <- ifelse(roc_1$direction == "<", ">=", "<=")
  
  sens_all_0 <- sens_all
  spec_all_0 <- spec_all
  
  res <- data.frame(
    sex = sex,
    tests = test,
    group = paste0(group_vec[2], " vs ", group_vec[1]),
    auc = paste0(round(roc_1$auc,2), " (",round(roc_1$ci[1],2),"-",round(roc_1$ci[3],2),")"),
    
    # Using cut-value from combine both sexes
    cut_all = paste(sign_PD, round(cut_all,2)),
    sens_all = round(sens_all,2),
    spec_all = round(spec_all,2),
    ppv_all = round(ppv_all,2),
    npv_all = round(npv_all,2),
    LR_p_all = round(sens_all_0/(1-spec_all_0),2),
    LR_n_all = round((1-sens_all_0)/spec_all_0,2),
    
    # Sex-specific
    cut = paste(sign_PD, round(roc_1$thresholds[ind],2)),
    sens = round(roc_1$sensitivities[ind],2),
    spec = round(roc_1$specificities[ind],2),
    ppv = round(ppv[ind],2),
    npv = round(npv[ind],2),
    LR_p = round(LR_p[ind],2),
    LR_n = round(LR_n[ind],2),
    type1 = type1,
    type2 = type2
  )
  
  return(list(dfPlt = dfPlt,
              res = res,
              cut_all = cut_all))
}