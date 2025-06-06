get_repeat_annotations <- function(
    default_or_extended = "default",
    keep_non_central = TRUE) {
    rmannShared <- read_csv(confALL$aref$rmann_shared)
    if (confALL$aref$update_ref_with_tldr$response == "yes") {
        if (confALL$aref$update_ref_with_tldr$per_sample == "yes") {
            rmannSamples <- list()
            for (sample in confALL$aref$samples) {
                df <- read_csv(sprintf("aref/%s/%s_annotations/%s_rmann_nonref.csv", default_or_extended, sample, sample))
                df$sample_name <- sample
                rmannSamples[[sample]] <- df
            }
            rmannnonref <- do.call(rbind, rmannSamples) %>% tibble()
            rmann <- bind_rows(rmannShared, rmannnonref)
            if (!keep_non_central) {
                rmann <<- rmann %>% filter(refstatus != "NonCentral")
            }
        } else if (confALL$aref$update_ref_with_tldr$per_sample == "no") {
            rmann <- rmannShared
            if (!keep_non_central) {
                rmann <<- rmann %>% filter(refstatus != "NonCentral")
            }
        }
    } else {
        rmann <- rmannShared
    }

    rmann <- rmann %>%
        mutate(req_integrative = factor(req_integrative, levels = c("Old Trnc", "Old FL", "Yng Trnc", "Yng FL", "Yng Intact"))) %>%
        mutate(ltr_viral_status = factor(ltr_viral_status, levels = c("Int (Has 5LTR)", "Int (No 5'LTR)", "5'LTR (FL Int)", "3'LTR (FL Int)", "5'LTR (Trnc Int)", "3'LTR (Trnc Int)", "LTR (Solo)", "Other")))

    return(rmann)
}
