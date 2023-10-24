library(tidyverse)

sumpp <- function(pp) rowSums(pp * (col(pp) - 1))

first_plot_data <- function(x,
    labels = c("Younger", "Older")) {
    ndif <- sum(x$flag)
    if (ndif == 0) 
        stop(paste(deparse(substitute(x))), " contains no items flagged for DIF")
    if (ndif == x$ni) 
        stop("all items in ", paste(deparse(substitute(x))), 
            " have been flagged for DIF")
    if (x$ng != length(labels)) 
        labels <- paste("Group", 1:x$ng)

    maxcat <- ncol(x$ipar.sparse)

    theta <- seq(x$options$minTheta, x$options$maxTheta, x$options$inc)
    difitems <- (1:x$ni)[x$flag]
    difselections <- x$selection[x$flag]
    itemnames <- row.names(x$ipar.sparse)
    gpar <- array(NA, c(ndif, maxcat, x$ng))
    cpar <- as.matrix(x$ipar.sparse[1:(x$ni - ndif), ])
    pp <- array(NA, c(length(theta), ndif, maxcat, x$ng))
    gtheta <- split(x$calib.sparse$theta, x$group)
    gdensity <- matrix(0, length(theta), x$ng)
    for (i in 1:x$ng) {
        gdensity[, i] <- density(unlist(gtheta[names(table(x$group))[i]]), 
            n = length(theta), from = x$options$minTheta, to = x$options$maxTheta, 
            bw = 0.25)$y
    }
    new_names <- paste0("V", 1:length(labels))
    names(new_names) <- labels
    dat1 <- as_tibble(gdensity) |>
        rename(all_of(new_names)) |>
        pivot_longer(everything())
    dat1
}

item_plot_data <- function(x, labels = c("Younger", "Older")) {

    theta <- seq(x$options$minTheta, x$options$maxTheta, x$options$inc)
    difitems <- (1:x$ni)[x$flag]
    difselections <- x$selection[x$flag]
    itemnames <- row.names(x$ipar.sparse)
    gpar <- array(NA, c(ndif, maxcat, x$ng))
    cpar <- as.matrix(x$ipar.sparse[1:(x$ni - ndif), ])
    pp <- array(NA, c(length(theta), ndif, maxcat, x$ng))
    gtheta <- split(x$calib.sparse$theta, x$group)
    gdensity <- matrix(0, length(theta), x$ng)
    for (i in 1:x$ng) {
        gdensity[, i] <- density(unlist(gtheta[names(table(x$group))[i]]), 
            n = length(theta), from = x$options$minTheta, to = x$options$maxTheta, 
            bw = 0.25)$y
    }

    ivalues <- NULL
    for (i in 1:length(difitems)) {
        ncat <- x$ncat[difitems[i]]
        table_info <- tibble(
            xlabs = c("theta", "theta", "theta", "theta"),
            ylabs = c("Item Score", "Item Score", "probability", "Size"),
            titles = c(
            paste0("Item True Score Functions - Item ", difselections[i]),
            "Differences in Item True Score Functions",
            "Item Response Functions",
            "Impact (Weighted by Density)"))
        # p1 is theta then columns with label names
        gvalues <- tibble(theta = theta, y = seq(0, ncat - 1, along.with = theta))

        for (g in 1:x$ng) {
            gpar[i, , g] <- unlist(x$ipar.sparse[which(itemnames == 
                paste0("I", difselections[i], ".", g)), ])
            if (x$options$model == "GPCM") 
                pp[, i, 1:ncat, g] <- probgpcm(theta, gpar[i, 
                  1, g], gpar[i, 2:ncat, g])
            else pp[, i, 1:ncat, g] <- probgrm(theta, gpar[i, 
                1, g], gpar[i, 2:ncat, g])
            labi <- labels[g]
            gvalues <- mutate(gvalues, "{labi}_p1" := sumpp(pp[, i, 1:ncat, g]))
        }

        chi12 <- paste(x$stats[difitems[i], "df12"], ")=", x$stats[difitems[i], 
            "chi12"], sep = "")
        pseudo12 <- x$stats[difitems[i], paste("pseudo12.", x$options$pseudo.R2, 
            sep = "")]
        beta12 <- round(x$stats[difitems[i], "beta12"], 4)
        chi13 <- paste(x$stats[difitems[i], "df13"], ")=", x$stats[difitems[i], 
            "chi13"], sep = "")
        pseudo13 <- x$stats[difitems[i], paste("pseudo13.", x$options$pseudo.R2, 
            sep = "")]
        chi23 <- paste(x$stats[difitems[i], "df23"], ")=", x$stats[difitems[i], 
            "chi23"], sep = "")
        pseudo23 <- x$stats[difitems[i], paste("pseudo23.", x$options$pseudo.R2, 
            sep = "")]

        text_tibble <- tibble(
            x = min(gvalues$theta),
            y = (ncat - 1) * c(1, .9, .8),
            label = c(
                substitute(paste("Pr(", chi[12]^2, 
                    ",", chi12, ",", R[12]^2, "=", pseudo12, ",", Delta, 
                    "(", beta[1], ")=", beta12, sep = "")),
                substitute(paste("Pr(", 
                    chi[13]^2, ",", chi13, ",", R[13]^2, "=", pseudo13, 
                    sep = "")),
                substitute(paste("Pr(", 
                    chi[23]^2, ",", chi23, ",", R[23]^2, "=", pseudo23, 
                    sep = ""))
            )
        )
        
        for (g in 2:x$ng) {
            labi <- labels[g]
            gvalues <- mutate(gvalues, "{labi}_p2" :=
                abs(sumpp(pp[, i, 1:ncat, 1]) - sumpp(pp[, i, 1:ncat, g])))
        }
        for (g in 1:x$ng) {
            for (k in 1:ncat) {
                labi <- labels[g]
                gvalues <- mutate(gvalues, "{labi}_p3_{k}" := pp[, i, k, g])
            }
        }
 
        g_values <- NULL
        gk_values <- NULL
        for (g in 1:x$ng) {
            g_values[[g]] <-paste(round(gpar[i, , g][!is.na(gpar[i, , g])], 2), collapse = ", ")
            k_values <- NULL
            for (k in 2:ncat) {
                if (!is.na(gpar[i, k, g])) 
                  k_values[[k]] <- gpar[i, k, g]
            }
            gk_values[[g]] <- k_values
        }

        for (g in 2:x$ng) {
            labi <- labels[g]
            gvalues <- mutate(gvalues, "{labi}_p4" :=
                (gdensity[, g] * abs(sumpp(pp[, i, 1:ncat, 1]) - sumpp(pp[, i, 1:ncat, g]))))
        }

        ivalues[[i]] <- list("data" = gvalues, "plot_names" = table_info,
            "chart2_text" = text_tibble, "chart_3_text" = g_values, 
            "chart_3_rug" = gk_values)
        
    }
    names(ivalues) <- difitems
    return(ivalues)


}

# # not working
# all_dif_data <- function(x, labels = c("Younger", "Older")) {
    
#     theta <- seq(x$options$minTheta, x$options$maxTheta, x$options$inc)
#     gpar <- array(NA, c(ndif, maxcat, x$ng))
#     cpar <- as.matrix(x$ipar.sparse[1:(x$ni - ndif), ])

#     plot_names <- list(xlab = "theta", ylab = "TCC", title = "All Items")

#     gvalues <- tibble(
#         y = (seq(0, sum(!is.na(x$ipar)) - x$ni, along = theta)),
#         x = theta
#     )
    
#     for (g in 1:x$ng) {
#         labi <- labels[g]
#         apar <- rbind(cpar, gpar[, , g])
#         gvalues <- mutate(gvalues, "{labi}" := pp[, i, k, g])
#         tcc(apar[, 1], apar[, -1, drop = F], theta, 
#             model = x$options$model)

#     }
#     return(list(data = gvalues, plot_names = plot_names))
   
# }


initial_purified_data <- function(x, labels = c("Younger", "Older")) {
    dat_boxplot <- tibble(
        "boxplot" = x$calib$theta - x$calib.sparse$theta,
    )

    dat = NULL
    for (i in 1:x$ng) {
        labi <- labels[i]
        xvalue <- x$calib$theta[x$group == as.numeric(names(table(x$group))[i])]
        yvalue <- dat_boxplot$boxplot[x$group == as.numeric(names(table(x$group))[i])]
        dati <- tibble(x = xvalue, y = yvalue)
        colnames(dati) <- c(paste0(labi, "_x"), paste0(labi, "_y"))
        dat[[labi]] <- dati
    }
    return(list(boxplot = dat_boxplot, data = dat, labels = c(x = "initial theta", y = "initial - purified")))

}
