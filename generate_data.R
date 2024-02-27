if (!all(c("CptPointProcess", "remotes")%in%installed.packages())){
    install.packages(c("remotes"))
    remotes::install_local("CptPointProcess_0.0.0.1.tar.gz")
}
require("CptPointProcess")
require("parallel")

generate_proc_pois_h <- function(lambda, nb_event) {
    cumsum(rexp(nb_event, rate = lambda))
}

generate_proc_pois_cpm <- function(lambdas, nb_event_vec) {
    event_times <- NULL
    df <- data.frame(lambdas = lambdas, nb_event_vec = nb_event_vec)
    for(row_id in seq_len(nrow(df))) {
        lambda <- df[row_id, 1]
        nb_event <- df[row_id, 2]

        # Génère les temps selon le lambda courant
        proc_pois_seg <- generate_proc_pois_h(lambda = lambda, 
            nb_event = nb_event)
        
        # Détecte le dernier temps du processus précédent
        last_time <- ifelse(!is.null(event_times[length(event_times)]), 
            event_times[length(event_times)], 0)
        # Allonge le vecteur des temps en ajoutant au vecteur les temps générer
        event_times <- c(event_times, proc_pois_seg + last_time)
    }
    return(event_times)
}

full_lambdas <- c(seq(2,20, by = 2))
set.seed(1234)
one_simulation <- function(lambdas, nb_event_vec){
    data <- generate_proc_pois_cpm(lambdas = lambdas, nb_event_vec = nb_event_vec)
    max_time <- data[length(data)]
    data <- data/max_time
    if (sum(nb_event_vec) <= 10){
        K <- sum(nb_event_vec) - 3
    } else {
        K <- 10
    }

    result <- CptPointProcess(ProcessData = data.frame(times = data[-length(data)]), Kmax = K)
    # Augmente la barre
    # pb$tick()
    K_hat <- result$K.est
    if (K_hat == length(lambdas)){
        estim_lambdas <- result$SegK$lambda * max_time
        ecart_lambdas <- sum((lambdas - estim_lambdas))^2
    } else {
        ecart_lambdas <- NA
    }
    return(data.frame(N_per_group = sum(nb_event_vec)/length(lambdas), real_K = length(lambdas), K_hat = K_hat, correct_estim = K_hat == length(lambdas), ecart_lambdas = ecart_lambdas))
}

generate_dataset <- function(full_lambdas, max_delta_N = 200, nb_rep = 10, nb_sample = 10) {
    delta_Ns <- seq(20, max_delta_N, by = 10)
    nb_lambdas <- seq(2, length(full_lambdas))

    conditions <- data.frame(delta_N = rep(delta_Ns, nb_rep * nb_sample * length(nb_lambdas)),
        nb_lambda = rep(nb_lambdas, each =  nb_rep * nb_sample * length(delta_Ns))
        )

    return(do.call("rbind", parallel::mclapply(seq_len(nrow(conditions)), function(row_id){
        # Extraction des conditions
        nb_lambda <- conditions[row_id, 2]
        delta_N <- conditions[row_id, 1]

        cat(paste("nb_lambda", nb_lambda, "delta_N", delta_N))

        # Construction des paramètres
        nb_event_vec <- rep(delta_N, nb_lambda)
        lambda_vec <- sample(full_lambdas, nb_lambda)
        res <- one_simulation(lambdas = lambda_vec, nb_event_vec = nb_event_vec)
        res$delta_N <- delta_N
        res$nb_lambda <- nb_lambda
        return(res)
    }, mc.cores = parallel::detectCores() - 1)))

}

data <- generate_dataset(full_lambdas = full_lambdas)
save(data, file = "msv-sem.Rds")