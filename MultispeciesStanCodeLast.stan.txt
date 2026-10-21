/* HMM model for seed production at Mt. Rainier (multispecies)
Eléonore Perret + Victor Van der Meersch + Elizabeth Wolkovich + Janneke Hille Ris Lambers

August 2026 (Last data modified)
July 2026 (code last modified)
October 2026 (comments updated by JHRL, anything in caps and parentheses e.g. (CHECK) denotes something JHRL wants to verify or a question)

Data: we have seed production by multiple species (7) in multiple stands (18) over multiple years (2009 - 2024) at MORA.

Data specification: 
- Not all stands were surveyed starting in the same year and some were abandoned before 2024, so we use an ragged vector to establish the start and end observation for a specific stand. 
- One stand (SUNR) had no data for 2023 (all samples lost), so we added missing observations for the 3 species in this plot (CHECK)
- We pooled (summed) seed counts from all traps present at a site which were then divided by the total area of seed trap (see transformed data + negative binomial model); this area varied by year / stand due to differing sizes of traps and missing traps

Model: 
- We assume seed production (also called emission) is distributed as a negative binomial (overdispersed counts). 
- We assume two latent states (the HMM); low seed production (state 1) and high seed production (state 2 - the masting state). 
- This can be described by a 2x2 matrix defined by two parameters; the probability of remaining in state 1 when in state 1 (theta1 with 1 - theta1 the probability of transitioning to state 2), and the probability of staying in state 2 when in state 2 (theta2 - with 1 - theta2 the probability of switching to state 1)
- This also requires parameter rho, the probability of starting in state 1 or state 2 (this is for the first year of data)

Partial pooling:
- We assume partial pooling for species and stand seed production means. Specifically:
	- Low state partially pooled = global low mean + species effect + species specific stand effect
	- High state partially pooled = global low mean + delta (dif masting vs low) species effect + species specific stand effect
- We assume partial pooling for theta1 and theta2. Specifically:
	- Theta 1 (prob staying in low state) = global theta1 + species effect (but no stand effect)
	- Theta 2 (prob staying in high state) = global theta 2 + species effect (but no stand effect)
- We assume partial pooling of rho by species only (CHECK)
- We assume partial pooling of variance terms describing mean seed production (variation among species, stands) (CHECK)

Model stability:
- We used visual plots to diagnose funnelling and centered accordingly
- We centered... (CHECK and specify)
	- area surveyed (which varied year to year) by a baseline, defined as the minimum area surveyed (see transformed data)
- To ensure that state 1 is always lower than state 2, we specified delta as the difference between the two and forced it to be >0 by fitting on the log scale (CHECK)


Inference: 
- We are interested in seed productions means (by species and stand) in low and high states (defining boom vs. bust, difference)
- We are also interested theta1 and theta2 (also by species and stand), the transitions from low to high and high to low (defining masting frequencies).
- We are also interested in synchrony (between and within species), which we can estimate from species and stand-specific latent states

 */


//DATA

data {
  int<lower=1> N;                                   // total number of observations
  int<lower=1> F;                                   // number of species x stand series
  int<lower=1> S;                                   // number of species
  int<lower=1> N_stands;                            // number of unique stands

  array[N] int<lower=0> y;                          // Seed counts

  array[F] int<lower=1> start_idxs;                 // defines start of seed data for a speciesxstand combo (CHECK)
  array[F] int<upper=N> end_idxs;                   // defines end of seed data for a speciesxstand combo (CHECK)
  array[F] int<lower=1, upper=N_stands> stand_id;   // numeric id of stand (CHECK)
  array[F] int<lower=1, upper=S> species_id;        // numeric id of species (CHECK)

  vector<lower=0>[N] area;                          // area of traps surveyed for a particular species-stand-year combination
  
  array [N] int <lower= 0, upper = 1> obs_missing;  // vector identifying whether observations are missing (CHECK)
}


//TRANSFORMED DATA

transformed data {                                  // for stability, everything fit as a function of area / min(area)
  real baseline_area = min(area);
  vector[N] log_area_ratio;
  for (t in 1:N)
    log_area_ratio[t] = log(area[t] / baseline_area);
}


//PARAMETERS (these are what STAN updates during sampling)

parameters {
  // Initial state probability  
  simplex[2] rho;                                  // 2x2 matrix to define initial state (QUESTION - by species/ stand same for all?)

  // Transitions (pooled by species and non-centered) 
  real                   grand_logit_theta1;       // global mean of theta1 (QUESTION - logit in name, or logit transformation happen here?)
  real                   grand_logit_theta2;       // global mean of theta2

  vector[S]              alpha_theta1_species_nc;  // species deviation from the global theta 1, non-centered (QUESTION - should this be centered?)
  vector[S]              alpha_theta2_species_nc;  // species deviation from the global theta 2, non-centered (QUESTION - should this be centered?)
  real<lower=0>          sigma_theta1_species;     // variation among species in theta 1
  real<lower=0>          sigma_theta2_species;     // variation among species in theta 2

  // Emission means (pooled by species and stand) 
  // Low state: grand mean + species effect + stand effect (centered)
  real                   mu_log_low;              // global mean of seed production
  vector[S]              alpha_low_species;       // species deviation from the low state mean
  real<lower=0>          sigma_low_species;       // variation among species in low state mean
  matrix [S, N_stands]   alpha_stand;             // species-specific stand deviation from the low state mean
  real <lower=0>         sigma_stand;             // variation among stands (QUESTION - species specific?)

  // Delta (log-scale gap between states,): grand mean + species effect + stand effect (non-centered)
  real                   mu_log_delta;            // difference between low and high state seed production
  vector[S]              log_delta_species_nc;    // deviation per species in this difference, non-centered
  real<lower=0>          sigma_log_delta_species; // variation among species
  matrix [S, N_stands]   log_delta_stand_nc;      // species-specific deviation per stands non-centered  
  real <lower=0>         sigma_delta_stand;       // variation among stands  (QUESTION - species specific?)

  // Dispersion: (for negative binomial seed production), per state and species
  vector [S]             log_phi_species;         // species specific overdispersion
  vector [2]             log_phi_state;           // state specific overdispersion (QUESTION - species + state?)
}


//TRANSFORMED PARAMETERS (these are direct transformations, often for stability (CHECK)

transformed parameters {                                  
  // Non-centered species effects on transition probabilities
  vector[S] alpha_theta1_species = sigma_theta1_species * alpha_theta1_species_nc; // centers the species-specific differences from global theta1
  vector[S] alpha_theta2_species = sigma_theta2_species * alpha_theta2_species_nc; // centers the species-specific differences from global theta2

  // Non-centered species and stand effects on delta seed production (difference between low and high seed production state)
  vector[S] log_delta_species    = sigma_log_delta_species * log_delta_species_nc;
  matrix[S, N_stands] log_delta_stand   = sigma_delta_stand * log_delta_stand_nc;
  
  // Per-species and stand emission means (seed production)
  vector[F] log_alpha_low;
  vector[F] log_alpha_high;

  for (f in 1:F) {
    int s  = species_id[f];                            // s is the species index - get from species_id (CHECK)
    int st = stand_id[f];                              // st is the stand index - get from stand_id (CHECK)

    log_alpha_low[f]  = mu_log_low                     // global mean seed production
                      + alpha_low_species[s]           // add the species effect
                      + alpha_stand[s,st];             // add the species-specific stand effect

    log_alpha_high[f] = log_alpha_low[f]               // seed production at low state
                  + log1p_exp(mu_log_delta)            // add the difference between low and high (between low and high, constrained to be > 0)
                  + log1p_exp(log_delta_species[s])    // add a species specific effect (QUESTION - should this be constrained to be > 0?)
                  + log1p_exp(log_delta_stand[s,st]);  // add a species specific stand effect (QUESTION - should this be constrained to be > 0?)
  }

  //  Per-series transition matrices
  vector<lower=0, upper=1>[F] theta1;
  vector<lower=0, upper=1>[F] theta2;
  array[F] matrix[2, 2] Gamma;                         // What is going on here?

  for (f in 1:F) {
    int s  = species_id[f];

    theta1[f] = inv_logit(grand_logit_theta1          // theta1 is global mean + species effect; note inv.logit constrains theta1 to be between 0 and 1
                         + alpha_theta1_species[s])
                         ;
    theta2[f] = inv_logit(grand_logit_theta2          // theta2 is global mean + species effect; note inv.logit constrains theta1 to be between 0 and 1
                         + alpha_theta2_species[s])
                      ;

    Gamma[f][1, 1] = theta1[f];                      // (QUESTION - why gamma?)
    Gamma[f][1, 2] = 1 - theta1[f];
    Gamma[f][2, 1] = 1 - theta2[f];
    Gamma[f][2, 2] = theta2[f];
  }

  // Emission log-likelihoods
  matrix[2, N] log_omega;                             // matrix denoting the probability of observing data if state 1 (col1) vs state 2 col2) (QUESTION)
  for (f in 1:F) {
    int s = species_id[f];  
    for (t in start_idxs[f]:end_idxs[f]) {
      if (obs_missing[t]){                           // Missing observations do not provide information 
        log_omega[1,t]= 0;
        log_omega [2,t] =0;
      } else {                                       // Likelihood of data given state 1, state 2? (QUESTION) 
        log_omega[1, t] = neg_binomial_2_log_lpmf(   // State 1 likelihood 
          y[t] | log_alpha_low[f]  + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[1])); // run through this (QUESTION) 
      log_omega[2, t] = neg_binomial_2_log_lpmf(     // State 2 likelihood
          y[t] | log_alpha_high[f] + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[2])); // run through this (QUESTION) 
      }
    }
  }
}

// NOTE FROM JANNEKE - up until here I mostly had time to comment / think, here is about where I ran out of time


//MODEL

model {
  // Initial state 
  rho                       ~ dirichlet(rep_vector(2.0, 2)); // repeat assumption of state 2 twice (QUESTION - not sure why twice??? or is this prior?)

  //Transitions (priors) (QUESTION)
  grand_logit_theta1        ~ normal(0.5, 1);    // prior for global theta1 (grand mean), on logit scale (set highish, prob of staying in low if low)
  grand_logit_theta2        ~ normal(-1, 0.5);   // prior for global theta2 (grand mean), on logit scale (set low, prob of staying in masting if nasting)

  alpha_theta1_species_nc   ~ normal(0, 1);
  alpha_theta2_species_nc   ~ normal(0, 1);
  sigma_theta1_species      ~ normal(0, 0.7);
  sigma_theta2_species      ~ normal(0, 0.7);

  // Emission means (mostly priors), (QUESTION)
  mu_log_low                ~ normal(2, 1.5);
  alpha_low_species         ~ normal(0, sigma_low_species);
  sigma_low_species         ~ normal(0.5, 1);
  for (s in 1:S) 
    alpha_stand [s]         ~ normal(0, sigma_stand);  // CHECK - need an end to for loop?
  sigma_stand               ~ normal(1, 0.5);

  mu_log_delta              ~ normal(1.5, 1.5);
  log_delta_species_nc      ~ normal(0, 1);
  sigma_log_delta_species   ~ normal(1.5, 0.7);
  for (s in 1:S)
    log_delta_stand_nc [s]  ~ normal(0, 1);
  sigma_delta_stand         ~ normal(1, 0.5); 

  // Dispersion
  log_phi_species           ~ normal(log(4), 0.6); 
  log_phi_state             ~ normal(0, 0.3); 

  // HMM marginal likelihood
  for (f in 1:F)
    target += hmm_marginal(
        log_omega[, start_idxs[f]:end_idxs[f]],
        Gamma[f],
        rho);
}


//GENERATED QUANTILES  (QUESTION) - probably could use more time or hand holding here

generated quantities {
  array[N] int<lower=0>          y_rep;
  array[N] int<lower=1, upper=2> state;

  for (f in 1:F) {
    int start_id = start_idxs[f];
    int end_id   = end_idxs[f];

    state[start_id:end_id] = hmm_latent_rng(
        log_omega[, start_id:end_id],
        Gamma[f],
        rho);

    for (t in start_id:end_id) {
  int s = species_id[f];  
  
  if (obs_missing[t]) {
    y_rep[t] = 0;
  } else if (state[t] == 1) {
    y_rep[t] = neg_binomial_2_log_rng(
        log_alpha_low[f]  + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[1]));
  }else{
    y_rep[t] = neg_binomial_2_log_rng(
        log_alpha_high[f] + log_area_ratio[t], exp(log_phi_species[s] + log_phi_state[2]));
      }
    }
  }
}

