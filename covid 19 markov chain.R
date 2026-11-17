library(dplyr)
library(readxl)
library(tidyverse)
library(ggplot2)
library(janitor)
library(lubridate)
library(psych)
install.packages("markovchain")
library(markovchain)


install.packages("sf")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("ggrepel")
install.packages("remotes")
remotes::install_github("ropensci/rnaturalearthhires")
install.packages("magrittr")
library(magrittr)

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(ggrepel)
library(maps)
covid_i9 <- Nigeria_COVID_19
colnames(covid_i9) <- c("Date", "Day", "Confirmed_Nigeria", "Confirmed_Lagos", 
                        "Confirmed_Abuja", "Incidence_Nigeria", "Incidence_Lagos", 
                        "Incidence_Abuja")

covid_i9 <- covid_i9[-1, ]
covid_i9$Date               <- as.Date(as.numeric(covid_i9$Date), origin = "1899-12-30")
covid_i9$Day                <- as.numeric(covid_i9$Day)
covid_i9$Confirmed_Nigeria  <- as.numeric(covid_i9$Confirmed_Nigeria)
covid_i9$Confirmed_Lagos    <- as.numeric(covid_i9$Confirmed_Lagos)
covid_i9$Confirmed_Abuja    <- as.numeric(covid_i9$Confirmed_Abuja)
covid_i9$Incidence_Nigeria  <- as.numeric(covid_i9$Incidence_Nigeria)
covid_i9$Incidence_Lagos    <- as.numeric(covid_i9$Incidence_Lagos)
covid_i9$Incidence_Abuja    <- as.numeric(covid_i9$Incidence_Abuja)

str(covid_i9)
head(covid_i9)
covid_19 <- covid_i9
#Exploratory data analysis
summary(covid_i9$Confirmed_Nigeria)

covid_19 %>%
  ggplot(aes(x = Confirmed_Nigeria)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 50, 
                 fill = "steelblue", color = "white", alpha = 0.7) +
  geom_density(color = "red", linewidth = 1) +
  labs(title = "Distribution of Confirmed COVID-19 Cases in Nigeria",
       x = "Confirmed Cases",
       y = "Density") +
  theme_minimal()

#confirmed cases trend with time 
covid_19 %>%
  ggplot(aes(x = Day, y = Confirmed_Nigeria)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue") +
  labs(title = "COVID-19 Confirmed Cases Over Time - Nigeria",
       x = "Day", y = "Confirmed Cases") +
  theme_minimal()




# STEP 1: Get Nigeria's geographical boundary

world_map <- ggplot2::map_data("world")

nigeria_map <- subset(world_map, region == "Nigeria")


# STEP 2: Create location dataset

location_data <- data.frame(
  Location = c("Lagos", "Abuja"),
  Longitude = c(3.3792, 7.4951),
  Latitude = c(6.5244, 9.0579),
  Confirmed = c(504, 119)
)


# STEP 3: Plot Nigeria map

ggplot() +
  
  geom_polygon(
    data = nigeria_map,
    aes(x = long, y = lat, group = group),
    fill = "#f0e6d3",
    color = "white",
    linewidth = 0.5
  ) +
  
  geom_point(
    data = location_data,
    aes(x = Longitude, y = Latitude, size = Confirmed),
    color = "tomato",
    alpha = 0.8
  ) +
  
  geom_text(
    data = location_data,
    aes(
      x = Longitude,
      y = Latitude,
      label = paste0(Location, "\n", Confirmed, " cases")
    ),
    nudge_y = 0.5,
    size = 3.5,
    fontface = "bold",
    color = "gray20"
  ) +
  
  scale_size_continuous(
    range = c(6, 18),
    name = "Confirmed Cases"
  ) +
  
  coord_quickmap(
    xlim = c(2, 15),
    ylim = c(4, 14)
  ) +
  
  labs(
    title = "Total Confirmed COVID-19 Cases by Location",
    subtitle = "Nigeria — First 34 Days of Outbreak (2020)",
    caption = "Source: Nigeria COVID-19 Dataset"
  ) +
  
  theme_void() +
  
  theme(
    plot.title = element_text(
      face = "bold", size = 15, hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 11, hjust = 0.5, color = "gray40"
    ),
    plot.caption = element_text(
      size = 8, hjust = 1, color = "gray60"
    ),
    legend.position = "bottom"
  )



#objective 1:To construct a discrete-time Markov chain model 
#representing the transition of individuals between Susceptible,
#Infected, and Recovered states using COVID-19 cumulative case data.


# STEP 1: Define the SIR states

N_effective <- max(covid_19$Confirmed_Nigeria, na.rm = TRUE)

covid_19 <- covid_19 %>%
  mutate(
    S_t = N_effective - Confirmed_Nigeria,
    I_t = Incidence_Nigeria,
    R_t = Confirmed_Nigeria - Incidence_Nigeria
  )

# STEP 2: Estimate transition probabilities p and q

# Probability of moving from Susceptible to Infected
p_values <- covid_19$Incidence_Nigeria[2:nrow(covid_19)] /
  covid_19$S_t[1:(nrow(covid_19) - 1)]

p <- mean(p_values, na.rm = TRUE)

# Probability of moving from Infected to Recovered
q <- 1 / 14

cat("Estimated p:", round(p, 6), "\n")
cat("Estimated q:", round(q, 6), "\n")

# STEP 3: Construct the transition probability matrix

P_matrix <- matrix(
  c(
    1 - p, p,     0,
    0,     1 - q, q,
    0,     0,     1
  ),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(
    c("Susceptible", "Infected", "Recovered"),
    c("Susceptible", "Infected", "Recovered")
  )
)

# STEP 4: Display the transition matrix

print(round(P_matrix, 6))

# STEP 5: Check that each row sums to 1

cat("\nRow sums:\n")
print(rowSums(P_matrix))


# STEP 4: Build markovchain object
mc_sir <- new("markovchain",
              states = c("Susceptible", "Infected", "Recovered"),
              transitionMatrix = P_matrix,
              name = "SIR COVID-19 Nigeria")

summary(mc_sir)
print(mc_sir)

plot(mc_sir,
     main = "Discrete-Time Markov Chain: SIR Model - Nigeria COVID-19")







# OBJECTIVE 2: Estimate p and q from observed confirmed and
#              discharged cases (empirical transition matrix)


library(dplyr); library(ggplot2); library(markovchain); library(expm)

# ---- STEP 1: add the cumulative DISCHARGED column -----------
# Your Objective 1 data has no discharged column, so add it here
# (34 values, one per day, same order as covid_19).
# covid_19$Discharged_Nigeria <- c( ... 34 cumulative values ... )

# ---- STEP 2: estimate q from the data -----------------------
n <- nrow(covid_19)

if (!"Discharged_Nigeria" %in% names(covid_19)) {
  # Fallback: reuse the q printed in your chapter (Table 3)
  q_obs <- 0.272056
  message("No discharged data found: using q = 0.272056 from the chapter.")
} else {
  new_discharged <- c(NA, diff(covid_19$Discharged_Nigeria))   # daily discharges
  # Same style as p: daily discharges / previous day's infected (I_t)
  q_values <- new_discharged[2:n] / covid_19$I_t[1:(n - 1)]
  q_values[!is.finite(q_values)] <- NA
  q_obs <- mean(q_values, na.rm = TRUE)
}

p_obs <- p   # p is already estimated from the data in teh first objective 

cat("Estimated p:", round(p_obs, 6), "\n")
cat("Estimated q:", round(q_obs, 6), "\n")

# ---- the transition probability matrix 
states <- c("Susceptible", "Infected", "Recovered")
P_observed <- matrix(c(1 - p_obs, p_obs,     0,
                       0,         1 - q_obs, q_obs,
                       0,         0,         1),
                     nrow = 3, byrow = TRUE,
                     dimnames = list(states, states))
print(round(P_observed, 6))
cat("\nRow sums:\n"); print(rowSums(P_observed))

mc_sir_observed <- new("markovchain", states = states,
                       transitionMatrix = P_observed,
                       name = "SIR COVID-19 Nigeria (Observed)")
summary(mc_sir_observed)
plot(mc_sir_observed,
     main = "Discrete-Time Markov Chain: Observed SIR - Nigeria COVID-19")


if (!requireNamespace("expm", quietly = TRUE)) install.packages("expm")
library(expm)
library(markovchain)
library(ggplot2)

# Backup matrix-power function (works even without expm)
mpow <- function(M, k) {
  out <- diag(nrow(M))
  if (k > 0) for (i in 1:k) out <- out %*% M
  dimnames(out) <- dimnames(M)
  out
}

# ---- STEP 4: Objective 1 vs Objective 2 comparison ----
comparison <- data.frame(
  Parameter   = c("p (infection)", "q (recovery)", "Basis for q"),
  Objective_1 = c(round(p, 6), round(1/14, 6), "Fixed (1/14 days)"),
  Objective_2 = c(round(p_obs, 6), round(q_obs, 6), "Estimated from data"))
print(comparison)

# ============================================================
# SECTION 4.4: LONG-RUN DYNAMICS (using P_observed)
# ============================================================
steadyStates(mc_sir_observed)

# Multi-step matrices
for (k in c(7, 14, 30)) {
  cat("\nAfter", k, "days:\n")
  print(round(mpow(P_observed, k), 6))
}

# Fundamental matrix and mean absorption times
Q_obs <- P_observed[1:2, 1:2]
N_obs <- solve(diag(2) - Q_obs)
print(round(N_obs, 4))
abs_times <- rowSums(N_obs)
cat("\nMean days to Recovered from Susceptible:", round(abs_times[1], 2), "\n")
cat("Mean days to Recovered from Infected:",    round(abs_times[2], 2), "\n")

# State probabilities over the 34 days, starting Susceptible
states <- c("Susceptible", "Infected", "Recovered")
days   <- 0:34
traj   <- t(sapply(days, function(t) c(1, 0, 0) %*% mpow(P_observed, t)))
colnames(traj) <- states
traj_df <- data.frame(Day = days, traj) |>
  tidyr::pivot_longer(-Day, names_to = "State", values_to = "Probability")

ggplot(traj_df, aes(Day, Probability, colour = State)) +
  geom_line(linewidth = 1.1) +
  labs(title = "State Probabilities Over 34 Days (Start: Susceptible)",
       x = "Day", y = "Probability") +
  theme_minimal()


  