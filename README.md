# Markov chain analysis of COVID-19 in Nigeria

<p align="center">
  <img src="images/covid19-thumbnail.png" alt="Markov Chain Analysis of COVID-19 in Nigeria" width="100%">
</p>


This project uses a discrete-time Markov chain to examine the early COVID-19 outbreak in Nigeria, covering the first 34 days of reported cases. The analysis organizes the data into Susceptible (S), Infected (I), and Recovered (R) states and constructs a transition probability matrix.

## Data

The analysis uses daily Nigerian COVID-19 records containing cumulative confirmed cases. Daily new infections are calculated as the difference between consecutive cumulative confirmed counts.

The working dataset includes these variables:

- `Day`: observation day, from 1 to 34.
- `Confirmed_Nigeria`: cumulative confirmed cases.
- `new_infections`: daily increase in confirmed cases.
- `S_t`: reconstructed susceptible count.
- `I_t`: daily new infections, used as the infected-state proxy in this analysis.
- `R_t`: reconstructed residual state.
- `new_recoveries`: change in the reconstructed `R_t` value.

The maximum cumulative confirmed count in the 34-day period is 873. This value is used as the analytical total (`N_effective`) for the reconstruction; it is not Nigeria's actual population.

## State construction

The states are reconstructed as follows:

```r
N_effective <- max(covid_19$Confirmed_Nigeria, na.rm = TRUE)

covid_19 <- covid_19 %>%
  mutate(
    new_infections = c(
      Confirmed_Nigeria[1],
      diff(Confirmed_Nigeria)
    ),
    S_t = N_effective - Confirmed_Nigeria,
    I_t = new_infections,
    R_t = Confirmed_Nigeria - new_infections
  )

covid_19 <- covid_19 %>%
  mutate(
    new_recoveries = c(NA, diff(R_t))
  )
```

## Objective 1: Transition matrix with an assumed recovery probability

For the first objective, the infection transition probability is estimated from daily new infections relative to the previous day's susceptible count. The recovery transition probability is set using an assumed 14-day average infectious period.

```r
p_values <- covid_19$new_infections[2:nrow(covid_19)] /
  covid_19$S_t[1:(nrow(covid_19) - 1)]

p <- mean(p_values, na.rm = TRUE)
q <- 1 / 14
```

The transition matrix uses the state order Susceptible, Infected, and Recovered:

```r
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
```

The matrix represents transitions from Susceptible to Infected, Infected to Recovered, and an absorbing Recovered state.

## Objective 2: Estimating both transition probabilities

The second objective is to estimate both infection and recovery transition probabilities from observed data rather than fixing the recovery probability at `1/14`.

**Data limitation:** the available dataset does not contain observed daily discharged or recovered cases. In the reconstructed table, `new_recoveries` is calculated as the change in `R_t`. Given the definition of `R_t`, this derived value equals the previous day's `I_t` (after the initial observation). Therefore, using

```r
q_t <- new_recoveries_t / I_previous
```

would mechanically produce a recovery probability of 1. It would not be an independently estimated recovery probability based on observed discharges.

To complete Objective 2 as an empirical estimation, daily discharge data from a reliable source are required. Any alternative estimate must be clearly identified as an assumption or proxy.

## R packages

The analysis uses the following R packages:

```r
install.packages(c("dplyr", "magrittr", "markovchain", "ggplot2"))
```

Load them with:

```r
library(dplyr)
library(magrittr)
library(markovchain)
library(ggplot2)
```

## Validation

Check that the transition probabilities are between 0 and 1 and that each row of the transition matrix sums to 1.

```r
round(P_matrix, 6)
rowSums(P_matrix)
```

## Limitations

- The analysis covers only the first 34 days of the outbreak.
- `N_effective = 873` is the final cumulative confirmed count in this dataset, not the national population.
- `I_t` is defined as daily new infections, rather than the conventional SIR count of all currently infectious individuals.
- `R_t` and `new_recoveries` are reconstructed variables, not observed discharge records.
- The 14-day recovery period in Objective 1 is a modelling assumption.
- The current data do not support an independent empirical estimate of the recovery transition probability for Objective 2.

## Project files

A suggested repository structure is:

```text
nigeria-covid-markov-analysis/
├── README.md
├── data/
│   └── Nigeria_COVID-19.xlsx
├── scripts/
│   ├── data_preparation.R
│   ├── objective_1_markov_chain.R
│   └── objective_2_markov_chain.R
└── outputs/
    └── figures/
```

Only upload the dataset if you have permission to redistribute it.

## Author

**Stanley Agbam**
