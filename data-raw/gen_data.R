library(haven)
library(here)
library(data.table)
df = read_dta(here("data-raw/stevenson.dta"))
df = as.data.table(df)

df[, judge_pre := 
  (1 * df$judge_pre_1) + (2 * df$judge_pre_2) +
  (3 * df$judge_pre_3) + (4 * df$judge_pre_4) + 
  (5 * df$judge_pre_5) + (6 * df$judge_pre_6) + 
  (7 * df$judge_pre_7) + (8 * df$judge_pre_8)
]

df[, t := 
  (1 * df$t1) + (2 * df$t2) +
  (3 * df$t3) + (4 * df$t4) + 
  (5 * df$t5) + (6 * df$t6)
]

df[, trial_time_of_day := fcase(
  morning == 1, "Morning",
  evening == 1, "Evening",
  grave == 1, "Grave",
  weekend == 1, "Weekend"
)]

df$judge_pre_1 = df$judge_pre_2 = df$judge_pre_3 = df$judge_pre_4 = df$judge_pre_5 = df$judge_pre_6 = df$judge_pre_7 = df$judge_pre_8 = NULL
df$t1 = df$t2 = df$t3 = df$t4 = df$t5 = df$t6 = NULL
df$morning = df$evening = df$grave = df$weekend = NULL
df$onePrior = df$threePriors = NULL 
df$sum =  df$F1 =  df$F2 =  df$F3 =  df$F =  df$M1 =  df$M2 =  df$M3 =  df$M = NULL


setcolorder(df, c("judge_pre", "jail3", "guilt", "bailDate", "white", "black", "prior_felChar", "priorCases", "prior_guilt"))

stevenson = as.data.frame(df)
usethis::use_data(stevenson, overwrite = TRUE)
