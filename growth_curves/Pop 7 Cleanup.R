library(gcplyr)
library(dplyr)
library(lubridate)

dat <- read_wides(list.files(path = "./Pop 7 Raw/", full.names = TRUE),
                  startrow = 41, startcol = "B",
                  metadata = list(starttime = c(7, "B")),
                  run_names_path = FALSE)

dat2 <- trans_wide_to_tidy(
  dat, id_cols = c("file", "starttime", "Time", "Tｰ 600"))
dat2 <- merge_dfs(dat2, collapse = TRUE, drop = TRUE)
                          
dat2 <- mutate(dat2, Time = as.numeric(hms(Time)))

pt1_rows <- which(dat2$file == "2024-02-02 William Growth Curve Pop 7 Part 1 Data")
pt2_rows <- which(dat2$file == "2024-02-02 William Growth Curve Pop 7 Part 2 Data")

dat2[pt2_rows, "Time"] <- dat2[pt2_rows, "Time"] +
  difftime(strptime(dat2[pt2_rows[1], "starttime"], format = "%I:%M:%S %p"),
           strptime(dat2[pt1_rows[1], "starttime"], format = "%I:%M:%S %p"),
           units = "secs")

write.csv(select(dat2, c("file", "Time", "Well", "Measurements")), 
          file = "2024-02-02 William Growth Curve Pop 7 Tidy.csv",
          row.names = FALSE)
