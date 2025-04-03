 

dat <- read.csv("/husky/cellbuster/pountain2024/GSM6725354_D1_lib5_counts.txt.gz",sep="\t")
#head(dat)
#colnames(dat)
#rownames(dat)
dat <- dat[,-1]



df <- data.frame(
  cnt=sort(colSums(dat), decreasing = TRUE)
)

df$rank <- 1:nrow(df)

ggplot(df, aes(log10(rank),1+log10(cnt))) + geom_point()

dim(dat)
