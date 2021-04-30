
# load data and functions

load("data/olfactory_bulb_11.Rdata")
source("functions/Boost_HP_function.R")

# filter data
filter_result <- filter_count(count, loc, min_total = 10, min_percentage = 0.1)
loc_f <- filter_result[[1]]
count_f <- filter_result[[2]]


start_time <- Sys.time()

# In this example, we just run first 100 genes in this dataset

detect_result <- Boost_Ising(count_f[,1:100], loc_f, normalization = 'tss', clustermethod = 'MGC')

end_time <- Sys.time()

print(end_time - start_time)

# display first 5 rows 
print(detect_result[1:5,])

# get SE genes
SE_gene <- rownames(detect_result)[which(detect_result$BF_neg > 150)
print(SE_gene)

# number of detected SE genes 
print(length(SE_gene))
