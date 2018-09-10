# Day 1 - Intro to R

my_numbers = c(1, 4, 6)

my_characters = c("M", "E", "H")

my_numbers[1] + 3

my_logicals = c(TRUE, FALSE, TRUE, TRUE)


## Making a list
my_list = list(MyNumbers=my_numbers,
               MyLetters=my_characters)

## Retrieving from a list
my_list[1]
my_list[[1]]
my_list[["MyNumbers"]]
my_list$MyNumber


## Make a dataframe
my_df = data.frame(MyLetters=my_characters,
                   MyNumbers=rep(my_numbers, 2))
my_df

## Retrieve from a dataframe
my_df$MyNumbers
my_df$MyNumbers[1]
my_df[1, 2]

rownames(my_df) = c("M", "E", "H", "M2", "E2", "H2")

my_df

my_df["M", "MyNumbers"]


## Making a matrix

### E1. Create a matrix with 4 rows, 5 columns
### ..and name both rows and columns
m <- matrix(nrow=4, ncol=5)
rownames(m) <- c("A", "B", "C", "D")
colnames(m) <- c("A", "B", "C", "D", "E")

my_mat = matrix(nrow=4, ncol=5,
                dimnames=list(Row=LETTERS[1:4],
                              Col=letters[1:5]))


## Generating numbers
my_lucky_numbers = seq(1, 300, by=30)
my_normal_numbers = rnorm(10, mean=10, sd=2)
my_other_lucky_numbers = seq(1, 1000, length.out=5)

## Create a character vector..
## ..with every 3rd engligh alphabet letter


my_letters = letters[seq(1,26,3)]
my_other_other_lucky_numbers = seq(1,26,3)

names(my_other_other_lucky_numbers) = my_letters

my_other_other_lucky_numbers["m"]


## Frequent 2d-array functions
dim(my_mat)
colnames(my_mat)
rownames(my_mat)


## Handling characters
paste("Sample", LETTERS, sep="-")
attached_uppercase_letters = paste(LETTERS, collapse="")

english_alphabet_uppercase_letters = unlist(
  strsplit(attached_uppercase_letters, "", TRUE))



### Excercise of slide 26 
gene_mat <- matrix(runif(400), 100, 4)
colnames(gene_mat) <- paste0("sample", LETTERS[1:4])
# Alternatively:
gene_mat <- matrix(runif(400), 100, 4,
                   dimnames=list(NULL,
                                 paste0("sample", LETTERS[1:4])))
gene_mat[,1] <- gene_mat[,1] + 2
gene_mat[,2] <- gene_mat[,2] * 4


mean(gene_mat[,1])
mean(gene_mat[,2])
mean(gene_mat[,3])
mean(gene_mat[,4])

for(i in 1:ncol(gene_mat)){
  mean_column_i = mean(gene_mat[, i])
  my_output_string = paste(
    "Mean of column", i, "is", mean_column_i)
  print(my_output_string)
}

column_means = apply(gene_mat, 2, mean)


max_mean_of_columns = max(column_means)
which(column_means == max_mean_of_columns)


which.max(column_means)


## Excercise of slide 26 part 6
### Which column has the largest value
which.max(apply(gene_mat, 2, max))


# Loading files into R
virchip_df = read.csv("VirchipSuppTable.tsv", header=TRUE,
                      sep="\t")

options(stringsAsFactors=TRUE)
mice_df = read.csv("mice.csv")
summary(mice_df)

## head, dim, str, summary
summary(mice_df)

my_genotypes_str = c("AA", "Aa", "aa")
my_genotypes = factor(c("AA", "Aa", "aa"))

my_genes = factor(c("BRCA1", "TP53", "VHL"),
                  levels=c("VHL", "BRCA1", "TP53"))

## Use a for loop and a if block ..
## ..to print Young or old if age is higher than 9
for(i in 1:nrow(mice_df)){
  if(mice_df$Age[i] > 9){
    print("Old")
  }else{
    print("Young")
  }
}
