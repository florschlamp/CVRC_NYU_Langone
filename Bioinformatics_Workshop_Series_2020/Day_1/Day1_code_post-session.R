## 
## Friday September 25th
## Day 1: “Introduction: First Steps on using R in Rstudio”
##



## How to navigate RStudio
# (console, environment panel, etc.)

## How to interact with the console
# Do some basic math examples

## Making / assigning variables
# Do some math examples using the variables instead

#### NOTE
# In the console you can run code by just pressing 'enter' in your keyboard
# In a file, lines need to be run by pressing "command/control"+"enter"
####


x <- 2 

x

x*3

x+5

y = x+5

y


Variable_A <- 45
Variable_B = 10
# what happens if you call a variable that doesn't exist?
Variable_C


## R scripts
# keeping records of everything you did
# (but also to use it again later)
# BUT lines need to be run by pressing "command/control"+"enter"

WT <- 3
KO <- 6

difference_WT_vs_KO <- KO - WT
difference_WT_vs_KO

# change values of WT and KO
# re-run formulas
# highlight multiple lines to run them all at same time
# (always using "Ctrl+Enter")

# Vectors! another type of variable:
my_numbers <- c(1,2,3,4,5)
my_numbers

my_numbers*10
my_numbers+1

my_numbers_multiplied_by_100 <- my_numbers*100
my_numbers_multiplied_by_100

# How to remove variables
rm(Variable_A) # one at a time
rm(Variable_B,x,y) # multiple ones
# how to clear the entire environment panel?



## Types of data
# numeric
A <- 0.01
class(A)

my_list <- c(1,2,3,4,5)
my_list
class(my_list)

# character
desired_color <- "blue"
desired_color
class(desired_color)

my_colors <- c("red","blue","yellow")
my_colors
class(my_colors)


# what happens when a number is coded as character?

x = 5
y = "5"
x+1
y+1
class(x)
class(y)
# re-assign a variable as numeric
y <- as.numeric(y)
y



# data frames (table)
my_table <- data.frame(fruit=c("apple","blueberry","banana"),
                       color=c("red","blue","yellow"))
my_table
class(my_table) # comment on this

## Ways to navigate data (vs. values)




## Functions
my_list

class(my_list)

length(my_list)
max(my_list)
min(my_list)
mean(my_list)
summary(my_list)

# let's make a more complex list of numbers
my_list <- c(22,35,7.8,0.01,365,12,65,78.98,13)
length(my_list)
max(my_list)
min(my_list)
mean(my_list)
summary(my_list)




my_list[1] 
# numbers inside square brackets "[]" represent index / position in the list
my_list[5]
my_list[10]
length(my_list) 

end_of_list <- length(my_list)
my_list[end_of_list]

my_list[9]

my_list[length(my_list)]

head(my_list)
tail(my_list)
head(my_list,5)
tail(my_list,1)


# Functions for data frames

length(my_table) # 2 ??
dim(my_table) # dim = dimension
# rows vs columns
nrow(my_table) # nrow = number of rows
ncol(my_table)

colnames(my_table)
rownames(my_table)
rownames(my_table) <- c("fruit1","fruit2","fruit3")
my_table




## <<< BREAK >>>>



# Set working directory to file location
setwd("~/Desktop/Day 1")


## Read a file into R
mammals <- read.csv("mammals.csv")
dim(mammals)
# or use 'Import Dataset' in the Environment Panel

head(mammals)
tail(mammals)

nrow(mammals)
ncol(mammals)

colnames(mammals)

mammals$species
mammals$body

list_of_body_weights <- mammals$body
mean(list_of_body_weights)


max(mammals$body)
min(mammals$body)

summary(mammals$body) # kg
summary(mammals$brain) # different units! # grams

# Let use better column names
colnames(mammals)
colnames(mammals) <- c("species","body_weight_kg","brain_weight_g")
head(mammals)

# Change units to kg
mammals$brain_weight_kg <- mammals$brain_weight_g/1000
head(mammals)

# create and delete columns
mammals$newcolumn <- 0
head(mammals)
mammals$newcolumn <- NULL

colnames(mammals)
dim(mammals)



# how to grab data from columns
mammals$body_weight_kg
mammals[,2]
mammals[,"body_weight_kg"]

# data.frame[desired_rows,desired_columns]

mammals[32,]

mammals$species == "Human"

mammals[mammals$species == "Human",]

is.human <- mammals$species == "Human"
mammals[is.human,]

# "==" is an 'if' function, "is it equal?"
# answer is either True or False 
# (also known as a boolean function)


mammals[c(14,55,61),]



my_animals_of_interest <- c("Lesser short-tailed shrew","Musk shrew","Tree shrew")

mammals[mammals$species %in% my_animals_of_interest,]

mammals[mammals$species %in% my_animals_of_interest,c("species","body_weight_kg")]
mammals[c(14,55,61),c(1,2)]




# How many mammals weigh more than 500 kg?
is.heavy <- mammals$body_weight_kg > 500
mammals[is.heavy,]

mammals[mammals$body_weight_kg > 500,]

# if you want just the names:
mammals[mammals$body_weight_kg > 500,1]
mammals[mammals$body_weight_kg > 500,'species']
mammals[mammals$body_weight_kg > 500,]$species

## Test what you have learned!!
# How many mammals have a brain smaller than 1 g?

# < try your own code here >
mammals[mammals$brain_weight_g<1,]
subset_table_small_brains <- mammals[mammals$brain_weight_g<1,]

# if you want to export data you just filtered:
write.csv(subset_table_small_brains, file="mammals_small_brain.csv") 

is.light <- mammals$brain_weight_g < 1
mammals[is.light,]


head(mammals)
mammals$newcolumn <- 0
rm(mammals$newcolumn)


## how to save the updated data frame
write.csv(mammals, file="mammals_updated.csv")

# read new data back again
new_table  <- read.csv("mammals_updated.csv")



## install packages
# first try through Packages tab in the right panel

# installing with code: 
install.packages("ggplot2")


# if there is no autocomplete check in bioconductor
# https://www.bioconductor.org/

# code for bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")

## THE END!
# Make sure you save your files!
# You don't need to save the 'workspace'

#### Homework for next session:
# install the following packages: DESeq2, edgeR, ggplot2


