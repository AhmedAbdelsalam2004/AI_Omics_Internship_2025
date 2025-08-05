# Creating sub-folders

dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")


#Load dataset

data <- read.csv(file.choose())

# Check structure

View(data)
str(data)
summary(data)


# Identify variables with incorrect or inconsistent data types.

# - The two BMIs corresponding to patient ids 16 and 17 
# - are integers not numerical, but it looks like R studio auto fixes this.  

#  Convert variables to appropriate data types where needed

data$gender = as.factor(data$gender)
data$diagnosis = as.factor(data$diagnosis)
data$smoker = as.factor(data$smoker)
str(data)


# Create a new variable for smoking status as a binary factor:
# 1 for "Yes", 0 for "No"

data$smoker_fac <- ifelse(data$smoker == "Yes", 1, 0)
str(data)

# Save the cleaned dataset in your clean_data folder with the name patient_info_clean.csv

write.csv(data, file = "C:/Users/a3abd/OneDrive/Desktop/AI_Omics_Internship_2025/Module_I/clean_data/patient_info_clean.csv")


# Save your R script in your script folder with name "class_Ib"

# Save the entire R workspace
save.image(file = "AhmedAbdelsalam_Class_Ib_Assignment.RData")

