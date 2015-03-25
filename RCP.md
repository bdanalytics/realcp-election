# US Elections: Republican classification
bdanalytics  

**  **    
**Date: (Thu) Mar 26, 2015**    

# Introduction:  

Data: RealClearPolitics
Source: 
    Training:   https://courses.edx.org/c4x/MITx/15.071x_2/asset/PollingData.csv  
    New:        <newdt_url>  
Time period: 



# Synopsis:

Based on analysis utilizing <> techniques, <conclusion heading>:  

### ![](<filename>.png)

## Potential next steps include:

# Analysis: 

```r
rm(list=ls())
set.seed(12345)
options(stringsAsFactors=FALSE)
source("~/Dropbox/datascience/R/mydsutils.R")
source("~/Dropbox/datascience/R/myplot.R")
source("~/Dropbox/datascience/R/mypetrinet.R")
# Gather all package requirements here
#suppressPackageStartupMessages(require())

#require(sos); findFn("pinv", maxPages=2, sortby="MaxScore")

# Analysis control global variables
glb_is_separate_newent_dataset <- FALSE    # or TRUE
glb_split_entity_newent_datasets <- TRUE   # or FALSE
glb_split_newdata_method <- "condition"          # or "sample"
glb_split_newdata_condition <- "Year >= 2012"
# eval(expression(Year >= 2012), glb_entity_df, parent.frame())
# eval(substitute("Year >= 2012"), glb_entity_df, parent.frame())
# do.call("subset", list(glb_entity_df, expression(Year >= 2012)))

glb_split_newdata_size <- 0.25               # > 0 & < 1
glb_split_sample.seed <- 1000               # or any integer 

glb_predct_var <- "Republican"           # or NULL
glb_predct_var_name <- paste0(glb_predct_var, ".predict")
glb_id_vars <- c("State", "Year")                # or NULL

glb_exclude_vars_as_features <- union(glb_id_vars, ".rnorm")     # or NULL                      
# List chrs converted into factors; num/int transformed  
# glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, 
#                                       c("<col_name>")     # or NULL
#                                       )
# List feats that shd be excluded due to known causation by prediction variable
# glb_exclude_vars_as_features <- union(glb_exclude_vars_as_features, 
#                                       c("<col_name>")     # or NULL
#                                       )

glb_mice_complete.seed <- 144               # or any integer
glb_is_regression <- FALSE; glb_is_classification <- TRUE

glb_mdl <- glb_sel_mdl <- glb_dmy_mdl <- NULL
glb_models_df <- data.frame()

script_df <- data.frame(chunk_label="import_data", chunk_step_major=1, chunk_step_minor=0)
print(script_df)
```

```
##   chunk_label chunk_step_major chunk_step_minor
## 1 import_data                1                0
```

## Step `1`: import data

```r
glb_entity_df <- myimport_data(
    url="https://courses.edx.org/c4x/MITx/15.071x_2/asset/PollingData.csv", 
    comment="glb_entity_df", force_header=TRUE,
    print_diagn=(glb_is_separate_newent_dataset | 
                !glb_split_entity_newent_datasets))
```

```
## [1] "Reading file ./data/PollingData.csv..."
## [1] "dimensions of data in ./data/PollingData.csv: 145 rows x 7 cols"
```

```r
# Sync w/ Recitation
print(table(glb_entity_df$Year))
```

```
## 
## 2004 2008 2012 
##   50   50   45
```

```r
summary(glb_entity_df)
```

```
##     State                Year        Rasmussen          SurveyUSA       
##  Length:145         Min.   :2004   Min.   :-41.0000   Min.   :-33.0000  
##  Class :character   1st Qu.:2004   1st Qu.: -8.0000   1st Qu.:-11.7500  
##  Mode  :character   Median :2008   Median :  1.0000   Median : -2.0000  
##                     Mean   :2008   Mean   :  0.0404   Mean   : -0.8243  
##                     3rd Qu.:2012   3rd Qu.:  8.5000   3rd Qu.:  8.0000  
##                     Max.   :2012   Max.   : 39.0000   Max.   : 30.0000  
##                                    NA's   :46         NA's   :71        
##    DiffCount           PropR          Republican    
##  Min.   :-19.000   Min.   :0.0000   Min.   :0.0000  
##  1st Qu.: -6.000   1st Qu.:0.0000   1st Qu.:0.0000  
##  Median :  1.000   Median :0.6250   Median :1.0000  
##  Mean   : -1.269   Mean   :0.5259   Mean   :0.5103  
##  3rd Qu.:  4.000   3rd Qu.:1.0000   3rd Qu.:1.0000  
##  Max.   : 11.000   Max.   :1.0000   Max.   :1.0000  
## 
```

```r
if (glb_is_separate_newent_dataset) {
    glb_newent_df <- myimport_data(
        url="<newdt_url>", 
        comment="glb_newent_df", force_header=TRUE, print_diagn=TRUE)
} else {
    if (!glb_split_entity_newent_datasets) {
        stop("Not implemented yet") 
        glb_newent_df <- glb_entity_df[sample(1:nrow(glb_entity_df),
                                          max(2, nrow(glb_entity_df) / 1000)),]                    
    } else      if (glb_split_newdata_method == "condition") {
            glb_newent_df <- do.call("subset", 
                list(glb_entity_df, parse(text=glb_split_newdata_condition)))
            glb_entity_df <- do.call("subset", 
                list(glb_entity_df, parse(text=paste0("!(", 
                                                      glb_split_newdata_condition,
                                                      ")"))))
        } else if (glb_split_newdata_method == "sample") {
                require(caTools)
                
                set.seed(glb_split_sample.seed)
                split <- sample.split(glb_entity_df[, glb_predct_var], 
                                      SplitRatio=(1-glb_split_newdata_size))
                glb_newent_df <- glb_entity_df[!split, ] 
                glb_entity_df <- glb_entity_df[split ,]
        } else stop("glb_split_newdata_method should be %in% c('condition', 'sample')")   

    comment(glb_newent_df) <- "glb_newent_df"
    myprint_df(glb_newent_df)
    str(glb_newent_df)

    if (glb_split_entity_newent_datasets) {
        myprint_df(glb_entity_df)
        str(glb_entity_df)        
    }
}         
```

```
##          State Year Rasmussen SurveyUSA DiffCount     PropR Republican
## 7      Arizona 2012         8        NA         4 0.8333333          1
## 10    Arkansas 2012        NA        NA         2 1.0000000          1
## 13  California 2012        NA       -14        -6 0.0000000          0
## 16    Colorado 2012         3        -2        -5 0.3076923          0
## 19 Connecticut 2012        -7       -13        -8 0.0000000          0
## 24     Florida 2012         2         0         6 0.6666667          0
##            State Year Rasmussen SurveyUSA DiffCount      PropR Republican
## 27       Georgia 2012        NA         8         4 1.00000000          1
## 63      Michigan 2012        -5        NA       -10 0.08333333          0
## 105     Oklahoma 2012        NA        NA         1 1.00000000          1
## 120 South Dakota 2012        NA        NA         1 1.00000000          1
## 123    Tennessee 2012        NA        NA         1 1.00000000          1
## 143    Wisconsin 2012         0        NA        -8 0.00000000          0
##             State Year Rasmussen SurveyUSA DiffCount     PropR Republican
## 126         Texas 2012        NA        NA         4 1.0000000          1
## 129          Utah 2012        NA        NA         1 1.0000000          1
## 134      Virginia 2012         2        NA        -4 0.3333333          0
## 137    Washington 2012       -13       -14        -8 0.0000000          0
## 140 West Virginia 2012        NA        NA         1 1.0000000          1
## 143     Wisconsin 2012         0        NA        -8 0.0000000          0
## 'data.frame':	45 obs. of  7 variables:
##  $ State     : chr  "Arizona" "Arkansas" "California" "Colorado" ...
##  $ Year      : int  2012 2012 2012 2012 2012 2012 2012 2012 2012 2012 ...
##  $ Rasmussen : int  8 NA NA 3 -7 2 NA NA NA NA ...
##  $ SurveyUSA : int  NA NA -14 -2 -13 0 8 NA NA NA ...
##  $ DiffCount : int  4 2 -6 -5 -8 6 4 -2 1 -5 ...
##  $ PropR     : num  0.833 1 0 0.308 0 ...
##  $ Republican: int  1 1 0 0 0 0 1 0 1 0 ...
##  - attr(*, "comment")= chr "glb_newent_df"
##     State Year Rasmussen SurveyUSA DiffCount PropR Republican
## 1 Alabama 2004        11        18         5     1          1
## 2 Alabama 2008        21        25         5     1          1
## 3  Alaska 2004        NA        NA         1     1          1
## 4  Alaska 2008        16        NA         6     1          1
## 5 Arizona 2004         5        15         8     1          1
## 6 Arizona 2008         5        NA         9     1          1
##         State Year Rasmussen SurveyUSA DiffCount     PropR Republican
## 4      Alaska 2008        16        NA         6 1.0000000          1
## 20   Delaware 2004        NA        NA        -2 0.0000000          0
## 46   Kentucky 2004        NA        21         3 1.0000000          1
## 73    Montana 2004        NA        NA         3 1.0000000          1
## 104  Oklahoma 2008        31        24         2 1.0000000          1
## 141 Wisconsin 2004        -1        NA         1 0.5333333          0
##             State Year Rasmussen SurveyUSA DiffCount     PropR Republican
## 138 West Virginia 2004         6        NA         6 1.0000000          1
## 139 West Virginia 2008         9        NA        11 1.0000000          1
## 141     Wisconsin 2004        -1        NA         1 0.5333333          0
## 142     Wisconsin 2008        -7       -16       -12 0.0000000          0
## 144       Wyoming 2004        NA        NA         1 1.0000000          1
## 145       Wyoming 2008        19        21         3 1.0000000          1
## 'data.frame':	100 obs. of  7 variables:
##  $ State     : chr  "Alabama" "Alabama" "Alaska" "Alaska" ...
##  $ Year      : int  2004 2008 2004 2008 2004 2008 2004 2008 2004 2008 ...
##  $ Rasmussen : int  11 21 NA 16 5 5 7 10 -11 -27 ...
##  $ SurveyUSA : int  18 25 NA NA 15 NA 5 NA -11 -24 ...
##  $ DiffCount : int  5 5 1 6 8 9 8 5 -8 -5 ...
##  $ PropR     : num  1 1 1 1 1 1 1 1 0 0 ...
##  $ Republican: int  1 1 1 1 1 1 1 1 0 0 ...
```

```r
script_df <- rbind(script_df,
                   data.frame(chunk_label="cleanse_data", 
                              chunk_step_major=max(script_df$chunk_step_major)+1, 
                              chunk_step_minor=0))
print(script_df)
```

```
##    chunk_label chunk_step_major chunk_step_minor
## 1  import_data                1                0
## 2 cleanse_data                2                0
```

## Step `2`: cleanse data

```r
script_df <- rbind(script_df, 
                   data.frame(chunk_label="inspect_explore_data", 
                              chunk_step_major=max(script_df$chunk_step_major), 
                              chunk_step_minor=1))
print(script_df)
```

```
##            chunk_label chunk_step_major chunk_step_minor
## 1          import_data                1                0
## 2         cleanse_data                2                0
## 3 inspect_explore_data                2                1
```

### Step `2`.`1`: inspect/explore data

```r
#print(str(glb_entity_df))
#View(glb_entity_df)

# List info gathered for various columns
# <col_name>:   <description>; <notes>
# PropR: proportion Republican: proportion of all polls leading up to the election that predicted a Republican winner.

# Create new features that help diagnostics
#   Convert factors to dummy variables
#   Build splines   require(splines); bsBasis <- bs(training$age, df=3)

add_new_diag_feats <- function(obs_df, obs_twin_df) {
    require(plyr)
    
    obs_df <- mutate(obs_df,
#         <col_name>.NA=is.na(<col_name>)

#         <col_name>.fctr=factor(<col_name>, 
#                     as.factor(union(obs_df$<col_name>, obs_twin_df$<col_name>))) 
#         <col_name>.fctr=relevel(factor(<col_name>, 
#                     as.factor(union(obs_df$<col_name>, obs_twin_df$<col_name>))),
#                                   "<max_n_val>") 

          # This doesn't work - use sapply instead
#         <col_name>.fctr_num=grep(<col_name>, levels(<col_name>.fctr)), 
#         
#         Date.my=as.Date(strptime(Date, "%m/%d/%y %H:%M")),
#         Year=year(Date.my),
#         Month=months(Date.my),
#         Weekday=weekdays(Date.my)

#         <col_name>.log=log(<col.name>)        
        .rnorm=rnorm(1)
                        )

    # If levels of a factor are different across obs_df & glb_newent_df; predict.glm fails  
    # Transformations not handled by mutate
#     obs_df$<col_name>.fctr.num <- sapply(1:nrow(obs_df), 
#         function(row_ix) grep(obs_df[row_ix, "<col_name>"],
#                               levels(obs_df[row_ix, "<col_name>.fctr"])))
    
    print(summary(obs_df))
    print(sapply(names(obs_df), function(col) sum(is.na(obs_df[, col]))))
    return(obs_df)
}

glb_entity_df <- add_new_diag_feats(glb_entity_df, glb_newent_df)
```

```
## Loading required package: plyr
```

```
##     State                Year        Rasmussen          SurveyUSA       
##  Length:100         Min.   :2004   Min.   :-41.0000   Min.   :-33.0000  
##  Class :character   1st Qu.:2004   1st Qu.:-10.0000   1st Qu.:-11.0000  
##  Mode  :character   Median :2006   Median :  1.0000   Median : -1.0000  
##                     Mean   :2006   Mean   :  0.2467   Mean   :  0.1579  
##                     3rd Qu.:2008   3rd Qu.:  9.0000   3rd Qu.: 12.0000  
##                     Max.   :2008   Max.   : 39.0000   Max.   : 30.0000  
##                                    NA's   :23         NA's   :43        
##    DiffCount          PropR          Republican       .rnorm      
##  Min.   :-19.00   Min.   :0.0000   Min.   :0.00   Min.   :0.6301  
##  1st Qu.: -5.25   1st Qu.:0.0000   1st Qu.:0.00   1st Qu.:0.6301  
##  Median :  1.00   Median :0.6458   Median :1.00   Median :0.6301  
##  Mean   : -0.88   Mean   :0.5383   Mean   :0.53   Mean   :0.6301  
##  3rd Qu.:  4.00   3rd Qu.:1.0000   3rd Qu.:1.00   3rd Qu.:0.6301  
##  Max.   : 11.00   Max.   :1.0000   Max.   :1.00   Max.   :0.6301  
##                                                                   
##      State       Year  Rasmussen  SurveyUSA  DiffCount      PropR 
##          0          0         23         43          0          0 
## Republican     .rnorm 
##          0          0
```

```r
glb_newent_df <- add_new_diag_feats(glb_newent_df, glb_entity_df)
```

```
##     State                Year        Rasmussen          SurveyUSA      
##  Length:45          Min.   :2012   Min.   :-19.0000   Min.   :-29.000  
##  Class :character   1st Qu.:2012   1st Qu.: -5.0000   1st Qu.:-13.000  
##  Mode  :character   Median :2012   Median :  0.0000   Median : -4.000  
##                     Mean   :2012   Mean   : -0.6818   Mean   : -4.118  
##                     3rd Qu.:2012   3rd Qu.:  5.2500   3rd Qu.:  5.000  
##                     Max.   :2012   Max.   : 14.0000   Max.   : 14.000  
##                                    NA's   :23         NA's   :28       
##    DiffCount           PropR          Republican         .rnorm       
##  Min.   :-16.000   Min.   :0.0000   Min.   :0.0000   Min.   :-0.2762  
##  1st Qu.: -6.000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:-0.2762  
##  Median : -2.000   Median :0.4000   Median :0.0000   Median :-0.2762  
##  Mean   : -2.133   Mean   :0.4985   Mean   :0.4667   Mean   :-0.2762  
##  3rd Qu.:  2.000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:-0.2762  
##  Max.   :  8.000   Max.   :1.0000   Max.   :1.0000   Max.   :-0.2762  
##                                                                       
##      State       Year  Rasmussen  SurveyUSA  DiffCount      PropR 
##          0          0         23         28          0          0 
## Republican     .rnorm 
##          0          0
```

```r
#pairs(subset(glb_entity_df, select=-c(col_symbol)))

#   Histogram of predictor in glb_entity_df & glb_newent_df
# Check for glb_newent_df & glb_entity_df features range mismatches

# Other diagnostics:
# print(subset(glb_entity_df, <col1_name> == max(glb_entity_df$<col1_name>, na.rm=TRUE) & 
#                         <col2_name> <= mean(glb_entity_df$<col1_name>, na.rm=TRUE)))

# print(glb_entity_df[which.max(glb_entity_df$<col_name>),])

# print(<col_name>_freq_glb_entity_df <- mycreate_tbl_df(glb_entity_df, "<col_name>"))
# print(which.min(table(glb_entity_df$<col_name>)))
# print(which.max(table(glb_entity_df$<col_name>)))
# print(which.max(table(glb_entity_df$<col1_name>, glb_entity_df$<col2_name>)[, 2]))
# print(table(glb_entity_df$<col1_name>, glb_entity_df$<col2_name>))
# print(table(is.na(glb_entity_df$<col1_name>), glb_entity_df$<col2_name>))
# print(mycreate_xtab(glb_entity_df, c(<col1_name>, <col2_name>)))
# print(<col1_name>_<col2_name>_xtab_glb_entity_df <- 
#   mycreate_xtab(glb_entity_df, c("<col1_name>", "<col2_name>")))
# <col1_name>_<col2_name>_xtab_glb_entity_df[is.na(<col1_name>_<col2_name>_xtab_glb_entity_df)] <- 0
# print(<col1_name>_<col2_name>_xtab_glb_entity_df <- 
#   mutate(<col1_name>_<col2_name>_xtab_glb_entity_df, 
#             <col3_name>=(<col1_name> * 1.0) / (<col1_name> + <col2_name>))) 

# print(<col2_name>_min_entity_arr <- 
#    sort(tapply(glb_entity_df$<col1_name>, glb_entity_df$<col2_name>, min, na.rm=TRUE)))
# print(<col1_name>_na_by_<col2_name>_arr <- 
#    sort(tapply(glb_entity_df$<col1_name>.NA, glb_entity_df$<col2_name>, mean, na.rm=TRUE)))

# Other plots:
print(myplot_histogram(glb_entity_df, glb_predct_var))
```

```
## stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust this.
```

![](RCP_files/figure-html/inspect_explore_data_1-1.png) 

```r
# print(myplot_box(df=glb_entity_df, ycol_names="<col1_name>"))
# print(myplot_box(df=glb_entity_df, ycol_names="<col1_name>", xcol_name="<col2_name>"))
# print(myplot_line(subset(glb_entity_df, Symbol %in% c("KO", "PG")), 
#                   "Date.my", "StockPrice", facet_row_colnames="Symbol") + 
#     geom_vline(xintercept=as.numeric(as.Date("2003-03-01"))) +
#     geom_vline(xintercept=as.numeric(as.Date("1983-01-01")))        
#         )
# print(myplot_scatter(glb_entity_df, "<col1_name>", "<col2_name>", smooth=TRUE))
# print(myplot_scatter(glb_entity_df, "<col1_name>", "<col2_name>", colorcol_name="<Pred.fctr>"))

script_df <- rbind(script_df, 
    data.frame(chunk_label="manage_missing_data", 
        chunk_step_major=max(script_df$chunk_step_major), 
        chunk_step_minor=script_df[nrow(script_df), "chunk_step_minor"]+1))
print(script_df)
```

```
##            chunk_label chunk_step_major chunk_step_minor
## 1          import_data                1                0
## 2         cleanse_data                2                0
## 3 inspect_explore_data                2                1
## 4  manage_missing_data                2                2
```

### Step `2`.`2`: manage missing data

```r
# print(sapply(names(glb_entity_df), function(col) sum(is.na(glb_entity_df[, col]))))
# print(sapply(names(glb_newent_df), function(col) sum(is.na(glb_newent_df[, col]))))
# glb_entity_df <- na.omit(glb_entity_df)
# glb_newent_df <- na.omit(glb_newent_df)
# df[is.na(df)] <- 0

# Impute missing data
if (!glb_is_separate_newent_dataset) {
    # Combine entity & newent
    union_df <- orderBy(~ State + Year, 
                       rbind(mutate(glb_entity_df, .src = "entity"),
                             mutate(glb_newent_df, .src = "newent")))
    union_imputed_df <- union_df[, setdiff(setdiff(names(glb_entity_df), 
                                                   glb_predct_var), 
                                           glb_exclude_vars_as_features)]
    print(summary(union_imputed_df))

    require(mice)
    set.seed(glb_mice_complete.seed)
    union_imputed_df <- complete(mice(union_imputed_df))
    print(summary(union_imputed_df))
    #print(all.equal(union_imputed_df, imputed))

    #union_df[, names(union_imputed_df)] <- union_imputed_df[, names(union_imputed_df)]
    # sync w/ Lecture
    tmp_imputed_df <- myimport_data(
    url="https://courses.edx.org/c4x/MITx/15.071x_2/asset/PollingData_Imputed.csv", 
    comment="tmp_imputed_df", force_header=TRUE, print_diagn=FALSE)
    print(summary(tmp_imputed_df))
   
    union_df <- merge(
        union_df[, union(c(glb_id_vars, ".src"), 
                         setdiff(names(union_df), names(tmp_imputed_df)))], 
                      tmp_imputed_df)
    #print(all.equal(union_df[, names(union_imputed_df)], imputed))
    print(summary(union_df))
    
    # Partition again
    glb_entity_df <- subset(union_df, .src == "entity", select=-.src)
    comment(glb_entity_df) <- "glb_entity_df"
    glb_newent_df <- subset(union_df, .src == "newent", select=-.src)
    comment(glb_newent_df) <- "glb_newent_df"
    
    # Generate summaries
    print(summary(glb_entity_df))
    print(sapply(names(glb_entity_df), function(col) sum(is.na(glb_entity_df[, col]))))
    print(summary(glb_newent_df))
    print(sapply(names(glb_newent_df), function(col) sum(is.na(glb_newent_df[, col]))))

} else stop("Not implemented yet")
```

```
##    Rasmussen          SurveyUSA          DiffCount           PropR       
##  Min.   :-41.0000   Min.   :-33.0000   Min.   :-19.000   Min.   :0.0000  
##  1st Qu.: -8.0000   1st Qu.:-11.7500   1st Qu.: -6.000   1st Qu.:0.0000  
##  Median :  1.0000   Median : -2.0000   Median :  1.000   Median :0.6250  
##  Mean   :  0.0404   Mean   : -0.8243   Mean   : -1.269   Mean   :0.5259  
##  3rd Qu.:  8.5000   3rd Qu.:  8.0000   3rd Qu.:  4.000   3rd Qu.:1.0000  
##  Max.   : 39.0000   Max.   : 30.0000   Max.   : 11.000   Max.   :1.0000  
##  NA's   :46         NA's   :71
```

```
## Loading required package: mice
## Loading required package: Rcpp
## Loading required package: lattice
## mice 2.22 2014-06-10
```

```
## 
##  iter imp variable
##   1   1  Rasmussen  SurveyUSA
##   1   2  Rasmussen  SurveyUSA
##   1   3  Rasmussen  SurveyUSA
##   1   4  Rasmussen  SurveyUSA
##   1   5  Rasmussen  SurveyUSA
##   2   1  Rasmussen  SurveyUSA
##   2   2  Rasmussen  SurveyUSA
##   2   3  Rasmussen  SurveyUSA
##   2   4  Rasmussen  SurveyUSA
##   2   5  Rasmussen  SurveyUSA
##   3   1  Rasmussen  SurveyUSA
##   3   2  Rasmussen  SurveyUSA
##   3   3  Rasmussen  SurveyUSA
##   3   4  Rasmussen  SurveyUSA
##   3   5  Rasmussen  SurveyUSA
##   4   1  Rasmussen  SurveyUSA
##   4   2  Rasmussen  SurveyUSA
##   4   3  Rasmussen  SurveyUSA
##   4   4  Rasmussen  SurveyUSA
##   4   5  Rasmussen  SurveyUSA
##   5   1  Rasmussen  SurveyUSA
##   5   2  Rasmussen  SurveyUSA
##   5   3  Rasmussen  SurveyUSA
##   5   4  Rasmussen  SurveyUSA
##   5   5  Rasmussen  SurveyUSA
##    Rasmussen         SurveyUSA         DiffCount           PropR       
##  Min.   :-41.000   Min.   :-33.000   Min.   :-19.000   Min.   :0.0000  
##  1st Qu.: -8.000   1st Qu.:-11.000   1st Qu.: -6.000   1st Qu.:0.0000  
##  Median :  3.000   Median :  1.000   Median :  1.000   Median :0.6250  
##  Mean   :  2.062   Mean   :  1.138   Mean   : -1.269   Mean   :0.5259  
##  3rd Qu.: 12.000   3rd Qu.: 16.000   3rd Qu.:  4.000   3rd Qu.:1.0000  
##  Max.   : 39.000   Max.   : 30.000   Max.   : 11.000   Max.   :1.0000  
## [1] "Reading file ./data/PollingData_Imputed.csv..."
## [1] "dimensions of data in ./data/PollingData_Imputed.csv: 145 rows x 7 cols"
##     State                Year        Rasmussen         SurveyUSA      
##  Length:145         Min.   :2004   Min.   :-41.000   Min.   :-33.000  
##  Class :character   1st Qu.:2004   1st Qu.:-10.000   1st Qu.:-11.000  
##  Mode  :character   Median :2008   Median :  3.000   Median :  1.000  
##                     Mean   :2008   Mean   :  2.048   Mean   :  1.359  
##                     3rd Qu.:2012   3rd Qu.: 12.000   3rd Qu.: 16.000  
##                     Max.   :2012   Max.   : 39.000   Max.   : 30.000  
##    DiffCount           PropR          Republican    
##  Min.   :-19.000   Min.   :0.0000   Min.   :0.0000  
##  1st Qu.: -6.000   1st Qu.:0.0000   1st Qu.:0.0000  
##  Median :  1.000   Median :0.6250   Median :1.0000  
##  Mean   : -1.269   Mean   :0.5259   Mean   :0.5103  
##  3rd Qu.:  4.000   3rd Qu.:1.0000   3rd Qu.:1.0000  
##  Max.   : 11.000   Max.   :1.0000   Max.   :1.0000  
##     State                Year          .src               .rnorm       
##  Length:145         Min.   :2004   Length:145         Min.   :-0.2762  
##  Class :character   1st Qu.:2004   Class :character   1st Qu.:-0.2762  
##  Mode  :character   Median :2008   Mode  :character   Median : 0.6301  
##                     Mean   :2008                      Mean   : 0.3488  
##                     3rd Qu.:2012                      3rd Qu.: 0.6301  
##                     Max.   :2012                      Max.   : 0.6301  
##    Rasmussen         SurveyUSA         DiffCount           PropR       
##  Min.   :-41.000   Min.   :-33.000   Min.   :-19.000   Min.   :0.0000  
##  1st Qu.:-10.000   1st Qu.:-11.000   1st Qu.: -6.000   1st Qu.:0.0000  
##  Median :  3.000   Median :  1.000   Median :  1.000   Median :0.6250  
##  Mean   :  2.048   Mean   :  1.359   Mean   : -1.269   Mean   :0.5259  
##  3rd Qu.: 12.000   3rd Qu.: 16.000   3rd Qu.:  4.000   3rd Qu.:1.0000  
##  Max.   : 39.000   Max.   : 30.000   Max.   : 11.000   Max.   :1.0000  
##    Republican    
##  Min.   :0.0000  
##  1st Qu.:0.0000  
##  Median :1.0000  
##  Mean   :0.5103  
##  3rd Qu.:1.0000  
##  Max.   :1.0000  
##     State                Year          .rnorm         Rasmussen     
##  Length:100         Min.   :2004   Min.   :0.6301   Min.   :-41.00  
##  Class :character   1st Qu.:2004   1st Qu.:0.6301   1st Qu.:-10.00  
##  Mode  :character   Median :2006   Median :0.6301   Median :  3.50  
##                     Mean   :2006   Mean   :0.6301   Mean   :  2.02  
##                     3rd Qu.:2008   3rd Qu.:0.6301   3rd Qu.: 12.00  
##                     Max.   :2008   Max.   :0.6301   Max.   : 39.00  
##    SurveyUSA        DiffCount          PropR          Republican  
##  Min.   :-33.00   Min.   :-19.00   Min.   :0.0000   Min.   :0.00  
##  1st Qu.:-11.00   1st Qu.: -5.25   1st Qu.:0.0000   1st Qu.:0.00  
##  Median :  1.50   Median :  1.00   Median :0.6458   Median :1.00  
##  Mean   :  1.32   Mean   : -0.88   Mean   :0.5383   Mean   :0.53  
##  3rd Qu.: 16.50   3rd Qu.:  4.00   3rd Qu.:1.0000   3rd Qu.:1.00  
##  Max.   : 30.00   Max.   : 11.00   Max.   :1.0000   Max.   :1.00  
##      State       Year     .rnorm  Rasmussen  SurveyUSA  DiffCount 
##          0          0          0          0          0          0 
##      PropR Republican 
##          0          0 
##     State                Year          .rnorm          Rasmussen      
##  Length:45          Min.   :2012   Min.   :-0.2762   Min.   :-28.000  
##  Class :character   1st Qu.:2012   1st Qu.:-0.2762   1st Qu.:-11.000  
##  Mode  :character   Median :2012   Median :-0.2762   Median :  2.000  
##                     Mean   :2012   Mean   :-0.2762   Mean   :  2.111  
##                     3rd Qu.:2012   3rd Qu.:-0.2762   3rd Qu.: 11.000  
##                     Max.   :2012   Max.   :-0.2762   Max.   : 34.000  
##    SurveyUSA         DiffCount           PropR          Republican    
##  Min.   :-30.000   Min.   :-16.000   Min.   :0.0000   Min.   :0.0000  
##  1st Qu.:-11.000   1st Qu.: -6.000   1st Qu.:0.0000   1st Qu.:0.0000  
##  Median :  0.000   Median : -2.000   Median :0.4000   Median :0.0000  
##  Mean   :  1.444   Mean   : -2.133   Mean   :0.4985   Mean   :0.4667  
##  3rd Qu.: 16.000   3rd Qu.:  2.000   3rd Qu.:1.0000   3rd Qu.:1.0000  
##  Max.   : 25.000   Max.   :  8.000   Max.   :1.0000   Max.   :1.0000  
##      State       Year     .rnorm  Rasmussen  SurveyUSA  DiffCount 
##          0          0          0          0          0          0 
##      PropR Republican 
##          0          0
```

```r
print(mycreate_xtab(df=glb_entity_df, xtab_col_names="Republican"))
```

```
## Loading required package: reshape2
```

```
##   Republican _n.length
## 1          0        47
## 2          1        53
```

```r
print(table(sign(glb_entity_df$Rasmussen)))
```

```
## 
## -1  0  1 
## 42  2 56
```

```r
print(xtabs(~ Republican + sign(Rasmussen), glb_entity_df))
```

```
##           sign(Rasmussen)
## Republican -1  0  1
##          0 42  1  4
##          1  0  1 52
```

```r
script_df <- rbind(script_df, 
    data.frame(chunk_label="encode_retype_data", 
        chunk_step_major=max(script_df$chunk_step_major), 
        chunk_step_minor=script_df[nrow(script_df), "chunk_step_minor"]+1))
print(script_df)
```

```
##            chunk_label chunk_step_major chunk_step_minor
## 1          import_data                1                0
## 2         cleanse_data                2                0
## 3 inspect_explore_data                2                1
## 4  manage_missing_data                2                2
## 5   encode_retype_data                2                3
```

### Step `2`.`3`: encode/retype data

```r
# map_<col_name>_df <- myimport_data(
#     url="<map_url>", 
#     comment="map_<col_name>_df", print_diagn=TRUE)
# map_<col_name>_df <- read.csv(paste0(getwd(), "/data/<file_name>.csv"), strip.white=TRUE)

# glb_entity_df <- mymap_codes(glb_entity_df, "<from_col_name>", "<to_col_name>", 
#     map_<to_col_name>_df, map_join_col_name="<map_join_col_name>", 
#                           map_tgt_col_name="<to_col_name>")
# glb_newent_df <- mymap_codes(glb_newent_df, "<from_col_name>", "<to_col_name>", 
#     map_<to_col_name>_df, map_join_col_name="<map_join_col_name>", 
#                           map_tgt_col_name="<to_col_name>")
    					
# glb_entity_df$<col_name>.fctr <- factor(glb_entity_df$<col_name>, 
#                     as.factor(union(glb_entity_df$<col_name>, glb_newent_df$<col_name>)))
# glb_newent_df$<col_name>.fctr <- factor(glb_newent_df$<col_name>, 
#                     as.factor(union(glb_entity_df$<col_name>, glb_newent_df$<col_name>)))

script_df <- rbind(script_df, 
                   data.frame(chunk_label="extract_features", 
                              chunk_step_major=max(script_df$chunk_step_major)+1, 
                              chunk_step_minor=0))
print(script_df)
```

```
##            chunk_label chunk_step_major chunk_step_minor
## 1          import_data                1                0
## 2         cleanse_data                2                0
## 3 inspect_explore_data                2                1
## 4  manage_missing_data                2                2
## 5   encode_retype_data                2                3
## 6     extract_features                3                0
```

## Step `3`: extract features

```r
# Create new features that help prediction
# <col_name>.lag.2 <- lag(zoo(glb_entity_df$<col_name>), -2, na.pad=TRUE)
# glb_entity_df[, "<col_name>.lag.2"] <- coredata(<col_name>.lag.2)
# <col_name>.lag.2 <- lag(zoo(glb_newent_df$<col_name>), -2, na.pad=TRUE)
# glb_newent_df[, "<col_name>.lag.2"] <- coredata(<col_name>.lag.2)
# 
# glb_newent_df[1, "<col_name>.lag.2"] <- glb_entity_df[nrow(glb_entity_df) - 1, 
#                                                    "<col_name>"]
# glb_newent_df[2, "<col_name>.lag.2"] <- glb_entity_df[nrow(glb_entity_df), 
#                                                    "<col_name>"]
                                                   
# glb_entity_df <- mutate(glb_entity_df,
#     <new_col_name>=
#                     )

# glb_newent_df <- mutate(glb_newent_df,
#     <new_col_name>=
#                     )

# print(summary(glb_entity_df))
# print(summary(glb_newent_df))

# print(sapply(names(glb_entity_df), function(col) sum(is.na(glb_entity_df[, col]))))
# print(sapply(names(glb_newent_df), function(col) sum(is.na(glb_newent_df[, col]))))

# print(myplot_scatter(glb_entity_df, "<col1_name>", "<col2_name>", smooth=TRUE))

script_df <- rbind(script_df, 
                   data.frame(chunk_label="select_features", 
                              chunk_step_major=max(script_df$chunk_step_major)+1, 
                              chunk_step_minor=0))
print(script_df)
```

```
##            chunk_label chunk_step_major chunk_step_minor
## 1          import_data                1                0
## 2         cleanse_data                2                0
## 3 inspect_explore_data                2                1
## 4  manage_missing_data                2                2
## 5   encode_retype_data                2                3
## 6     extract_features                3                0
## 7      select_features                4                0
```

## Step `4`: select features

```r
print(glb_feats_df <- 
    myselect_features(glb_entity_df, glb_exclude_vars_as_features, glb_predct_var))
```

```
##                  id     cor.y cor.y.abs
## PropR         PropR 0.9484204 0.9484204
## SurveyUSA SurveyUSA 0.8101645 0.8101645
## DiffCount DiffCount 0.8092777 0.8092777
## Rasmussen Rasmussen 0.7929252 0.7929252
```

```r
script_df <- rbind(script_df, 
    data.frame(chunk_label="remove_correlated_features", 
        chunk_step_major=max(script_df$chunk_step_major),
        chunk_step_minor=script_df[nrow(script_df), "chunk_step_minor"]+1))        
print(script_df)
```

```
##                  chunk_label chunk_step_major chunk_step_minor
## 1                import_data                1                0
## 2               cleanse_data                2                0
## 3       inspect_explore_data                2                1
## 4        manage_missing_data                2                2
## 5         encode_retype_data                2                3
## 6           extract_features                3                0
## 7            select_features                4                0
## 8 remove_correlated_features                4                1
```

### Step `4`.`1`: remove correlated features

```r
print(glb_feats_df <- orderBy(~-cor.y, merge(glb_feats_df, 
          mydelete_cor_features(glb_feats_df, glb_entity_df, glb_predct_var, 
                                glb_exclude_vars_as_features), 
          all.x=TRUE)))
```

```
##               PropR SurveyUSA DiffCount Rasmussen
## PropR     1.0000000 0.8616478 0.8273785 0.8431180
## SurveyUSA 0.8616478 1.0000000 0.5222585 0.9365837
## DiffCount 0.8273785 0.5222585 1.0000000 0.5109169
## Rasmussen 0.8431180 0.9365837 0.5109169 1.0000000
##               PropR SurveyUSA DiffCount Rasmussen
## PropR     0.0000000 0.8616478 0.8273785 0.8431180
## SurveyUSA 0.8616478 0.0000000 0.5222585 0.9365837
## DiffCount 0.8273785 0.5222585 0.0000000 0.5109169
## Rasmussen 0.8431180 0.9365837 0.5109169 0.0000000
## [1] "cor(SurveyUSA, Rasmussen)=0.9366"
```

![](RCP_files/figure-html/remove_correlated_features-1.png) 

```
## [1] "cor(Republican, SurveyUSA)=0.8102"
## [1] "cor(Republican, Rasmussen)=0.7929"
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : pseudoinverse used at -0.005
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : neighborhood radius 1.005
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : reciprocal condition number 0
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : There are other near singularities as well. 1.01
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : pseudoinverse used at -0.005
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : neighborhood radius 1.005
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : reciprocal condition number 0
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : There are other near singularities as well. 1.01
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : pseudoinverse used at -0.005
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : neighborhood radius 1.005
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : reciprocal condition number 0
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : There are other near singularities as well. 1.01
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : pseudoinverse used at -0.005
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : neighborhood radius 1.005
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : reciprocal condition number 0
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : There are other near singularities as well. 1.01
```

```
## Warning in mydelete_cor_features(glb_feats_df, glb_entity_df,
## glb_predct_var, : Dropping Rasmussen as a feature
```

![](RCP_files/figure-html/remove_correlated_features-2.png) 

```
##                  id     cor.y cor.y.abs
## PropR         PropR 0.9484204 0.9484204
## SurveyUSA SurveyUSA 0.8101645 0.8101645
## DiffCount DiffCount 0.8092777 0.8092777
##               PropR SurveyUSA DiffCount
## PropR     1.0000000 0.8616478 0.8273785
## SurveyUSA 0.8616478 1.0000000 0.5222585
## DiffCount 0.8273785 0.5222585 1.0000000
##               PropR SurveyUSA DiffCount
## PropR     0.0000000 0.8616478 0.8273785
## SurveyUSA 0.8616478 0.0000000 0.5222585
## DiffCount 0.8273785 0.5222585 0.0000000
## [1] "cor(PropR, SurveyUSA)=0.8616"
```

![](RCP_files/figure-html/remove_correlated_features-3.png) 

```
## [1] "cor(Republican, PropR)=0.9484"
## [1] "cor(Republican, SurveyUSA)=0.8102"
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : pseudoinverse used at -0.005
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : neighborhood radius 1.005
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : reciprocal condition number 0
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : There are other near singularities as well. 1.01
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : pseudoinverse used at -0.005
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : neighborhood radius 1.005
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : reciprocal condition number 0
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : There are other near singularities as well. 1.01
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : pseudoinverse used at -0.005
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : neighborhood radius 1.005
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : reciprocal condition number 0
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : There are other near singularities as well. 1.01
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : pseudoinverse used at -0.005
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : neighborhood radius 1.005
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : reciprocal condition number 0
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : There are other near singularities as well. 1.01
```

```
## Warning in mydelete_cor_features(glb_feats_df, glb_entity_df,
## glb_predct_var, : Dropping SurveyUSA as a feature
```

![](RCP_files/figure-html/remove_correlated_features-4.png) 

```
##                  id     cor.y cor.y.abs
## PropR         PropR 0.9484204 0.9484204
## DiffCount DiffCount 0.8092777 0.8092777
##               PropR DiffCount
## PropR     1.0000000 0.8273785
## DiffCount 0.8273785 1.0000000
##               PropR DiffCount
## PropR     0.0000000 0.8273785
## DiffCount 0.8273785 0.0000000
## [1] "cor(PropR, DiffCount)=0.8274"
```

![](RCP_files/figure-html/remove_correlated_features-5.png) 

```
## [1] "cor(Republican, PropR)=0.9484"
## [1] "cor(Republican, DiffCount)=0.8093"
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : pseudoinverse used at -0.005
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : neighborhood radius 1.005
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : reciprocal condition number 0
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : There are other near singularities as well. 1.01
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : pseudoinverse used at -0.005
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : neighborhood radius 1.005
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : reciprocal condition number 0
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : There are other near singularities as well. 1.01
```

```
## geom_smooth: method="auto" and size of largest group is <1000, so using loess. Use 'method = x' to change the smoothing method.
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : pseudoinverse used at -0.005
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : neighborhood radius 1.005
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : reciprocal condition number 0
```

```
## Warning in simpleLoess(y, x, w, span, degree, parametric, drop.square,
## normalize, : There are other near singularities as well. 1.01
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : pseudoinverse used at -0.005
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : neighborhood radius 1.005
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : reciprocal condition number 0
```

```
## Warning in predLoess(y, x, newx, s, weights, pars$robust, pars$span,
## pars$degree, : There are other near singularities as well. 1.01
```

```
## Warning in mydelete_cor_features(glb_feats_df, glb_entity_df,
## glb_predct_var, : Dropping DiffCount as a feature
```

![](RCP_files/figure-html/remove_correlated_features-6.png) 

```
##          id     cor.y cor.y.abs
## PropR PropR 0.9484204 0.9484204
##          id     cor.y cor.y.abs cor.low
## 2     PropR 0.9484204 0.9484204       1
## 4 SurveyUSA 0.8101645 0.8101645      NA
## 1 DiffCount 0.8092777 0.8092777      NA
## 3 Rasmussen 0.7929252 0.7929252      NA
```

```r
script_df <- rbind(script_df, 
                   data.frame(chunk_label="run_models", 
                              chunk_step_major=max(script_df$chunk_step_major)+1, 
                              chunk_step_minor=0))
print(script_df)
```

```
##                  chunk_label chunk_step_major chunk_step_minor
## 1                import_data                1                0
## 2               cleanse_data                2                0
## 3       inspect_explore_data                2                1
## 4        manage_missing_data                2                2
## 5         encode_retype_data                2                3
## 6           extract_features                3                0
## 7            select_features                4                0
## 8 remove_correlated_features                4                1
## 9                 run_models                5                0
```

## Step `5`: run models

```r
max_cor_y_x_var <- subset(glb_feats_df, cor.low == 1)[1, "id"]

#   Regression:
if (glb_is_regression) {
    #   Linear:
    myrun_mdl_fn <- myrun_mdl_lm
}    

#   Classification:
if (glb_is_classification) {
    #   Logit Regression:
    myrun_mdl_fn <- myrun_mdl_glm
}    
    
# Add dummy model - random variable
#   Potential Enhancements:
#       For classifiers, it shd generate proba/outcomes that mimics the freq
#           distribution of glb_predct_var values; Right now it always generates
#           0 (most frequent ?)
ret_lst <- myrun_mdl_fn(indep_vars_vctr=".rnorm",
                        lcl_predct_var=glb_predct_var, 
                        lcl_predct_var_name=glb_predct_var_name,
                        fit_df=glb_entity_df, OOB_df=glb_newent_df)
```

```
## Loading required package: ROCR
## Loading required package: gplots
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```
## 
## Call:
## glm(formula = reformulate(indep_vars_vctr, response = lcl_predct_var), 
##     family = "binomial", data = fit_df)
## 
## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -1.229  -1.229   1.127   1.127   1.127  
## 
## Coefficients: (1 not defined because of singularities)
##             Estimate Std. Error z value Pr(>|z|)
## (Intercept)   0.1201     0.2004     0.6    0.549
## .rnorm            NA         NA      NA       NA
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 138.27  on 99  degrees of freedom
## Residual deviance: 138.27  on 99  degrees of freedom
## AIC: 140.27
## 
## Number of Fisher Scoring iterations: 3
## 
##    feats n.fit R.sq.fit R.sq.OOB Adj.R.sq.fit  SSE.fit SSE.OOB  AIC.fit
## 1 .rnorm   100       NA       NA           NA 401.4452      NA 140.2692
##   auc.fit auc.OOB
## 1     0.5     0.5
```

```r
glb_dmy_mdl <- glb_mdl

# Highest cor.y
ret_lst <- myrun_mdl_fn(indep_vars_vctr=max_cor_y_x_var,
                        lcl_predct_var=glb_predct_var, 
                        lcl_predct_var_name=glb_predct_var_name,
                        fit_df=glb_entity_df, OOB_df=glb_newent_df)
```

```
## 
## Call:
## glm(formula = reformulate(indep_vars_vctr, response = lcl_predct_var), 
##     family = "binomial", data = fit_df)
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -2.22880  -0.06541   0.10260   0.10260   1.37392  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -6.146      1.977  -3.108 0.001882 ** 
## PropR         11.390      3.153   3.613 0.000303 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 138.269  on 99  degrees of freedom
## Residual deviance:  15.772  on 98  degrees of freedom
## AIC: 19.772
## 
## Number of Fisher Scoring iterations: 8
## 
##    feats n.fit R.sq.fit R.sq.OOB Adj.R.sq.fit  SSE.fit SSE.OOB   AIC.fit
## 2  PropR   100       NA       NA           NA 272.8743      NA  19.77232
## 1 .rnorm   100       NA       NA           NA 401.4452      NA 140.26922
##     auc.fit   auc.OOB
## 2 0.9955841 0.9990079
## 1 0.5000000 0.5000000
```

```r
glb_sel_mdl <- glb_mdl

# Enhance Highest cor.y model with additions of interaction terms that were 
#   dropped due to high correlations
if (nrow(subset(glb_feats_df, is.na(cor.low))) > 0)
    ret_lst <- myrun_mdl_fn(indep_vars_vctr=c(max_cor_y_x_var, 
        paste(max_cor_y_x_var, 
              subset(glb_feats_df, is.na(cor.low))[, "id"], sep=":")),
                        glb_predct_var, glb_predct_var_name,
                            fit_df=glb_entity_df, OOB_df=glb_newent_df)    
```

```
## 
## Call:
## glm(formula = reformulate(indep_vars_vctr, response = lcl_predct_var), 
##     family = "binomial", data = fit_df)
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -1.68559  -0.11220   0.00006   0.00181   1.82804  
## 
## Coefficients:
##                 Estimate Std. Error z value Pr(>|z|)  
## (Intercept)      -5.0650     1.9664  -2.576   0.0100 *
## PropR            10.7834     5.2599   2.050   0.0404 *
## PropR:SurveyUSA   1.0233     1.1094   0.922   0.3563  
## PropR:DiffCount   0.4275     0.8457   0.505   0.6132  
## PropR:Rasmussen  -0.6211     0.6752  -0.920   0.3576  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 138.269  on 99  degrees of freedom
## Residual deviance:  13.666  on 95  degrees of freedom
## AIC: 23.666
## 
## Number of Fisher Scoring iterations: 11
## 
##                                                      feats n.fit R.sq.fit
## 3 PropR, PropR:SurveyUSA, PropR:DiffCount, PropR:Rasmussen   100       NA
## 2                                                    PropR   100       NA
## 1                                                   .rnorm   100       NA
##   R.sq.OOB Adj.R.sq.fit  SSE.fit SSE.OOB   AIC.fit   auc.fit   auc.OOB
## 3       NA           NA 164.9796      NA  23.66631 0.9967884 1.0000000
## 2       NA           NA 272.8743      NA  19.77232 0.9955841 0.9990079
## 1       NA           NA 401.4452      NA 140.26922 0.5000000 0.5000000
```

```r
# Low correlated X
ret_lst <- myrun_mdl_fn(indep_vars_vctr=subset(glb_feats_df, 
                                               cor.low == 1)[, "id"],
                        glb_predct_var, glb_predct_var_name,
                        fit_df=glb_entity_df, OOB_df=glb_newent_df)
```

```
## 
## Call:
## glm(formula = reformulate(indep_vars_vctr, response = lcl_predct_var), 
##     family = "binomial", data = fit_df)
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -2.22880  -0.06541   0.10260   0.10260   1.37392  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -6.146      1.977  -3.108 0.001882 ** 
## PropR         11.390      3.153   3.613 0.000303 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 138.269  on 99  degrees of freedom
## Residual deviance:  15.772  on 98  degrees of freedom
## AIC: 19.772
## 
## Number of Fisher Scoring iterations: 8
## 
##                                                      feats n.fit R.sq.fit
## 3 PropR, PropR:SurveyUSA, PropR:DiffCount, PropR:Rasmussen   100       NA
## 2                                                    PropR   100       NA
## 4                                                    PropR   100       NA
## 1                                                   .rnorm   100       NA
##   R.sq.OOB Adj.R.sq.fit  SSE.fit SSE.OOB   AIC.fit   auc.fit   auc.OOB
## 3       NA           NA 164.9796      NA  23.66631 0.9967884 1.0000000
## 2       NA           NA 272.8743      NA  19.77232 0.9955841 0.9990079
## 4       NA           NA 272.8743      NA  19.77232 0.9955841 0.9990079
## 1       NA           NA 401.4452      NA 140.26922 0.5000000 0.5000000
```

```r
# All X that is not user excluded
ret_lst <- myrun_mdl_fn(indep_vars_vctr=setdiff(setdiff(names(glb_entity_df),
                                                        glb_predct_var),
                                                glb_exclude_vars_as_features),
                        glb_predct_var, glb_predct_var_name,
                        fit_df=glb_entity_df, OOB_df=glb_newent_df)
```

```
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
```

```
## 
## Call:
## glm(formula = reformulate(indep_vars_vctr, response = lcl_predct_var), 
##     family = "binomial", data = fit_df)
## 
## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -1.642   0.000   0.000   0.000   1.239  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)
## (Intercept)  32.2548    25.4764   1.266    0.205
## Rasmussen    -1.1915     0.9141  -1.303    0.192
## SurveyUSA     3.6408     2.8720   1.268    0.205
## DiffCount     5.0957     4.4035   1.157    0.247
## PropR       -60.5197    50.9710  -1.187    0.235
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 138.2692  on 99  degrees of freedom
## Residual deviance:   4.9286  on 95  degrees of freedom
## AIC: 14.929
## 
## Number of Fisher Scoring iterations: 14
## 
##                                                      feats n.fit R.sq.fit
## 3 PropR, PropR:SurveyUSA, PropR:DiffCount, PropR:Rasmussen   100       NA
## 2                                                    PropR   100       NA
## 4                                                    PropR   100       NA
## 5                   Rasmussen, SurveyUSA, DiffCount, PropR   100       NA
## 1                                                   .rnorm   100       NA
##   R.sq.OOB Adj.R.sq.fit  SSE.fit SSE.OOB   AIC.fit   auc.fit   auc.OOB
## 3       NA           NA 164.9796      NA  23.66631 0.9967884 1.0000000
## 2       NA           NA 272.8743      NA  19.77232 0.9955841 0.9990079
## 4       NA           NA 272.8743      NA  19.77232 0.9955841 0.9990079
## 5       NA           NA 118.2244      NA  14.92861 0.9995986 0.9841270
## 1       NA           NA 401.4452      NA 140.26922 0.5000000 0.5000000
```

```r
# User specified
ret_lst <- myrun_mdl_fn(indep_vars_vctr=c("SurveyUSA", "DiffCount"),
                        glb_predct_var, glb_predct_var_name,
                        fit_df=glb_entity_df, OOB_df=glb_newent_df)
```

```
## 
## Call:
## glm(formula = reformulate(indep_vars_vctr, response = lcl_predct_var), 
##     family = "binomial", data = fit_df)
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -2.04741  -0.00977   0.00561   0.03751   1.32999  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)  
## (Intercept)  -0.6827     1.0468  -0.652   0.5143  
## SurveyUSA     0.3309     0.2226   1.487   0.1371  
## DiffCount     0.6619     0.3663   1.807   0.0708 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 138.269  on 99  degrees of freedom
## Residual deviance:  11.154  on 97  degrees of freedom
## AIC: 17.154
## 
## Number of Fisher Scoring iterations: 9
## 
##                                                      feats n.fit R.sq.fit
## 3 PropR, PropR:SurveyUSA, PropR:DiffCount, PropR:Rasmussen   100       NA
## 2                                                    PropR   100       NA
## 4                                                    PropR   100       NA
## 6                                     SurveyUSA, DiffCount   100       NA
## 5                   Rasmussen, SurveyUSA, DiffCount, PropR   100       NA
## 1                                                   .rnorm   100       NA
##   R.sq.OOB Adj.R.sq.fit  SSE.fit SSE.OOB   AIC.fit   auc.fit   auc.OOB
## 3       NA           NA 164.9796      NA  23.66631 0.9967884 1.0000000
## 2       NA           NA 272.8743      NA  19.77232 0.9955841 0.9990079
## 4       NA           NA 272.8743      NA  19.77232 0.9955841 0.9990079
## 6       NA           NA 181.4397      NA  17.15390 0.9977921 0.9940476
## 5       NA           NA 118.2244      NA  14.92861 0.9995986 0.9841270
## 1       NA           NA 401.4452      NA 140.26922 0.5000000 0.5000000
```

```r
# Simplify a model
# fit_df <- glb_entity_df; glb_mdl <- step(<complex>_mdl)

plot_models_df <- mutate(glb_models_df, feats.label=substr(feats, 1, 20))
if (glb_is_regression)
    print(myplot_scatter(plot_models_df, "Adj.R.sq.fit", "R.sq.OOB") + 
          geom_text(aes(label=feats.label), data=plot_models_df, color="NavyBlue", 
                    size=3.5, angle=45))

if (glb_is_classification) {
    # Lower AIC is better
    plot_models_df[, "inv.AIC.fit"] <- 1.0 / plot_models_df[, "AIC.fit"] 
    print(myplot_scatter(plot_models_df, "inv.AIC.fit", "auc.OOB") + 
          geom_text(aes(label=feats.label), data=plot_models_df, color="NavyBlue", 
                    size=3.5, angle=45))
}
```

![](RCP_files/figure-html/run_models-1.png) 

```r
script_df <- rbind(script_df, 
                   data.frame(chunk_label="fit_training.all", 
                              chunk_step_major=max(script_df$chunk_step_major)+1, 
                              chunk_step_minor=0))
print(script_df)
```

```
##                   chunk_label chunk_step_major chunk_step_minor
## 1                 import_data                1                0
## 2                cleanse_data                2                0
## 3        inspect_explore_data                2                1
## 4         manage_missing_data                2                2
## 5          encode_retype_data                2                3
## 6            extract_features                3                0
## 7             select_features                4                0
## 8  remove_correlated_features                4                1
## 9                  run_models                5                0
## 10           fit_training.all                6                0
```

## Step `6`: fit training.all

```r
print(mdl_feats_df <- myextract_mdl_feats(glb_sel_mdl, glb_entity_df))
```

```
##       Estimate Std. Error  z value         Pr.z    id fit.feat
## PropR 11.39044    3.15252 3.613122 0.0003025325 PropR     TRUE
```

```r
if (glb_is_regression) {
    ret_lst <- myrun_mdl_lm(indep_vars_vctr=mdl_feats_df$id,
                        glb_predct_var, glb_predct_var_name, fit_df=glb_entity_df)
    glb_sel_mdl <- glb_mdl    
#     print(glb_models_df[nrow(glb_models_df), ])
    glb_entity_df[, glb_predct_var_name] <- predict(glb_sel_mdl, newdata=glb_entity_df)
    print(myplot_scatter(glb_entity_df, glb_predct_var, glb_predct_var_name, 
                         smooth=TRUE))
    glb_entity_df[, paste0(glb_predct_var_name, ".err")] <- 
        abs(glb_entity_df[, glb_predct_var_name] - glb_entity_df[, glb_predct_var])
    print(head(orderBy(reformulate(c("-", paste0(glb_predct_var_name, ".err"))), 
                       glb_entity_df)))                             
}    

if (glb_is_classification) {
    ret_lst <- myrun_mdl_glm(indep_vars_vctr=mdl_feats_df$id,
                        glb_predct_var, glb_predct_var_name, fit_df=glb_entity_df)
    glb_sel_mdl <- glb_mdl        
#     print(glb_models_df[nrow(glb_models_df), ])
    glb_entity_df[, paste0(glb_predct_var_name, ".proba")] <- 
        predict(glb_sel_mdl, newdata=glb_entity_df, type="response")

    require(ROCR)
    ROCRpred <- prediction(glb_entity_df[, paste0(glb_predct_var_name, ".proba")],
                           glb_entity_df[, glb_predct_var])
    ROCRperf <- performance(ROCRpred, "tpr", "fpr")
    plot(ROCRperf, colorize=TRUE, print.cutoffs.at=seq(0, 1, 0.1), text.adj=c(-0.2,1.7))
    
    # 0 & 1 does not generate outcomes for certain categories
    thresholds_df <- data.frame(threshold=seq(0.0, 1.0, 0.1))
    thresholds_df$f.score <- sapply(1:nrow(thresholds_df), function(row_ix) 
        mycompute_classifier_f.score(glb_sel_mdl, glb_entity_df, 
                                     thresholds_df[row_ix, "threshold"], 
                                     glb_predct_var, glb_predct_var_name))
    print(thresholds_df)
    print(myplot_line(thresholds_df, "threshold", "f.score"))
    
    glb_clf_proba_threshold <- thresholds_df[which.max(thresholds_df$f.score), 
                                             "threshold"]
    # This should change to maximize f.score.OOB ???
    print(sprintf("Classifier Probability Threshold: %0.4f to maximize f.score.fit",
                  glb_clf_proba_threshold))

    glb_entity_df[, glb_predct_var_name] <- 
        (glb_entity_df[, paste0(glb_predct_var_name, ".proba")] >= 
             glb_clf_proba_threshold) * 1.0
    print(mycreate_xtab(glb_entity_df, c(glb_predct_var, glb_predct_var_name)))
    print(sprintf("f.score=%0.4f", 
        mycompute_classifier_f.score(glb_sel_mdl, glb_entity_df, 
                                     glb_clf_proba_threshold, 
                                     glb_predct_var, glb_predct_var_name)))    
}    
```

```
## 
## Call:
## glm(formula = reformulate(indep_vars_vctr, response = lcl_predct_var), 
##     family = "binomial", data = fit_df)
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -2.22880  -0.06541   0.10260   0.10260   1.37392  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -6.146      1.977  -3.108 0.001882 ** 
## PropR         11.390      3.153   3.613 0.000303 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 138.269  on 99  degrees of freedom
## Residual deviance:  15.772  on 98  degrees of freedom
## AIC: 19.772
## 
## Number of Fisher Scoring iterations: 8
## 
##                                                      feats n.fit R.sq.fit
## 3 PropR, PropR:SurveyUSA, PropR:DiffCount, PropR:Rasmussen   100       NA
## 2                                                    PropR   100       NA
## 4                                                    PropR   100       NA
## 6                                     SurveyUSA, DiffCount   100       NA
## 5                   Rasmussen, SurveyUSA, DiffCount, PropR   100       NA
## 1                                                   .rnorm   100       NA
## 7                                                    PropR   100       NA
##   R.sq.OOB Adj.R.sq.fit  SSE.fit SSE.OOB   AIC.fit   auc.fit   auc.OOB
## 3       NA           NA 164.9796      NA  23.66631 0.9967884 1.0000000
## 2       NA           NA 272.8743      NA  19.77232 0.9955841 0.9990079
## 4       NA           NA 272.8743      NA  19.77232 0.9955841 0.9990079
## 6       NA           NA 181.4397      NA  17.15390 0.9977921 0.9940476
## 5       NA           NA 118.2244      NA  14.92861 0.9995986 0.9841270
## 1       NA           NA 401.4452      NA 140.26922 0.5000000 0.5000000
## 7       NA           NA 272.8743      NA  19.77232 0.9955841        NA
```

![](RCP_files/figure-html/fit_training.all-1.png) 

```
##    threshold   f.score
## 1        0.0 0.6928105
## 2        0.1 0.9724771
## 3        0.2 0.9724771
## 4        0.3 0.9724771
## 5        0.4 0.9532710
## 6        0.5 0.9622642
## 7        0.6 0.9622642
## 8        0.7 0.9523810
## 9        0.8 0.9514563
## 10       0.9 0.9411765
## 11       1.0 0.0000000
```

![](RCP_files/figure-html/fit_training.all-2.png) 

```
## [1] "Classifier Probability Threshold: 0.1000 to maximize f.score.fit"
##   Republican Republican.predict.0 Republican.predict.1
## 1          0                   44                    3
## 2          1                   NA                   53
## [1] "f.score=0.9725"
```

```r
print(glb_feats_df <- mymerge_feats_Pr.z(glb_feats_df, glb_sel_mdl, glb_entity_df))
```

```
##          id     cor.y cor.y.abs cor.low         Pr.z
## 2     PropR 0.9484204 0.9484204       1 0.0003025325
## 1 DiffCount 0.8092777 0.8092777      NA           NA
## 3 Rasmussen 0.7929252 0.7929252      NA           NA
## 4 SurveyUSA 0.8101645 0.8101645      NA           NA
```

```r
# Most of this code is used again in predict_newdata chunk
glb_analytics_diag_plots <- function(obs_df) {
    for (var in subset(glb_feats_df, Pr.z < 0.1)$id) {
        plot_df <- melt(obs_df, id.vars=var, 
                        measure.vars=c(glb_predct_var, glb_predct_var_name))
#         if (var == "<feat_name>") print(myplot_scatter(plot_df, var, "value", 
#                                              facet_colcol_name="variable") + 
#                       geom_vline(xintercept=<divider_val>, linetype="dotted")) else     
            print(myplot_scatter(plot_df, var, "value", facet_colcol_name="variable"))
    }
    
    if (glb_is_regression) {
        plot_vars_df <- subset(glb_feats_df, Pr.z < 0.1)
        print(myplot_prediction_regression(obs_df, 
                    ifelse(nrow(plot_vars_df) > 1, plot_vars_df$id[2], ".rownames"), 
                                           plot_vars_df$id[1],
                    glb_predct_var, glb_predct_var_name)
#               + facet_wrap(reformulate(plot_vars_df$id[2])) # if [1,2] is a factor                                                         
#               + geom_point(aes_string(color="<col_name>.fctr")) #  to color the plot
              )
    }    
    
    if (glb_is_classification) {
        if (nrow(plot_vars_df <- subset(glb_feats_df, Pr.z < 0.1)) == 0)
            warning("No coefficients in selected model are statistically significant")
        else print(myplot_prediction_classification(obs_df, 
                    ifelse(nrow(plot_vars_df) > 1, plot_vars_df$id[2], ".rownames"),
                                               plot_vars_df$id[1],
                    glb_predct_var, glb_predct_var_name, glb_id_vars)
#               + geom_hline(yintercept=<divider_val>, linetype = "dotted")
                )
    }    
}
glb_analytics_diag_plots(obs_df=glb_entity_df)
```

![](RCP_files/figure-html/fit_training.all-3.png) 

```
##         State Year    .rnorm Rasmussen SurveyUSA DiffCount     PropR
## 28     Hawaii 2004 0.6300986        12         4         2 0.7500000
## 38    Indiana 2008 0.6300986         3         0         2 0.6250000
## 141 Wisconsin 2004 0.6300986        -1        -5         1 0.5333333
##     Republican Republican.predict.proba Republican.predict .rownames
## 28           0                0.9165725                  1        28
## 38           0                0.7256888                  1        38
## 141          0                0.4821911                  1       141
##     Republican.fctr Republican.predict.accurate         .label
## 28                0                       FALSE    Hawaii:2004
## 38                0                       FALSE   Indiana:2008
## 141               0                       FALSE Wisconsin:2004
```

![](RCP_files/figure-html/fit_training.all-4.png) 

```r
script_df <- rbind(script_df, 
                   data.frame(chunk_label="predict_newdata", 
                              chunk_step_major=max(script_df$chunk_step_major)+1, 
                              chunk_step_minor=0))
print(script_df)
```

```
##                   chunk_label chunk_step_major chunk_step_minor
## 1                 import_data                1                0
## 2                cleanse_data                2                0
## 3        inspect_explore_data                2                1
## 4         manage_missing_data                2                2
## 5          encode_retype_data                2                3
## 6            extract_features                3                0
## 7             select_features                4                0
## 8  remove_correlated_features                4                1
## 9                  run_models                5                0
## 10           fit_training.all                6                0
## 11            predict_newdata                7                0
```

## Step `7`: predict newdata

```r
if (glb_is_regression)
    glb_newent_df[, glb_predct_var_name] <- predict(glb_sel_mdl, 
                                        newdata=glb_newent_df, type="response")

if (glb_is_classification) {
    # Compute selected model predictions
    glb_newent_df[, paste0(glb_predct_var_name, ".proba")] <- 
        predict(glb_sel_mdl, newdata=glb_newent_df, type="response")
    glb_newent_df[, glb_predct_var_name] <- 
        (predict(glb_sel_mdl, newdata=glb_newent_df, type="response") >= 
            glb_clf_proba_threshold) * 1.0

    # Compute dummy model predictions
    glb_newent_df[, paste0(glb_predct_var, ".preddmy.proba")] <- 
        predict(glb_dmy_mdl, newdata=glb_newent_df, type="response")
    glb_newent_df[, paste0(glb_predct_var, ".preddmy")] <- 
        (predict(glb_dmy_mdl, newdata=glb_newent_df, type="response") >= 
            glb_clf_proba_threshold) * 1.0
}
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```r
myprint_df(glb_newent_df[, c(glb_id_vars, glb_predct_var, glb_predct_var_name)])
```

```
##          State Year Republican Republican.predict
## 7      Arizona 2012          1                  1
## 10    Arkansas 2012          1                  1
## 13  California 2012          0                  0
## 16    Colorado 2012          0                  0
## 19 Connecticut 2012          0                  0
## 24     Florida 2012          0                  1
##            State Year Republican Republican.predict
## 16      Colorado 2012          0                  0
## 19   Connecticut 2012          0                  0
## 33         Idaho 2012          1                  1
## 69   Mississippi 2012          1                  1
## 78      Nebraska 2012          1                  1
## 120 South Dakota 2012          1                  1
##             State Year Republican Republican.predict
## 126         Texas 2012          1                  1
## 129          Utah 2012          1                  1
## 134      Virginia 2012          0                  0
## 137    Washington 2012          0                  0
## 140 West Virginia 2012          1                  1
## 143     Wisconsin 2012          0                  0
```

```r
if (glb_is_regression) {
    print(sprintf("Total SSE: %0.4f", 
                  sum((glb_newent_df[, glb_predct_var_name] - 
                        glb_newent_df[, glb_predct_var]) ^ 2)))
    print(sprintf("RMSE: %0.4f", 
                  (sum((glb_newent_df[, glb_predct_var_name] - 
                        glb_newent_df[, glb_predct_var]) ^ 2) / nrow(glb_newent_df)) ^ 0.5))                        
    print(myplot_scatter(glb_newent_df, glb_predct_var, glb_predct_var_name, 
                         smooth=TRUE))
                         
    glb_newent_df[, paste0(glb_predct_var_name, ".err")] <- 
        abs(glb_newent_df[, glb_predct_var_name] - glb_newent_df[, glb_predct_var])
    print(head(orderBy(reformulate(c("-", paste0(glb_predct_var_name, ".err"))), 
                       glb_newent_df)))                                                      

#     glb_newent_df[, "<Output Pred variable>"] <- func(glb_newent_df[, glb_pred_var_name])                         
}                         

if (glb_is_classification) {
    ROCRpred <- prediction(glb_newent_df[, paste0(glb_predct_var_name, ".proba")],
                           glb_newent_df[, glb_predct_var])
    print(sprintf("auc=%0.4f", auc <- as.numeric(performance(ROCRpred, "auc")@y.values)))   
    print(sprintf("probability threshold=%0.4f", glb_clf_proba_threshold))
    
    print(mycreate_xtab(glb_newent_df, c(glb_predct_var, glb_predct_var_name)))
    print(sprintf("f.score.sel=%0.4f", 
        mycompute_classifier_f.score(mdl=glb_sel_mdl, obs_df=glb_newent_df, 
                                     proba_threshold=glb_clf_proba_threshold, 
                                     lcl_predct_var=glb_predct_var, 
                                     lcl_predct_var_name=glb_predct_var_name)))
    
    print(mycreate_xtab(glb_newent_df, c(glb_predct_var, paste0(glb_predct_var, ".preddmy"))))
    print(sprintf("f.score.dmy=%0.4f", 
        mycompute_classifier_f.score(mdl=glb_dmy_mdl, obs_df=glb_newent_df, 
                                     proba_threshold=glb_clf_proba_threshold, 
                                     lcl_predct_var=glb_predct_var, 
                                     lcl_predct_var_name=paste0(glb_predct_var, ".preddmy"))))
}    
```

```
## [1] "auc=0.9990"
## [1] "probability threshold=0.1000"
##   Republican Republican.predict.0 Republican.predict.1
## 1          0                   22                    2
## 2          1                   NA                   21
## [1] "f.score.sel=0.9545"
##   Republican Republican.preddmy.1
## 1          0                   24
## 2          1                   21
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```
## [1] "f.score.dmy=0.6364"
```

```r
glb_analytics_diag_plots(glb_newent_df)
```

![](RCP_files/figure-html/predict_newdata-1.png) 

```
##      State Year     .rnorm Rasmussen SurveyUSA DiffCount     PropR
## 24 Florida 2012 -0.2761841         2         0         6 0.6666667
## 42    Iowa 2012 -0.2761841         1        -2        -2 0.4000000
##    Republican Republican.predict.proba Republican.predict
## 24          0                0.8096071                  1
## 42          0                0.1693852                  1
##    Republican.preddmy.proba Republican.preddmy .rownames Republican.fctr
## 24                     0.53                  1        24               0
## 42                     0.53                  1        42               0
##    Republican.predict.accurate       .label
## 24                       FALSE Florida:2012
## 42                       FALSE    Iowa:2012
```

![](RCP_files/figure-html/predict_newdata-2.png) 

Null Hypothesis ($\sf{H_{0}}$): mpg is not impacted by am_fctr.  
The variance by am_fctr appears to be independent. 

```r
# print(t.test(subset(cars_df, am_fctr == "automatic")$mpg, 
#              subset(cars_df, am_fctr == "manual")$mpg, 
#              var.equal=FALSE)$conf)
```
We reject the null hypothesis i.e. we have evidence to conclude that am_fctr impacts mpg (95% confidence). Manual transmission is better for miles per gallon versus automatic transmission.


```
## R version 3.1.3 (2015-03-09)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.10.2 (Yosemite)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] ROCR_1.0-6      gplots_2.16.0   reshape2_1.4.1  mice_2.22      
##  [5] lattice_0.20-30 Rcpp_0.11.5     plyr_1.8.1      doBy_4.5-13    
##  [9] survival_2.38-1 ggplot2_1.0.1  
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6        caTools_1.17.1      colorspace_1.2-6   
##  [4] digest_0.6.8        evaluate_0.5.5      formatR_1.0        
##  [7] gdata_2.13.3        grid_3.1.3          gtable_0.1.2       
## [10] gtools_3.4.1        htmltools_0.2.6     KernSmooth_2.23-14 
## [13] knitr_1.9           labeling_0.3        MASS_7.3-39        
## [16] Matrix_1.1-5        munsell_0.4.2       nnet_7.3-9         
## [19] proto_0.3-10        randomForest_4.6-10 rmarkdown_0.5.1    
## [22] rpart_4.1-9         scales_0.2.4        splines_3.1.3      
## [25] stringr_0.6.2       tools_3.1.3         yaml_2.1.13
```
