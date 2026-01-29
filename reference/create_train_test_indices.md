# Create train and test splits for cross-validation

Creates train and test splits for cross-validation by handling multiple
data types and supporting k-fold, leave-one-out (LOO), and
leave-percentage-out (LPO) methods. Handles missing values and maintains
data structure across multiple datasets.

## Usage

``` r
create_train_test_indices(
  data_list,
  cv_type = c("k-fold", "loo", "lpo"),
  k = 5,
  percentage = 20,
  number_folds = 10
)
```

## Arguments

- data_list:

  A list of datasets, one per likelihood. Each dataset can be a
  data.frame, SpatialPointsDataFrame, or metric_graph_data object

- cv_type:

  Type of cross-validation: "k-fold", "loo", or "lpo". Default is
  "k-fold"

- k:

  Number of folds for k-fold CV. Default is 5

- percentage:

  Training data percentage for LPO CV (1-99). Default is 20

- number_folds:

  Number of folds for LPO CV. Default is 10

## Value

A list where each element contains:

- train:

  Indices for training data mapped to original datasets

- test:

  Indices for test data mapped to original datasets

## Details

The function handles NA values by removing rows with any missing values
before creating splits. For multiple datasets, indices are mapped back
to their original positions in each dataset.
