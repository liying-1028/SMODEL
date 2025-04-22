# SMODEL


## Usage
#### Clone this repo.
```matlab
git clone https://github.com/liying-1028/SMODEL.git
cd SMODEL
```

#### Code description

- ```Data_processing.py```: Data preprocessing  
- ```Run_SMODEL_Simulated_data.m```: SMODEL for simulated spatial triple-omics  datasets
- ```Run_SMODEL_real_data.m```: SMODEL  for real datasets
- ```visualization.ipynb```: Clustering visualization

#### Example 

Take the dataset "simulated  spatial triple-omics" as an example

- Step 1: You can preprocess the data and prepare the results (poolsget) of the base clustering method  as input for the model training:

  ```python
  Data_processing.py
  ```

- Step 2: Main running pipeline:

  ```matlab
  Run_SMODEL_Simulated_data.m
  ```

- Step 3: Visualization:

  ```python
  visualization.ipynb
  ```

  #### Note: 

  - To reduce your time, we have uploaded our preprocessed data into the folder 'data'. You can perform the corresponding steps selectively.
  - To reproduce the result, you should use the default parameters. For new datasets, we suggest using the default parameter configuration (λ = 10, β = 0.2, anchor = 100, KNN = 12; refer to:```Run_SMODEL_real_data.m```).  Additionally, you may adjust these parameters according to the unique characteristics of your dataset to attain the best possible outcomes.
