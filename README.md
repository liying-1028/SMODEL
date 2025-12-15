# SMODEL

## Requirements

Before running SMODEL, please ensure the following MATLAB toolboxes are installed:

- **Statistics and Machine Learning Toolbox**
- **Deep Learning Toolbox**

## Usage
### Clone this repo
```matlab
git clone https://github.com/liying-1028/SMODEL.git
cd SMODEL
```
  
### Code description

- ```Data_processing.py```: Data preprocessing  
- ```Run_SMODEL_Simulated_data.m```: SMODEL for simulated spatial triple-omics  datasets
- ```Run_SMODEL_real_data.m```: SMODEL  for real datasets
- ```visualization.ipynb```: Clustering visualization

### Example 

Take the dataset "simulated  spatial triple-omics" as an example

- Step 1: Prepare the data. The SMODEL takes the preprocessed data and the base clustering results (`poolsget`) as input. For data preprocessing, please refer to the following script:

  ```python
  Data_processing.py
  ```
After preprocessing, the data should be organized and saved into a MATLAB `.mat` file, which will be directly used by SMODEL. We provide an example MATLAB script, `DataMat.m`, for generating this file.


- Step 2: Main running pipeline:

  ```matlab
  Run_SMODEL_Simulated_data.m
  ```

- Step 3: Visualization:

  ```python
  visualization.ipynb
  ```


### Note: 

- To reduce the time required for reproduction, we have uploaded our preprocessed data into the folder `data`. Users may directly use these files and skip Step 1 if desired.
- To reproduce the result, you should use the default parameters. For new datasets, we recommend starting with the default configuration (see `Run_SMODEL_real_data.m`).  Additionally, you may adjust the parameters according to the unique characteristics of your dataset to attain the best possible outcomes.
