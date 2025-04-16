📘 Project Overview
This project applies machine learning to identify bioactive compounds that inhibit Acetylcholinesterase (AChE), an important target in treating neurological disorders like Alzheimer’s disease and myasthenia gravis. It includes both a machine learning pipeline and a web application for predicting the bioactivity (pIC50) of chemical compounds using their SMILES notations.

🔄 Step-by-Step Workflow
1. Data Collection
Retrieved bioactivity data from the ChEMBL database using chembl_webresource_client.

Filtered compounds targeting human Acetylcholinesterase.

2. Bioactivity Evaluation
Converted IC50 values to pIC50 (logarithmic scale) for easier interpretation.

Classified compounds as active, inactive, or intermediate based on pIC50 thresholds.

3. Data Preprocessing
Removed missing and duplicate entries.

Ensured high-quality and clean dataset.

4. Exploratory Data Analysis
Performed Mann-Whitney U tests on Lipinski descriptors (molecular weight, logP, H-bond donors/acceptors).

Visualized key trends using boxplots and scatter plots with Seaborn.

5. Feature Engineering
Generated molecular fingerprints using RDKit and PaDEL-Descriptor.

Used fingerprints like PubChem, Klekota-Roth, and CDK Extended.

6. Feature Selection
Reduced dimensionality using variance thresholding to remove low-variance features.

7. Model Training and Selection
Trained models: Random Forest, Support Vector Regressor (SVR), XGBoost.

SVR gave the best performance with R² = 0.7512 on test data.

8. Model Evaluation
Evaluated models using RMSE and MAE.

Created scatter plots to compare predicted vs actual pIC50 values.

🚀 Bioactivity Prediction App
How it Works
Users input SMILES strings of chemical compounds.

The app generates molecular descriptors and predicts pIC50 using the trained model.

Predictions are displayed and downloadable as a CSV file.

Key Features
Predicts pIC50 values for up to 100 compounds at a time.

Uses combined molecular fingerprints for better accuracy.

Built using Streamlit for a simple and interactive UI.

Outputs are instantly downloadable for further analysis.

🛠 Technologies Used
Python

Machine Learning (SVR, Random Forest, XGBoost)

Streamlit (Web App)

RDKit and PaDEL (Descriptor Generation)

Pickle (Model Serialization)


