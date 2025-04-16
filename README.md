🌟 Drug Discovery Using Machine Learning + CholinEase App 🚀
Welcome to my project where Machine Learning meets Drug Discovery! This repository includes a comprehensive workflow for predicting the bioactivity of small molecules targeting Acetylcholinesterase (AChE) — a key enzyme in neurological disorders like Alzheimer’s disease.

📌 Project Highlights
🧪 1. Objective
Use machine learning models to predict pIC50 values (a measure of bioactivity) for potential drug-like compounds that inhibit Acetylcholinesterase.

🔍 Step-by-Step Overview
📥 2. Data Collection
Data sourced from ChEMBL using chembl_webresource_client.

Focused only on bioactivity data for human Acetylcholinesterase.

⚙️ 3. Data Preprocessing
Cleaned and curated dataset by:

Removing duplicates and missing values.

Converting IC50 values to pIC50 for better interpretability.

Categorized compounds as active, inactive, or intermediate.

📊 4. Exploratory Data Analysis
Visualized distribution of bioactivity classes and molecular properties.

Used Mann-Whitney U Test to find significant differences between active and inactive compounds across:

Molecular weight

logP

Hydrogen bond donors/acceptors

🤖 5. Machine Learning Workflow
🧬 Feature Engineering
Generated molecular descriptors and fingerprints using RDKit and PaDEL-Descriptor.

📉 Feature Selection
Applied variance thresholding to remove low-informative features.

🏗️ Model Training
Models trained:

Random Forest

Support Vector Regressor (SVR)

XGBoost

SVR gave the best results with R² = 0.7512

📈 Model Evaluation
Evaluated with:

RMSE

MAE

Plotted predicted vs actual pIC50 values for visual inspection.

🚀 Bioactivity Prediction App - CholinEase
💻 Built with:
Python, Streamlit, PaDELpy, pickle

🧠 Key Features:
Input: SMILES notations for up to 100 compounds

Output: Predicted pIC50 values in seconds

Uses three fingerprinting methods (PubChem, Klekota-Roth, CDK Extended)

Results downloadable in CSV format

📸 App Preview:
<p align="center"> <img src="CholinEase - Bioactivity Prediction App.pdf" alt="App Screenshot" width="600"/> </p>
