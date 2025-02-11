import os
import pickle
import base64
import subprocess
import pandas as pd
import streamlit as st
from padelpy import padeldescriptor

os.environ["JAVA_HOME"] = r"C:\Program Files\Java\jdk1.8.0_101"
os.environ["PATH"] += os.pathsep + os.path.join(os.environ["JAVA_HOME"], "bin")

# ✅ Set Streamlit Page Config
st.set_page_config(page_title="CholinEase - Bioactivity Prediction App", layout="centered")

# ✅ Stylish Header
st.markdown(
    """
    <h1 style='text-align: center; font-size: 50px; color: #4CAF50;'>CholinEase</h1>
    <h3 style='text-align: center; color: #666;'>Bioactivity Prediction App</h3>
    """, 
    unsafe_allow_html=True
)
st.image("drug_discovery1.jpg", use_container_width=True)
st.markdown(""" 
- App built in `Python` + `Streamlit` by Miron Shah  
---
""", unsafe_allow_html=True)

st.header("Predict pIC50")

# ✅ Function to Check Java Installation
def check_java():
    java_path = subprocess.getoutput("where java").strip()
    if "not found" in java_path.lower() or not java_path:
        st.error("❌ Java not found! Install Java and restart the app.")
        return False
    
    java_home = os.path.dirname(os.path.dirname(java_path))
    os.environ["JAVA_HOME"] = java_home
    os.environ["PATH"] += os.pathsep + os.path.join(java_home, "bin")
    
    st.success(f"✅ Java detected: {java_home}")
    return True

java_installed = check_java()

# ✅ Descriptor Calculation Function
def desc_calc():
    fingerprints = ['PubChem', 'KlekotaRoth', 'CDKextended']
    fingerprint_descriptortypes = [
        'PubchemFingerprinter.xml',
        'KlekotaRothFingerprinter.xml',
        'ExtendedFingerprinter.xml'
    ]
    output_files = {fp: f"{fp}_app_data.csv" for fp in fingerprints}

    for fp, xml_file in zip(fingerprints, fingerprint_descriptortypes):
        try:
            print(f"Generating descriptors for: {fp}")
            padeldescriptor(
                mol_dir='molecule.smi',  
                d_file=output_files[fp],  
                descriptortypes=xml_file,  
                detectaromaticity=True,
                standardizenitro=True,
                standardizetautomers=True,
                threads=2,
                removesalt=True,
                log=True,
                fingerprints=True
            )
            print(f"{fp} descriptors saved: {output_files[fp]}")
        except Exception as e:
            st.error(f"❌ Error generating {fp} descriptors: {e}")

# ✅ Function to Download CSV File
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# ✅ Model Prediction Function
def build_model(input_data, manual_data):
    try:
        load_model = pickle.load(open('acetylcholinesterase_model.pkl', 'rb'))
    except FileNotFoundError:
        st.error("❌ Model file 'acetylcholinesterase_model.pkl' not found! Upload the correct file.")
        return
    
    prediction = load_model.predict(input_data)
    st.header('**Predicted Output**')

    prediction_output = pd.Series(prediction, name='pIC50')
    notation = pd.Series(manual_data['smiles'], name='SMILES Notation')
    df = pd.concat([notation, prediction_output], axis=1)

    df.index = df.index + 1
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)

# ✅ User Input Section
smiles_input = st.text_area(
    "Please enter SMILES notations below (One per Line, Up to 100 Compounds)",
    height=250,
    placeholder="e.g. C1CCCCC1\nCC(=O)O\nO=C(C)O"
)

if st.button("Predict"):
    if not java_installed:
        st.error("❌ Java is not installed. Install Java and restart the app.")
    elif smiles_input:
        smiles_list = smiles_input.strip().split('\n')[:100]
        if len(smiles_list) > 0:
            manual_data = pd.DataFrame({'smiles': smiles_list})

            with st.spinner("🔄 Predicting..."):
                with open('molecule.smi', 'w') as f:
                    f.write('\n'.join(smiles_list))
                
                desc_calc()

                # ✅ Check if descriptor files exist before proceeding
                descriptor_files = [
                    'PubChem_app_data.csv',
                    'KlekotaRoth_app_data.csv',
                    'CDKextended_app_data.csv'
                ]
                
                if not all(os.path.exists(f) for f in descriptor_files):
                    st.error("❌ Descriptor files are missing! Ensure PaDEL-Descriptor ran correctly.")
                else:
                    descriptors = pd.concat([pd.read_csv(file) for file in descriptor_files], axis=1)
                    descriptors.to_csv("Combined_PubChem_CDK_Klekota_app_data.csv", index=False)

                    desc = pd.read_csv('Combined_PubChem_CDK_Klekota_app_data.csv')
                    try:
                        Xlist = list(pd.read_csv('Combined_PubChem_CDK_Klekota_modified.csv').columns)
                        desc_subset = desc[Xlist]
                        build_model(desc_subset, manual_data)
                    except FileNotFoundError:
                        st.error("❌ Descriptor reference file is missing! Ensure it is in the correct directory.")
        else:
            st.warning("⚠️ Please enter at least one SMILES notation.")
    else:
        st.warning("⚠️ Please enter SMILES strings to start the prediction!")
