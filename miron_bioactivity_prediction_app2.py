import os
import base64
import pickle
import pandas as pd
import streamlit as st
import subprocess
from padelpy import padeldescriptor
import requests

# ✅ Set Streamlit page config FIRST
st.set_page_config(page_title="CholinEase - Bioactivity Prediction App", layout="centered")

# ✅ Check if Java is installed and set the environment variable
def check_java():
    java_path = subprocess.getoutput("which java")
    if not java_path:
        st.error("❌ Java not found! Install Java and restart the app.")
        return False
    os.environ["JAVA_HOME"] = os.path.dirname(os.path.dirname(java_path))
    os.environ["PATH"] += os.pathsep + os.path.join(os.environ["JAVA_HOME"], "bin")
    return True

java_installed = check_java()

# ✅ App Header with Stylish Title
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
            print(f"{fp} descriptors saved in: {output_files[fp]}")
        except Exception as e:
            st.error(f"Error generating {fp} descriptors: {e}")

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
    "Please enter SMILES n

