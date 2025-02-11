import os
import base64
import pickle
import pandas as pd
import streamlit as st
from padelpy import padeldescriptor

# Ensure Java is set before calling PaDEL-Descriptor
os.environ["JAVA_HOME"] = "/usr/lib/jvm/java-11-openjdk-amd64"  # Linux
# os.environ["JAVA_HOME"] = "C:\\Program Files\\Java\\jre1.8.0_xx"  # Windows

# Ensure Java is in the PATH
os.environ["PATH"] += os.pathsep + os.path.join(os.environ["JAVA_HOME"], "bin")

# Check Java availability in Streamlit
java_check = os.popen("java -version").read()
st.text(f"Java Check:\n{java_check}")

# App Header
st.set_page_config(page_title="CholinEase - Bioactivity Prediction App", layout="centered")

st.markdown(
    """
    <h1 style='text-align: center; font-size: 50px; color: #4CAF50;'>CholinEase</h1>
    <h3 style='text-align: center; color: #666;'>Bioactivity Prediction App</h3>
    """,
    unsafe_allow_html=True
)

try:
    st.image("drug_discovery1.jpg", use_container_width=True)
except FileNotFoundError:
    st.error("Error: 'drug_discovery1.jpg' not found. Please ensure the file exists.")

st.markdown(""" 
- App built in `Python` + `Streamlit` by Miron Shah
--- 
""", unsafe_allow_html=True)

st.header("Predict pIC50")

def desc_calc():
    """Calculate molecular descriptors using PaDEL-Descriptor."""
    fingerprints = ['PubChem', 'KlekotaRoth', 'CDKextended']
    fingerprint_descriptortypes = [
        'PubchemFingerprinter.xml',
        'KlekotaRothFingerprinter.xml',
        'ExtendedFingerprinter.xml'
    ]
    output_files = {fp: f"{fp}_app_data.csv" for fp in fingerprints}
    
    for fp, xml_file in zip(fingerprints, fingerprint_descriptortypes):
        try:
            st.write(f"Generating descriptors for: {fp}")
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
            st.success(f"{fp} descriptors saved successfully!")
        except Exception as e:
            st.error(f"Error generating descriptors for {fp}: {e}")

def filedownload(df):
    """Generate a downloadable link for CSV output."""
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

def build_model(input_data, manual_data):
    """Load ML model and make predictions."""
    try:
        load_model = pickle.load(open('acetylcholinesterase_model.pkl', 'rb'))
        prediction = load_model.predict(input_data)

        st.header('**Predicted Output**')
        prediction_output = pd.Series(prediction, name='pIC50')
        notation = pd.Series(manual_data['smiles'], name='SMILES Notation')
        df = pd.concat([notation, prediction_output], axis=1)

        df.index = df.index + 1
        st.write(df)
        st.markdown(filedownload(df), unsafe_allow_html=True)
    except FileNotFoundError:
        st.error("Error: Model file 'acetylcholinesterase_model.pkl' not found.")
    except Exception as e:
        st.error(f"Error in model prediction: {e}")

smiles_input = st.text_area(
    "Please enter SMILES notations below (One per Line, Up to 100 Compounds)",
    height=250,
    placeholder="e.g. C1CCCCC1\nCC(=O)O\nO=C(C)O"
)

if st.button("Predict"):
    if smiles_input:
        smiles_list = smiles_input.strip().split('\n')[:100]
        if len(smiles_list) > 0:
            manual_data = pd.DataFrame({'smiles': smiles_list})

            with st.spinner("Predicting..."):
                with open('molecule.smi', 'w') as f:
                    f.write('\n'.join(smiles_list))

                desc_calc()

                try:
                    descriptors_to_combine = [
                        'PubChem_app_data.csv',
                        'KlekotaRoth_app_data.csv',
                        'CDKextended_app_data.csv'
                    ]
                    descriptors = pd.concat([pd.read_csv(file) for file in descriptors_to_combine], axis=1)
                    descriptors.to_csv("Combined_PubChem_CDK_Klekota_app_data.csv", index=False)

                    desc = pd.read_csv('Combined_PubChem_CDK_Klekota_app_data.csv')
                    Xlist = list(pd.read_csv('Combined_PubChem_CDK_Klekota_modified.csv').columns)
                    desc_subset = desc[Xlist]

                    build_model(desc_subset, manual_data)
                except FileNotFoundError as e:
                    st.error(f"File error: {e}")
                except Exception as e:
                    st.error(f"Unexpected error: {e}")
        else:
            st.warning("Please enter at least one SMILES notation.")
    else:
        st.warning("Please enter SMILES strings to start the prediction!")




