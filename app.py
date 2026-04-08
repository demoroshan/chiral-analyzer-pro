import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io
import requests
from stmol import showmol
import py3Dmol

st.set_page_config(page_title="Chiral Analyzer Pro", page_icon="🧪", layout="wide")

# Helper function to get SMILES and IUPAC from Name via PubChem
@st.cache_data
def fetch_compound_data(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES,IUPACName/JSON"
        response = requests.get(url)
        data = response.json()['PropertyTable']['Properties'][0]
        return data.get('CanonicalSMILES'), data.get('IUPACName', 'N/A')
    except:
        return None, None

st.title("🧪 Molecule Chiral Center Analyzer Pro")
st.markdown("Analyze compounds from PubChem, explore 2D/3D structures, learn stereochemistry with guided theory, and test your knowledge!")

# Compound Selection Setup
st.subheader("Select a Molecule")
compound_list = [
    "Ibuprofen", 
    "Thalidomide", 
    "Cholesterol", 
    "D-Glucose", 
    "Penicillin G",
    "Carvone", 
    "Limonene", 
    "Menthol", 
    "Camphor", 
    "Epinephrine",
    "Morphine", 
    "Paclitaxel", 
    "Amphetamine", 
    "Custom (Type below)"
]

selected_option = st.selectbox("Choose a famous chiral compound to analyze:", compound_list)

# Handle Custom Input
if selected_option == "Custom (Type below)":
    compound_name = st.text_input("Enter Compound Name:", "Caffeine")
else:
    compound_name = selected_option

# Fetch Data
smiles, iupac_name = fetch_compound_data(compound_name)

if smiles:
    mol = Chem.MolFromSmiles(smiles)
    mol_3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
    
    st.header(f"Compound Analysis: {compound_name.capitalize()}")
    st.caption(f"**IUPAC Name:** {iupac_name}")

    # Calculations
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.MolWt(mol)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    # Layout: Tabs for Organization
    tab1, tab2, tab3, tab4 = st.tabs(["📊 2D Structure & Filters", "🧊 3D Interactive Model", "📖 Stereochemistry Theory", "🧠 Knowledge Quiz"])

    with tab1:
        col1, col2 = st.columns([1.5, 1])
        
        with col2:
            st.subheader("Filter Chiral Centers")
            view_option = st.radio("Select which stereocenters to highlight:", ["All (R/S)", "Only R", "Only S"])
            
            st.subheader("Molecular Metadata")
            st.success(f"**Formula:** {formula}")
            st.metric("Molecular Weight", f"{mw:.2f} g/mol")
            st.write(f"**Chiral Centers Found:** {len(chiral_centers)}")
            
            # Setup Coloring Logic based on user filter
            atom_colors = {}
            for idx, config in chiral_centers:
                if config == 'R' and view_option in ["All (R/S)", "Only R"]:
                    atom_colors[idx] = (1.0, 0.7, 0.7) # Light Red
                    st.error(f"Carbon {idx}: **R** (Rectus)")
                elif config == 'S' and view_option in ["All (R/S)", "Only S"]:
                    atom_colors[idx] = (0.7, 0.7, 1.0) # Light Blue
                    st.info(f"Carbon {idx}: **S** (Sinister)")

        with col1:
            st.subheader("2D Visual Analysis")
            # Create High-Quality Drawing
            d2d = rdMolDraw2D.MolDraw2DCairo(600, 600)
            d2d.drawOptions().addAtomIndices = True 
            d2d.DrawMolecule(mol, highlightAtoms=list(atom_colors.keys()), highlightAtomColors=atom_colors)
            d2d.FinishDrawing()
            img_2d = Image.open(io.BytesIO(d2d.GetDrawingText()))
            st.image(img_2d, use_container_width=True)
            st.caption("🔴 Red = R Configuration | 🔵 Blue = S Configuration")

    with tab2:
        st.subheader("3D Stereochemical Model")
        spin = st.radio("3D Molecule Spin:", ["Spin", "Do Not Spin"], horizontal=True)
        
        # Create py3Dmol view
        mblock = Chem.MolToMolBlock(mol_3d)
        view = py3Dmol.view(width=800, height=500)
        view.addModel(mblock, 'mol')
        view.setStyle({'stick': {}, 'sphere': {'radius': 0.3}})
        
        if spin == "Spin":
            view.spin(True)
        else:
            view.spin(False)
            
        view.zoomTo()
        showmol(view, height=500, width=800)
        st.info("Drag to rotate manually, scroll to zoom in/out.")

    with tab3:
        st.markdown("### 🧬 What is a Chiral Center?")
        st.write("""
        A **chiral center** (or stereocenter) is an atom—most commonly carbon—that is bonded to **four unique groups**. 
        Because of its tetrahedral 3D geometry, the molecule becomes asymmetric. This means the molecule cannot be perfectly overlaid onto its mirror image, much like your left and right hands.
        """)
        
        st.markdown("### 🔄 The R/S Configuration System")
        st.write("""
        Chemists use the **Cahn-Ingold-Prelog (CIP) Priority Rules** to name these centers so everyone knows exactly which 3D arrangement is being discussed.
        
        **How it works:**
        1. **Assign Priority:** Look at the four atoms directly attached to the chiral center. Rank them 1 to 4 based on their atomic number (higher atomic number = higher priority).
        2. **Orient the Molecule:** Visualize the molecule so that the lowest priority group (usually Hydrogen) is pointing away from you (dashed line).
        3. **Trace the Path:** Draw a curve from Priority 1 → Priority 2 → Priority 3.
        """)
        
        col_r, col_s = st.columns(2)
        with col_r:
            st.error("**R (Rectus)**\n\nIf the path travels **clockwise** (to the right), the center is assigned the **R** configuration.")
        with col_s:
            st.info("**S (Sinister)**\n\nIf the path travels **counter-clockwise** (to the left), the center is assigned the **S** configuration.")

        st.divider()
        st.markdown(f"### 🧪 About {compound_name.capitalize()}")
        st.write(f"**IUPAC Name:** {iupac_name}")
        st.write(f"This compound has **{len(chiral_centers)}** identified chiral centers. In pharmacology and biology, the specific R/S configuration at each of these centers dictates how the molecule interacts with enzymes and receptors in the body. An incorrect stereocenter can render a drug completely inactive or highly toxic!")

    with tab4:
        st.subheader("🧠 Stereochemistry Knowledge Check")
        
        # Simple state-based quiz
        q1 = st.radio("1. What determines the priority of groups in the Cahn-Ingold-Prelog (CIP) system?", 
                      ["Alphabetical order of the element names", "Atomic mass", "Atomic number", "Electronegativity"], index=None)
        q2 = st.radio("2. If the lowest priority group is pointing TOWARDS you, and the path from 1 -> 2 -> 3 is clockwise, what is the actual configuration?", 
                      ["R", "S", "Achiral", "Meso"], index=None)
        
        if st.button("Check Answers"):
            if q1 == "Atomic number":
                st.success("Question 1 Correct! Priority is based on atomic number.")
            elif q1:
                st.error("Question 1 Incorrect. Hint: Think about protons.")
                
            if q2 == "S":
                st.success("Question 2 Correct! If the lowest priority group is pointing towards you, you must REVERSE the result. Clockwise is usually R, so reversed it becomes S.")
            elif q2:
                st.error("Question 2 Incorrect. Hint: Remember the rule of reversing the result when the lowest priority group is on a wedge!")

else:
    st.error("Could not find that compound. Please check the spelling or try entering a specific chemical name.")
