# 🧬 Molecular AI Dashboard

**Molecular AI Dashboard** is an intelligent platform designed to accelerate early-stage drug discovery by combining molecular property prediction with AI-assisted molecule generation.

Starting from a **SMILES** representation of a molecule, the platform predicts key physicochemical properties, evaluates drug-likeness according to Lipinski's Rule of Five, and assists researchers in identifying promising drug candidates through an interactive and user-friendly dashboard.

The project also explores the chemical space of the **ZINC250k** dataset and demonstrates how machine learning models can efficiently support molecular analysis and candidate selection.

---

# 🎯 Project Objectives

The main objectives of this project are to:

- Predict three important molecular properties: **logP**, **QED**, and **SAS**.
- Explore the chemical space of the ZINC250k molecular dataset.
- Generate novel drug-like molecules for virtual screening.
- Evaluate molecular drug-likeness using Lipinski's Rule of Five.
- Provide an interactive dashboard for real-time molecular analysis and visualization.

---

# 📌 Project Poster

The following scientific poster summarizes the entire project, including the motivation, objectives, data preparation, machine learning pipeline, experimental evaluation, chemical space exploration, candidate generation, and the final web application.

<p align="center">
  <img src="images/Poster.png" width="900">
</p>

---

# 🖥️ Interactive Dashboard

A complete web dashboard was developed to simplify molecular analysis and provide an intuitive interface for researchers.

The application allows users to analyze molecules in real time by simply entering a **SMILES** string, making the prediction process both fast and accessible.

## 🔹 Property Prediction Interface

Users can:

- Enter a SMILES molecule.
- Generate molecular fingerprints automatically.
- Predict **logP**, **QED**, and **SAS** simultaneously.
- Evaluate the molecule using Lipinski's Rule of Five.
- Display the molecular structure instantly.
- Monitor the prediction process through an execution log.

<p align="center">
  <img src="images/imagesdashboard_prediction.png" width="900">
</p>

---

## 🔹 Prediction Results

After processing the input molecule, the dashboard displays:

- Predicted values for **logP**, **QED**, and **SAS**.
- A visualization of the molecular structure.
- Drug-likeness assessment based on Lipinski's criteria.
- Real-time prediction logs.
- A clean and interactive interface for rapid decision-making.

<p align="center">
  <img src="images/dashboard_results.png" width="900">
</p>

---

# 📈 Key Results

The experimental evaluation demonstrates that the proposed workflow can accurately predict molecular properties while efficiently exploring the chemical space.

The developed platform successfully:

- Predicts molecular properties with high accuracy.
- Identifies promising drug-like candidates.
- Facilitates rapid molecular screening through an interactive interface.
- Provides an end-to-end workflow, from molecular input to property prediction and candidate analysis.

---

# 👩‍💻 Author

**Amina Bouazza**

Master's in Artificial Intelligence

Faculty of Sciences Ben M'Sick

Hassan II University of Casablanca
