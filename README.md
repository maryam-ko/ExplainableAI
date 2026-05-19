# Computational Reconstruction of Human Phospho-Regulatory Networks Using Machine Learning approaches 

## Overview
This MSc Bioinformatics research project investigates the reconstruction of human phospho-regulatory networks using machine learning approaches applied to large-scale phosphoproteomics datasets generated through mass spectrometry.

Phosphorylation plays a critical role in cellular signalling and regulation, yet many phospho-regulatory interactions remain poorly characterised. This project explores how machine learning methods can infer biologically meaningful phosphorylation relationships from high-dimensional phosphoproteomic data.

---

## Research Objectives
- Curate and process large-scale human phosphoproteomics datasets
- Reduce dimensionality and sparsity in phosphosite-level data
- Train machine learning models to predict phospho-regulatory interactions
- Evaluate predictive performance using cross-validation frameworks
- Compare predicted interactions against known biological interaction databases

---

## Dataset
- Integrated phosphoproteomics data from approximately **68 publicly available studies**
- Generated a phosphosite-level matrix containing:
  - **126,000+ phosphosites**
  - **11,000+ proteins**
- Data derived from mass spectrometry-based phosphoproteomics studies

---

## Methodology

### Data Processing
- Data cleaning and preprocessing
- Handling sparse phosphoproteomic datasets
- Feature selection using **Fisher Score**
- Dimensionality reduction and clustering approaches

### Machine Learning Models
#### Supervised Learning
- Linear Regression
- XGBoost

#### Unsupervised Techniques
- Clustering methods
- Dimensionality reduction

### Model Evaluation
- Cross-validation framework
- Performance assessment using:
  - R² coefficient values
  - Threshold-based model comparisons

---

## Biological Validation
To assess biological relevance, reconstructed interactions within:

- Mitogen-Activated Protein Kinase (MAPK) pathway
- Extracellular Signal-Regulated Kinase (ERK) signalling cascades

Predicted interactions were benchmarked against established protein interaction databases including:

- BioGRID
- OmniPath

This enabled identification of:
- Previously known regulatory interactions
- Potential novel phosphorylation relationships

---

## Technologies & Tools
- Python
- pandas
- NumPy
- scikit-learn
- XGBoost
- Matplotlib
- seaborn
- Jupyter Notebook
- Unix/Linux environments

---

## Key Outcomes
- Demonstrated that machine learning can effectively reconstruct phospho-regulatory networks from phosphoproteomics data
- Reduced high-dimensional phosphosite complexity using feature selection and clustering approaches
- Identified biologically relevant phosphorylation interactions supported by existing databases
- Established a computational framework for future phosphoproteomics and signalling network studies

---

## Future Improvements
- Integration of additional omics datasets
- Application of deep learning architectures
- Expansion of explainable AI methods 
- Experimental validation of predicted interactions

---

## Author
[Maryam Koddus]  
MSc Bioinformatics — Queen Mary University of London
