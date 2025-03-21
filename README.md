# GNN-PINN

## 3-Stage Flow: Thermodynamic Property Prediction via GNN & PINN

This repository implements a multi-stage training pipeline to predict thermodynamic properties using graph neural networks (GNNs) combined with physics-informed neural networks (PINNs). The model leverages both non-equilibrium and equilibrium datasets through a three-stage process, enabling robust prediction performance with custom loss functions that enforce physical consistency.

## Overview

The training pipeline consists of three distinct stages:

## Training Pipeline

1. **Pretraining (Non-equilibrium Data):**  
   Train the base model using non-equilibrium data with an MSE loss function to obtain graph embeddings.

2. **Saturation Pressure Model Training (Equilibrium Data):**  
   Using the previously obtained graph embeddings and equilibrium data, train a separate saturation pressure model with an MSE loss function to perform predictions.

3. **PINN Finetuning (Equilibrium Data):**  
   Finetune the base model using a custom PINN loss that incorporates both MSE and a PDE-based term. In this stage, the predicted saturation pressure from Stage 2 is used as an input to enforce thermodynamic consistency.

## Key Features

- **Data Preparation & Normalization:**  
  - Loads a CSV dataset containing molecular and thermodynamic data.
  - Splits the dataset into equilibrium and non-equilibrium subsets.
  - Applies z-score normalization to the key variables. 

- **Custom GNN Architecture:**  
  - **GCNLayer:** Implements graph convolution operations with optional dimensionality reduction.
  - **GRLayer:** Performs global pooling (mean, max, min, variance, sum, log-sum-exp) to aggregate node features.
  - Combines graph-based features with thermodynamic inputs via dense layers.

- **Custom PINN Loss Function:**  
  - Computes MSE for property predictions.
  - Adds a PDE penalty term for equilibrium data to enforce physical constraints.

- **Multi-Stage Training Process:**  
  - **Stage 1:** Pretrain the base model on non-equilibrium data.
  - **Stage 2:** Train the saturation pressure model on equilibrium data.
  - **Stage 3:** Finetune the base model using the predicted saturation pressure from Stage 2, with PINN loss.


## ⚙️Requirements

numpy==1.24.3

pandas==2.0.3

rdkit==2024.3.5

tensorflow==2.10.0

matplotlib==3.7.5

Please install CUDA and TensorFlow according to your system configuration. Ensure you follow the official guidelines to match your hardware and operating system requirements.

```bash
pip install -r requirements.txt
```

├── dataset.csv               # Input CSV file containing the dataset
├── README.md                 # This file
└──  GNN_PINN.py              # Main script 

