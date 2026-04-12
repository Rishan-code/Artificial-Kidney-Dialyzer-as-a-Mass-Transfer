import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import os
from sklearn.metrics import r2_score

def create_model_plots():
    # Paths
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(base_dir, 'data', 'dialyzer_dataset.csv')
    model_path = os.path.join(base_dir, 'data', 'rf_surrogate.pkl')
    output_dir = os.path.join(base_dir, 'plots')
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load Data and Model
    if not os.path.exists(data_path) or not os.path.exists(model_path):
        print("Data or model file missing. Please ensure src/ml_model.py has been run.")
        return

    df = pd.read_csv(data_path)
    with open(model_path, 'rb') as f:
        model = pickle.load(f)
        
    X = df[['Qb_ml_min', 'Qd_ml_min', 'Quf_ml_min']]
    y = df[['Clearance_Urea', 'Clearance_B12']]
    
    # Preds for whole dataset
    y_pred = model.predict(X)
    
    sns.set_theme(style="whitegrid")
    
    # 1. Feature Importance
    plt.figure(figsize=(8, 5))
    importances = model.feature_importances_
    features = ['Blood Flow (Qb)', 'Dialysate Flow (Qd)', 'Ultrafiltration (Quf)']
    sns.barplot(x=importances, y=features, palette="viridis")
    plt.title('Feature Importances in Predicting Clearance')
    plt.xlabel('Importance')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'feature_importance.png'), dpi=300)
    plt.close()
    
    # 2. Actual vs Predicted - Urea
    plt.figure(figsize=(7, 6))
    sns.scatterplot(x=y['Clearance_Urea'], y=y_pred[:, 0], alpha=0.5, color='blue')
    plt.plot([y['Clearance_Urea'].min(), y['Clearance_Urea'].max()],
             [y['Clearance_Urea'].min(), y['Clearance_Urea'].max()],
             'r--', lw=2)
    plt.title(f'Urea Clearance: Actual vs Predicted\nR² = {r2_score(y["Clearance_Urea"], y_pred[:, 0]):.4f}')
    plt.xlabel('Actual Clearance (mL/min)')
    plt.ylabel('Predicted Clearance (mL/min)')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'actual_vs_predicted_urea.png'), dpi=300)
    plt.close()

    # 3. Actual vs Predicted - VitB12
    plt.figure(figsize=(7, 6))
    sns.scatterplot(x=y['Clearance_B12'], y=y_pred[:, 1], alpha=0.5, color='green')
    plt.plot([y['Clearance_B12'].min(), y['Clearance_B12'].max()],
             [y['Clearance_B12'].min(), y['Clearance_B12'].max()],
             'r--', lw=2)
    plt.title(f'Vitamin B12 Clearance: Actual vs Predicted\nR² = {r2_score(y["Clearance_B12"], y_pred[:, 1]):.4f}')
    plt.xlabel('Actual Clearance (mL/min)')
    plt.ylabel('Predicted Clearance (mL/min)')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'actual_vs_predicted_b12.png'), dpi=300)
    plt.close()

    # 4. Error Distribution
    errors_urea = y['Clearance_Urea'] - y_pred[:, 0]
    errors_b12 = y['Clearance_B12'] - y_pred[:, 1]
    
    plt.figure(figsize=(8, 5))
    sns.histplot(errors_urea, kde=True, color='blue', label='Urea Error', alpha=0.5, bins=30)
    sns.histplot(errors_b12, kde=True, color='green', label='VitB12 Error', alpha=0.5, bins=30)
    plt.title('Distribution of Prediction Errors (Residuals)')
    plt.xlabel('Error (Actual - Predicted) (mL/min)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'error_distribution.png'), dpi=300)
    plt.close()
    
    print(f"Plots successfully generated in {output_dir}")

if __name__ == '__main__':
    create_model_plots()
