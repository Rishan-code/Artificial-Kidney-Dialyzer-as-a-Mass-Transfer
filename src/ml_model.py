# src/ml_model.py
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import root_mean_squared_error, r2_score, mean_absolute_error, mean_absolute_percentage_error
import pickle
import os
import sys

# Ensure imports work from src
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from solver import solve_solute_profile
from parameters import SOLUTES

def generate_synthetic_data(n_samples=1000):
    """
    Generates a dataset by randomly sampling dialyzer operating conditions
    and using the rigorous ODE solver to calculate the Clearance.
    """
    print(f"Generating synthetic dataset with {n_samples} samples using ODE solver...")
    
    data = []
    
    for _ in range(n_samples):
        # Randomly sample operational parameters bounds simulating real scenarios
        Qb_ml_min = np.random.uniform(200, 500)   # Blood flow: 200 - 500 mL/min
        Qd_ml_min = np.random.uniform(300, 800)   # Dialysate flow: 300 - 800 mL/min
        Quf_ml_min = np.random.uniform(0, 20)     # Ultrafiltration: 0 - 20 mL/min
        
        # Convert to m^3/s for the solver
        qb = Qb_ml_min * (1e-6 / 60)
        qd = Qd_ml_min * (1e-6 / 60)
        quf = Quf_ml_min * (1e-6 / 60)
        
        try:
            # Solve for Urea and VitB12 (Albumin is trivial as clearance is 0)
            res_urea = solve_solute_profile("Urea", qb, qd, quf)
            res_b12 = solve_solute_profile("VitB12", qb, qd, quf)
            
            data.append({
                'Qb_ml_min': Qb_ml_min,
                'Qd_ml_min': Qd_ml_min,
                'Quf_ml_min': Quf_ml_min,
                'Clearance_Urea': res_urea['Clearance_ml_min'],
                'Clearance_B12': res_b12['Clearance_ml_min']
            })
        except Exception as e:
            # Skip if BVP fails to converge (rare with these bounds)
            continue
            
        if len(data) % 100 == 0:
            print(f"  ... {len(data)} samples generated")
            
    df = pd.DataFrame(data)
    # Save the dataset to data/ folder
    os.makedirs(os.path.join(os.path.dirname(__file__), '..', 'data'), exist_ok=True)
    df.to_csv(os.path.join(os.path.dirname(__file__), '..', 'data', 'dialyzer_dataset.csv'), index=False)
    print("Dataset saved to data/dialyzer_dataset.csv")
    return df

def train_surrogate_model(df=None):
    """
    Trains a Random Forest Regressor on the synthetic dataset to predict
    Urea and VitB12 clearances.
    """
    if df is None:
        csv_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'dialyzer_dataset.csv')
        if not os.path.exists(csv_path):
            df = generate_synthetic_data(1000)
        else:
            df = pd.read_csv(csv_path)

    print("\nTraining Random Forest Surrogate ('Digital Twin') Model...")
    X = df[['Qb_ml_min', 'Qd_ml_min', 'Quf_ml_min']]
    y = df[['Clearance_Urea', 'Clearance_B12']]
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Random Forest is highly effective for fast non-linear physiological regressions
    rf_model = RandomForestRegressor(n_estimators=100, random_state=42)
    rf_model.fit(X_train, y_train)
    
    # Predictions
    y_pred = rf_model.predict(X_test)
    
    # Evaluate
    y_true_urea = y_test['Clearance_Urea']
    y_true_b12 = y_test['Clearance_B12']
    pred_urea = y_pred[:, 0]
    pred_b12 = y_pred[:, 1]
    
    r2_urea = r2_score(y_true_urea, pred_urea)
    rmse_urea = root_mean_squared_error(y_true_urea, pred_urea)
    mae_urea = mean_absolute_error(y_true_urea, pred_urea)
    mape_urea = mean_absolute_percentage_error(y_true_urea, pred_urea) * 100
    
    r2_b12 = r2_score(y_true_b12, pred_b12)
    rmse_b12 = root_mean_squared_error(y_true_b12, pred_b12)
    mae_b12 = mean_absolute_error(y_true_b12, pred_b12)
    mape_b12 = mean_absolute_percentage_error(y_true_b12, pred_b12) * 100

    print("\n--- Model Performance ---")
    print("Urea Prediction:")
    print(f"  R^2  : {r2_urea:.4f}")
    print(f"  RMSE : {rmse_urea:.2f} mL/min")
    print(f"  MAE  : {mae_urea:.2f} mL/min")
    print(f"  MAPE : {mape_urea:.2f}%")
    
    print("\nVitB12 Prediction:")
    print(f"  R^2  : {r2_b12:.4f}")
    print(f"  RMSE : {rmse_b12:.2f} mL/min")
    print(f"  MAE  : {mae_b12:.2f} mL/min")
    print(f"  MAPE : {mape_b12:.2f}%")
    
    # Save the model
    model_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'rf_surrogate.pkl')
    with open(model_path, 'wb') as f:
        pickle.dump(rf_model, f)
    print(f"\nModel strictly saved to {model_path}")
    
    return rf_model

if __name__ == "__main__":
    df = generate_synthetic_data(n_samples=1000)
    train_surrogate_model(df)
