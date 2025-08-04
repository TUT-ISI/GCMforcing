import os, sys
import netCDF4 as nc
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, PredictionErrorDisplay, root_mean_squared_error
from sklearn.inspection import permutation_importance, PartialDependenceDisplay
from sklearn.model_selection import learning_curve
import seaborn as sns    
import joblib
import numpy as np
import math

"""
by  Atte Laakso / Aalto University

Evaluates given ML model and saves different metrics and FI values to outputfiles
"""

def evaluate(best_model, x_test, y_test,output_dir,rs,modelsel):
    # evaluate on test set
    y_pred = best_model.predict(x_test)
    
    # ensure that y_test and y_pred are 1D pandas Series
    if isinstance(y_test, pd.DataFrame):
        y_test_series = y_test.iloc[:, 0]
    elif isinstance(y_test, pd.Series):
        y_test_series = y_test
    else:
        y_test_series = pd.Series(y_test)
    y_pred_series = pd.Series(y_pred, index=y_test_series.index)
    
    r2 = r2_score(y_test_series, y_pred_series)
    rmse = root_mean_squared_error(y_test_series, y_pred_series)

    # save metrics to a text file
    with open(os.path.join(output_dir, "metrics.txt"), "w") as f:
        f.write(f"R2: {r2}\n")
        f.write(f"RMSE: {rmse}\n")
    print("Saved evaluation metrics.")

    # plot prediction error as histogram
    fig, ax = plt.subplots(figsize=(6, 5))
    y_test_ar = np.asarray(y_test).flatten() # make these arrays for histogram plot
    y_pred_ar = np.asarray(y_pred).flatten()
    # 2D histogram using log scale
    h = ax.hist2d( y_test_ar, y_pred_ar, bins=50, cmap='viridis', norm=plt.matplotlib.colors.LogNorm())
    cbar = plt.colorbar(h[3], ax=ax) # add colorbar
    cbar.set_label('log10(N points)')
    
    min_val = min(min(y_test_ar), min(y_pred_ar))
    max_val = max(max(y_test_ar), max(y_pred_ar))
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', label='Ideal') # ideal 1:1 line

    # labels and title
    ax.set_xlabel('True Values')
    ax.set_ylabel('Predicted Values')
    ax.set_title(f'2D Histogram of Prediction Error for {modelsel}')
    ax.legend()
    # save
    out_path = os.path.join(output_dir, "prediction_error_hist2d.png")
    plt.savefig(out_path, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Prediction error 2D histogram saved to: {out_path}")

    # calculate permutation feature importances due to consistency
    results = permutation_importance(best_model, x_test, y_test)
    importances = results.importances_mean
    std_importances = results.importances_std
    feature_names = x_test.columns.tolist()

    # create DataFrame with importances and their std deviation for easier handling
    feature_importances = pd.DataFrame({
        'feature': feature_names,
        'importance': importances,
        'std': std_importances
    }).sort_values(by='importance', ascending=False)

    print("Feature importances:")
    print(feature_importances)

    # save to CSV
    feature_importances.to_csv(os.path.join(output_dir, "feature_importances.csv"), index=False)

    # plot results and std as error bars
    plt.figure(figsize=(10, 6))
    plt.barh(
        feature_importances["feature"],
        feature_importances["importance"],
        xerr=feature_importances["std"],
        color=sns.color_palette("viridis", len(feature_importances))
    )
    plt.xlabel("Permutation importance")
    plt.title(f"Feature Importances with Standard Deviation for {modelsel}")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "feature_importances_plot.png"))
    plt.close()

    # form correlation matrix and save it
    correlation_matrix = x_test.corr()
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0, vmin=-1, vmax=1)
    plt.title("Correlation Matrix Heatmap")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "correlation_matrix.png"))
    plt.close()

    print("Making learning curve")
    # learning curve plot
    X = x_test
    y = y_test_series
    
    # train_sizes, train_scores, test_scores = learning_curve(
    #     best_model, X, y, cv=3, scoring='r2', n_jobs=-1,
    #     train_sizes=np.linspace(0.1, 1.0, 4), random_state=rs
    # )
    
    # train_scores_mean = np.mean(train_scores, axis=1)
    # test_scores_mean = np.mean(test_scores, axis=1)

    # plt.figure(figsize=(8, 6))
    # plt.plot(train_sizes, train_scores_mean, 'o-', color="r", label="Training score")
    # plt.plot(train_sizes, test_scores_mean, 'o-', color="g", label="Cross-validation score")
    # plt.xlabel("Training examples")
    # plt.ylabel("R2 score")
    # plt.title(f"Learning Curve for {modelsel}")
    # plt.legend(loc="best")
    # plt.grid(True)
    # plt.tight_layout()
    # plt.savefig(os.path.join(output_dir, "learning_curve.png"))
    # plt.close()

    # make PDP
    common_params = {
        "subsample": 50,
        "n_jobs": 1,
        "grid_resolution": 20,
        "random_state": 0,
    }
    # numbers
    n_features = len(feature_names)
    n_cols = 3
    n_rows = math.ceil(n_features / n_cols)
    # create the correct number of subplots
    fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(12, 4 * n_rows), constrained_layout=True)
    ax = ax.flatten()

    # only pass the number of axes needed
    display = PartialDependenceDisplay.from_estimator(
        best_model,
        X,
        feature_names,
        ax=ax[:n_features],
        **common_params,
    )
    display.figure_.suptitle(
        (
            "Partial dependence of ΔADRE \n"
            f"for the AeroCom dataset with {modelsel}"
        ),
        fontsize=16,
    )
    # save
    plt.savefig(os.path.join(output_dir, "PDP.png"))
    plt.close()

    print("All ready!")
    return
