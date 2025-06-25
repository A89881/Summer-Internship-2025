import pandas as pd
# Reads data-set and formatts the columnns correctly
# Converts the units for distance vectors to lattice intergers
def format_data(url: str):
    i_df = pd.read_csv(url, sep=r"\s+", index_col=False)
    for col in i_df.columns[2:5]:
        i_df[col] = round(i_df[col] / 2.71)
    i_df = i_df.drop(["Jij"], axis=1)
    i_df.to_csv(r"Data\formatted-data.csv", sep=";", index=False)
    print("Done: The string url is: r'Data \ formatted-data.csv' ")
    return r"Data\formatted-data.csv"