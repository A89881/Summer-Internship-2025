import pandas as pd
# Reads data-set and formatts the columnns correctly
df = pd.read_csv(r"Data\raw-data-set.csv", sep=r"\s+", index_col=False)
# Converts the units for distance vectors to lattice intergers
for col in df.columns[2:5]:
    df[col] = round(df[col] / 2.71)
df = df.drop("Jij", axis=1)
df.to_csv(r"Data\formatted-data.csv", sep=";", index=False)