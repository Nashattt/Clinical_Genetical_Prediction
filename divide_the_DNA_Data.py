import pandas as pd

# Load your file
df = pd.read_csv('dna.csv')

# # Randomly sample 30% of the rows
subset = df.sample(frac=0.001, random_state=42)

# # Save to a new CSV file

subset.to_csv('dna_60.csv', index=False)
print(df['gene'].unique())