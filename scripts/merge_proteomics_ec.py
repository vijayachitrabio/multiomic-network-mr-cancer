import pandas as pd
import os
import gc

def main():
    print("Loading Endometrial Cancer data...")
    ec_path = "EC_data_final.csv"
    ec_df = pd.read_csv(ec_path)
    
    ec_ids = set(ec_df['Participant ID'].dropna().astype(int))
    print(f"Total participants in EC data: {len(ec_ids)}")
    
    olink_path = "/Users/vijayachitramodhukur/Library/Mobile Documents/com~apple~CloudDocs/ECLAI/uterine_fibroids/olink_data.csv"
    output_path = "EC_data_with_proteomics.csv"
    
    print("Loading Olink Proteomics data (this may take a moment)...")
    # Using chunksize to avoid memory overload just in case, though 24GB RAM should handle 1GB CSV fine
    chunk_list = []
    chunk_size = 10000
    for chunk in pd.read_csv(olink_path, chunksize=chunk_size):
        # Filter chunk to only include IDs present in EC data
        filtered_chunk = chunk[chunk['eid'].isin(ec_ids)]
        chunk_list.append(filtered_chunk)
        
    olink_filtered = pd.concat(chunk_list)
    print(f"Total matched participants in Proteomics data: {len(olink_filtered)}")
    
    del chunk_list
    gc.collect()
    
    print("Merging datasets...")
    merged_df = pd.merge(
        ec_df, 
        olink_filtered, 
        left_on='Participant ID', 
        right_on='eid', 
        how='inner'
    )
    
    print(f"Merged dataset shape: {merged_df.shape}")
    print(f"Saving merged data to {output_path}...")
    merged_df.to_csv(output_path, index=False)
    print("Merge complete!")

if __name__ == "__main__":
    main()
