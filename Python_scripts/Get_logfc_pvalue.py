import os
import pandas as pd
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='''Generate logFC and pValue matrices from edgeR/Deseq2 DEG results 
                                     get DEG results from Calculate_deg.R''')
    parser.add_argument('-g', '--gene_list', required=True, help='File containing list of genes')
    parser.add_argument('-d', '--deg_list', required=True, help='File listing sample names and DEG file paths(sample TAB file path)')
    parser.add_argument('-o', '--output_prefix', default='deg_matrix', help='Prefix for output files')
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    with open(args.gene_list) as f:
        genes = [line.strip() for line in f if line.strip()]
    
    deg_files = []
    with open(args.deg_list) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split()
                sample = parts[0]
                path = parts[1]
                deg_files.append((sample, path))
    
    logfc_df = pd.DataFrame(index=genes)
    pvalue_df = pd.DataFrame(index=genes)
    
    for sample, file_path in deg_files:
        if not os.path.exists(file_path):
            print(f"Warning: File not found: {file_path}")
            continue
            
        try:
            df = pd.read_csv(file_path, sep=',', header=0, comment='#')
            df = df.rename(columns={'Unnamed: 0': 'Gene'})
            df.set_index('Gene', inplace=True)
            
            logfc_df[sample] = df.reindex(genes).get('log2FoldChange', df.get('logFC'))
            pvalue_df[sample] = df.reindex(genes).get('pvalue', df.get('PValue'))
            
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            import traceback
            traceback.print_exc()
    
    logfc_output = f"{args.output_prefix}_logFC.txt"
    pvalue_output = f"{args.output_prefix}_pValue.txt"
    
    logfc_df.to_csv(logfc_output, sep='\t', na_rep='NA')
    pvalue_df.to_csv(pvalue_output, sep='\t', na_rep='NA')
    
    print(f"Successfully created:")
    print(f" - LogFC matrix: {logfc_output}")
    print(f" - pValue matrix: {pvalue_output}")
if __name__ == "__main__":
    main()