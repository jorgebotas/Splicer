import argparse
import pandas as pd 
from pathlib import Path
from .variant_extractor import VariantExtractor, TxEnricher


def path(relative_path: str) -> Path:
    """
    Return path relative to module parent
    """
    return Path(__file__).parent / relative_path

def parse_args() -> argparse.Namespace:
    """
    Parses the arguments.
    """
    parser = argparse.ArgumentParser(description="Splicer")
    parser.add_argument("--vcf", "-v",
                        required=True,
                        help=f"Path to VCF containing genetic information. \
                        Should be annotated")
    parser.add_argument("--samples", "-s",
                        required=True,
                        help="Path to sample tsv/csv. sampleid\tstatus")

    parser.add_argument("--ref", 
                        nargs="?", 
                        default="hg38", 
                        choices=("hg19", "hg38"), 
                        help="Reference human genome")

    parser.add_argument("--noX", 
                        default=False,
                        action="store_true",
                        help="Do not run analysis on chromosome X")

    parser.add_argument("--output", "-o", default=None,
        help="save path for output")


    parser.add_argument("--minSpliceAI", type=float, default=0.5,
        help="SpliceAI threshold")

    parser.add_argument("--cores", type=int, default=1,
        help="number of CPUs to use for multiprocessing")

    # # GENES?
    # parser.add_argument("--gene", "-g", required=True, 
    #     help="gene from which to extract variants")

    return parser.parse_args()



def main():

    args = parse_args()

    ref = pd.DataFrame()
    if args.ref == "hg19":
        ref = pd.read_csv(path("refs/ENSEMBL-lite_GRCh37.v75.txt"), 
                          delimiter="\t", header=0, index_col="gene")
    elif args.ref == "hg38":
        if args.noX:
            ref = pd.read_csv(path("refs/ENSEMBL-lite_GRCh38.v94.noX.txt"), 
                              delimiter="\t", header=0, index_col="gene")
        else:
            ref = pd.read_csv(path("refs/ENSEMBL-lite_GRCh38.v94.txt"), 
                              delimiter="\t", header=0, index_col="gene")

    

    # Sample file
    samples = pd.read_csv(args.samples, header=None, index_col=0, sep="\t")
    # Controls 0
    controls = samples[samples.iloc[:,0] == 0].index.astype(str).tolist()
    # Cases 1
    cases = samples[samples.iloc[:,0] == 1].index.astype(str).tolist()
    total_samples = samples.index.astype(str).tolist()


    
    extractor = TxEnricher(args.vcf, ref, total_samples, cases, controls)
    variants = extractor.transcript_enrichment()
    #variants = extractor.extract_variants()

    if args.output:
        variants.to_csv(args.output, sep="\t", index=False)
        

if __name__ == "__main__":
    main()
