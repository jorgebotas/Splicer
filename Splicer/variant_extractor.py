#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import namedtuple, defaultdict
import numpy as np 
import pandas as pd 
from pathlib import Path
from pysam import VariantFile, VariantRecord
import re
import scipy.sparse as sp
from scipy.stats import fisher_exact
from typing import Tuple
from tqdm import tqdm
from joblib import Parallel, delayed

SPLICEAI_EFFECTS = [ "AG", "AL", "DG", "GL" ] 
SpliceAI = namedtuple("SpliceAI",
    "DP_AG DP_AL DP_DG DP_GL DS_AG DS_AL DS_DG DS_DL SYMBOL")

Transcripts = namedtuple("Transcripts", "case ctrl")

def _get_datum(datum):
    """
    Get datum from VEP annotated VCF
    """
    # for fields that could return either direct value 
    # or tuple depending on header
    if type(datum) == tuple:
        return datum[0]
    else:
        return datum


class VariantExtractor:
    VARIANT_COLUMNS = [
        "CSQ", "GENE", "AC_1", "AC_0", "AN_1", "AN_0", "OR", "PVAL", "AC", 
        "AF", "POS", "REF", "ALT", "EXON", "AA_POS", "AA_CHANGE", "AA_LEN", 
        "ENSP", "EA", "SpliceAI_effect", "SpliceAI_delta"
    ]
    def __init__(self,
        vcf_path: Path,
        genes: pd.DataFrame,
        samples: list,
        cases: list, 
        controls: list, 
        min_spliceai: float = 0.5,
        ) -> None:

        self.vcf = vcf_path

        self.genes = genes

        self.samples = samples
        self.cases = cases
        self.controls = controls

        # SpliceAI threshold
        self.min_spliceai = min_spliceai

    def _get_contig(self, gene, vcf):
        contig = f"chr{str(gene.chrom)}" \
                if re.search(r"chr", next(vcf).chrom) \
                else str(gene.chrom)
        return contig

    def _get_ensp(self, c: str):
        return c.split("|")[7]

    def _convert_zygosity(self, genotype: tuple) -> int:
        """
        Convert a genotype tuple to a zygosity integer
        Args:
            genotype (tuple): The genotype of a variant for a sample
        Returns:
            int: The zygosity of the variant (0/1/2)
        """
        if genotype in [(1, 0), (0, 1)]:
            zygo = 1
        elif genotype == (1, 1):
            zygo = 2
        else:
            zygo = 0
        return zygo

    def _encode_genotypes(self, rec: VariantRecord) -> pd.Series:
        return pd.Series([self._convert_zygosity(rec.samples[sample]["GT"])\
                    for sample in self.samples], 
                index=self.samples, dtype=int)

    def _get_ac_case_ctrl(self, gts: pd.Series) -> Tuple[int, int]:
        case_gt = gts[gts.index.isin(self.cases)]
        ctrl_gt = gts[gts.index.isin(self.controls)] 
        ac_case = sum(case_gt)
        ac_ctrl = sum(ctrl_gt)
        return ac_case, ac_ctrl

    def _fisher_test(self, gts: pd.Series) -> Tuple:
        ac_case, ac_ctrl = self._get_ac_case_ctrl(gts)
        an_case = 2 * len(self.cases)
        an_ctrl = 2 * len(self.controls)
        fisher_matrix = [[ ac_case, ac_ctrl ],
                         [ an_case - ac_case, 
                           an_ctrl - ac_ctrl ]]
        odds_ratio, pval = fisher_exact(fisher_matrix)
        return ac_case, ac_ctrl, odds_ratio, pval

    def _aa_change(self, csq_str: str) -> Tuple[int, str, int]:
        """
        Extract amino acid change from CSQ field in VCF
        """
        csq = csq_str.split("|")
        aa_pos, aa_len = csq[5].split("/")
        aa_change = csq[6]
        return int(aa_pos), aa_change, int(aa_len)

    def _validate_ea(self, ea: float) -> float:
        """
        Checks for valid EA score
        Args:
            ea (str/float/None): EA score as string
        Returns:
            float: EA score between 0-100 if valid, otherwise returns NaN
        """
        try:
            ea = float(ea)
        except ValueError:
            if type(ea) == str and (ea == "fs-indel" or "STOP" in ea):
                ea = 100
            else:
                ea = np.nan
        except TypeError:
            ea = np.nan
        return ea


    def _get_ea(self,
        all_ea: tuple, 
        canon_ensp: tuple, 
        all_ensp: tuple, 
        csq: tuple
        ) -> float:
        """
        Return EA for VEP canonical transcript (ENSEMBL)
        """
        if "stop_gained" in csq or\
           "frameshift_variant" in csq or\
           "stop_lost" in csq:
            return 100
        #elif "splice_donor_variant" in csq or "splice_acceptor_variant" in csq:# \
             #or "splice_region_variant" in csq:
            #return 100
        try:
            canon_idx = all_ensp.index(canon_ensp)
        except ValueError:
            return np.nan
        else:
            return self._validate_ea(all_ea[canon_idx])

    def _get_spliceai(self, csq: list) -> SpliceAI:
        """
        Return consequence data associated to canonical transcript
        """

        idx = 10
        
        scores = [ float(s) if s else np.nan for s in csq[idx: idx + 8] ]
        symbol = str(csq[idx + 8])
        return SpliceAI(*scores, symbol)

    def _get_splice_effect(self, spliceai: SpliceAI) -> Tuple[str, float]:
        """
        Return SpliceAI effect with hightest delta score
        AG, DG, AL, DL
        """
        deltas = spliceai[4:-1]
        max_delta = max(deltas) # Get maximum score
        effect = SPLICEAI_EFFECTS[deltas.index(max_delta)]
        return effect, max_delta

    def _check_spliceai(self, spliceai: SpliceAI) -> bool:
        """
        Return whether any SpliceAI delta score surpasses predefined threshold 
        `self.min_spliceai`
        """
        return spliceai is not None and (
                spliceai.DS_AG >= self.min_spliceai or \
                spliceai.DS_AL >= self.min_spliceai or \
                spliceai.DS_DG >= self.min_spliceai or \
                spliceai.DS_DL >= self.min_spliceai)

    def _process_csq(self, csq_list) -> Tuple[str, float, set]:
        max_spliceai = ("", -1.0, set())
        for csq_str in csq_list:
            csq = csq_str.split("|")
            spliceai = self._get_spliceai(csq)
            splice_effect, splice_delta = self._get_splice_effect(spliceai)

            tx = spliceai.SYMBOL
            if splice_delta == max_spliceai[1]:
                max_spliceai[2].add(tx)
            elif splice_delta > max_spliceai[1]:
                max_spliceai = (splice_effect, splice_delta, set([tx]))

        return max_spliceai


    def _update_gene_variants(self,
            variants: list,
            rec: VariantRecord,
            gene_symbol) -> None:

        # Obtain genotypes and convert ot 0,1,2 encoding
        gts = self._encode_genotypes(rec)
        
        # Fisher's exact test
        ac_case, ac_ctrl, odds_ratio, pval = self._fisher_test(gts)

        csq_list = rec.info.get("CSQ", [])

        max_spliceai = self._process_csq(csq_list)

        variants.append([
            _get_datum(rec.info["Consequence"]),
            gene_symbol,
            ac_case,
            ac_ctrl,
            2 * len(self.cases),
            2 * len(self.controls),
            float(odds_ratio),
            pval,
            _get_datum(rec.info["AC"]),
            _get_datum(rec.info["AF"]),
            rec.pos,
            rec.ref,
            _get_datum(rec.alts),
            _get_datum(rec.info["EXON"]),
            "aa_pos",
            "aa_change",
            "aa_len",
            _get_datum(rec.info["ENSP"]), # cannon_ensp
            "ea",
            max_spliceai[0], # Effect
            max_spliceai[1], # Delta
            ",".join(map(str, max_spliceai[2])) # SYMBOLs
        ])
    
    def _extract_gene_variants(self,
          gene: pd.Series,
        ) -> list:
        """

        """
        # VCF file has to be read each time in pickled parallelization
        vcf = VariantFile(str(self.vcf))
        vcf.subset_samples(self.samples)

        variants = []
        contig = self._get_contig(gene=gene, vcf=vcf)
        start, stop = gene.start, gene.end
        for rec in vcf.fetch(contig=contig, start=start, stop=stop):
            gene_symbol = _get_datum(rec.info["SYMBOL"])
            # all_ensp = rec.info.get("Ensembl_proteinid", (canon_ensp,))
            # all_ea = rec.info.get("EA", (None,))

            # Assert the gene is the same as the gene name of interest
            if gene_symbol != gene.name:
                continue

            self._update_gene_variants(variants, rec, gene_symbol)
            
        return variants

    def extract_variants(self) -> pd.DataFrame:
        """

        """
        variants = Parallel(n_jobs=10)(delayed(self._extract_gene_variants)\
                (pd.Series({ "name": idx, **gene })) \
                for idx, gene in tqdm(list(self.genes.iterrows())))

        variant_df = pd.DataFrame(variants, columns=self.VARIANT_COLUMNS)

        return variant_df



class TxEnricher(VariantExtractor):
    VARIANT_COLUMNS = [
        "GENE", "TRANSCRIPT", "SCORE"
    ]
    def __init__(self,
        vcf_path: Path,
        genes: pd.DataFrame,
        samples: list,
        cases: list, 
        controls: list, 
        min_spliceai: float = 0.5,
        ) -> None:

        super().__init__(
            genes=genes,
            vcf_path=vcf_path,
            samples=samples,
            cases=cases,
            controls=controls,
            min_spliceai=min_spliceai)

        self.transcripts = defaultdict(lambda: defaultdict(int))

    def _get_tx_delta(self, ac_delta: int, splice_effect: str) -> int:
        if splice_effect in ("AG", "DG"):
            return ac_delta
        else:
            return - ac_delta

    def _update_gene_transcripts(self, rec: VariantRecord) -> None:
        # Obtain genotypes and convert ot 0,1,2 encoding
        gts = self._encode_genotypes(rec)
        
        # AC for case and control
        ac_case, ac_ctrl = self._get_ac_case_ctrl(gts)

        all_txs = set()
        tx_effects = { eff: set() for eff in SPLICEAI_EFFECTS }
        for csq_str in rec.info.get("CSQ", []):
            csq = csq_str.split("|")
            # Obtain SpliceAI prediction
            spliceai = self._get_spliceai(csq)
            # Effect with maximum delta
            splice_effect, splice_delta = self._get_splice_effect(spliceai)
            # Keep track of transcripts
            all_txs.add(spliceai.SYMBOL)
            # Keep track of transcripts predicted to have a splice variant
            if splice_delta >= self.min_spliceai:
                print(spliceai)
                tx_effects[splice_effect].add(spliceai.SYMBOL)

        gene_symbol = _get_datum(rec.info["SYMBOL"])
        print(gene_symbol, tx_effects)
        ac_delta = ac_case - ac_ctrl
        for splice_effect, txs in tx_effects:
            # If loss of splice site, penalize transcripts
            # If gain of splice site, reward transcripts
            if splice_effect in ("AL", "DL"):
                ac_delta *= -1
            # Update transcripts associated to variant
            for tx in txs:
                self.transcripts[gene_symbol][tx] += ac_delta
            # Update transcripts NOT associated to variant (opposite effect)
            for tx in all_txs:
                if tx not in txs:
                    self.transcripts[gene_symbol][tx] += -ac_delta

    def _gene_tx_enrichment(self,    
          gene: pd.Series,
        ) -> defaultdict:
        """

        """
        # VCF file has to be read each time in pickled parallelization
        vcf = VariantFile(str(self.vcf))
        vcf.subset_samples(self.samples)

        contig = self._get_contig(gene=gene, vcf=vcf)
        start, stop = gene.start, gene.end
        for rec in vcf.fetch(contig=contig, start=start, stop=stop):
            gene_symbol = _get_datum(rec.info["SYMBOL"])

            # Assert the gene is the same as the gene name of interest
            if gene_symbol != gene.symbol:
                continue

            self._update_gene_transcripts(rec)

        print(self.transcripts[str(gene.name)])
        return self.transcripts[str(gene.name)]

    def _format_variants(self) -> pd.DataFrame:
        variants = []
        for gene, txs in self.transcripts.items():
            for tx, score in txs.items():
                variants.append([gene, tx, score])
        return pd.DataFrame(variants, columns=self.VARIANT_COLUMNS)

    def transcript_enrichment(self) -> pd.DataFrame:
        # Parallelize transcript enrichment by gene
        Parallel(n_jobs=10)(delayed(self._gene_tx_enrichment)\
                (pd.Series({ "symbol": idx, **gene })) \
                for idx, gene in tqdm(list(self.genes.iterrows())))
        
        # Return formatted scores
        return self._format_variants()
