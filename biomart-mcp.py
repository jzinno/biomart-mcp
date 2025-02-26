"""
A MCP server to interface with Biomart.
"""

import sys
import time

import pybiomart

from functools import lru_cache
from mcp.server.fastmcp import FastMCP

DEFAULT_HOST = "http://www.ensembl.org"
MAX_RETRIES = 3
RETRY_DELAY = 2  # seconds

COMMON_ATTRIBUTES = [
    "ensembl_gene_id",
    "external_gene_name",
    "hgnc_symbol",
    "hgnc_id",
    "gene_biotype",
    "ensembl_transcript_id",
    "ensembl_peptide_id",
    "ensembl_exon_id",
    "description",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "band",
    "transcript_start",
    "transcript_end",
    "transcription_start_site",
    "transcript_length",
]

mcp = FastMCP("Biomart")


@lru_cache(maxsize=10)
def get_server():
    """Create and cache a server connection with error handling"""
    try:
        return pybiomart.Server(host=DEFAULT_HOST)
    except Exception as e:
        print(f"Error connecting to Biomart server: {str(e)}", file=sys.stderr)
        raise


@mcp.tool()
def list_marts():
    """
    List all available biomart marts
    Biomart has a hierarchy: MART -> DATASET
    Datasets can be queried for attributes with filters
    """
    try:
        server = get_server()
        return server.list_marts().to_csv(index=False).replace("\r", "")
    except Exception as e:
        print(f"Error listing marts: {str(e)}", file=sys.stderr)
        return f"Error: {str(e)}"


@mcp.tool()
def list_datasets(mart: str):
    """
    List all available biomart datasets for a given mart
    Valid values for mart are:
        - ENSEMBL_MART_ENSEMBL
        - ENSEMBL_MART_MOUSE
        - ENSEMBL_MART_ONTOLOGY
        - ENSEMBL_MART_GENOMIC
        - ENSEMBL_MART_SNP
        - ENSEMBL_MART_FUNCGEN
    """
    try:
        server = get_server()
        return server[mart].list_datasets().to_csv(index=False).replace("\r", "")
    except Exception as e:
        print(f"Error listing datasets for mart {mart}: {str(e)}", file=sys.stderr)
        return f"Error: {str(e)}"


@mcp.tool()
def list_common_attributes(mart: str, dataset: str):
    """List common available biomart attributes for a given dataset"""
    server = pybiomart.Server(host=DEFAULT_HOST)
    df = server[mart][dataset].list_attributes()
    df = df[df["name"].isin(COMMON_ATTRIBUTES)]
    return df.to_csv(index=False).replace("\r", "")


@mcp.tool()
def list_all_attributes(mart: str, dataset: str):
    """
    List all available biomart attributes for a given dataset
    This function is unstable and may cause errors
    Prefer list_common_attributes instead
    """
    server = pybiomart.Server(host=DEFAULT_HOST)
    df = server[mart][dataset].list_attributes()
    df = df[~df["name"].str.contains("_homolog_", na=False)]
    df = df[~df["name"].str.contains("dbass", na=False)]
    df = df[~df["name"].str.contains("affy_", na=False)]
    df = df[~df["name"].str.contains("agilent_", na=False)]
    return df.to_csv(index=False).replace("\r", "")


@mcp.tool()
def list_filters(mart: str, dataset: str):
    """List all available biomart filters for a given dataset"""
    server = pybiomart.Server(host=DEFAULT_HOST)
    return server[mart][dataset].list_filters().to_csv(index=False).replace("\r", "")


@mcp.tool()
def get_data(mart: str, dataset: str, attributes: list[str], filters: dict[str, str]):
    """Get data from biomart"""
    for attempt in range(MAX_RETRIES):
        try:
            server = get_server()
            return (
                server[mart][dataset]
                .query(attributes=attributes, filters=filters)
                .to_csv(index=False)
                .replace("\r", "")
            )
        except Exception as e:
            print(
                f"Error getting data (attempt {attempt+1}/{MAX_RETRIES}): {str(e)}",
                file=sys.stderr,
            )
            if attempt < MAX_RETRIES - 1:
                print(f"Retrying in {RETRY_DELAY} seconds...", file=sys.stderr)
                time.sleep(RETRY_DELAY)
            else:
                return f"Error: {str(e)}"


@mcp.tool()
@lru_cache(maxsize=1000)
def get_translation(mart: str, dataset: str, from_attr: str, to_attr: str, target: str):
    """
    Get translation from one id to another

    Example:
    get_translation("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", "hgnc_symbol", "ensembl_gene_id", "TP53")
    >>> "ENSG00000141510"
    """
    try:
        server = get_server()
        dataset = server[mart][dataset]
        df = dataset.query(attributes=[from_attr, to_attr])
        result_dict = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
        if target not in result_dict:
            print(f"Target '{target}' not found in translation", file=sys.stderr)
            return f"Error: Target '{target}' not found"
        return result_dict[target]
    except Exception as e:
        print(f"Error in translation: {str(e)}", file=sys.stderr)
        return f"Error: {str(e)}"


if __name__ == "__main__":
    mcp.run()
