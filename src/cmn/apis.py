"""Small clients for the UniProt and KEGG REST APIs."""

import csv
import io
import re
import xml.etree.ElementTree as ET

import requests

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
KEGG = "https://rest.kegg.jp"
HUMAN = 9606


def uniprot_human(gene, fields=("accession", "gene_primary", "protein_name", "length", "xref_geneid", "xref_kegg", "go_p")):
    """Reviewed (Swiss-Prot) human UniProt entry of a gene, as a dict of the requested TSV fields."""
    params = {
        "query": f"gene_exact:{gene} AND organism_id:{HUMAN} AND reviewed:true",
        "format": "tsv",
        "fields": ",".join(fields),
    }
    r = requests.get(UNIPROT, params=params, timeout=60)
    r.raise_for_status()
    rows = list(csv.DictReader(io.StringIO(r.text), delimiter="\t"))
    return rows[0] if rows else None


def uniprot_by_entrez(entrez, fields=("accession", "gene_primary", "protein_name", "go_p")):
    """Reviewed human UniProt entries of an Entrez Gene ID (a gene can map to more than one entry)."""
    params = {
        "query": f"xref:geneid-{entrez} AND organism_id:{HUMAN} AND reviewed:true",
        "format": "tsv",
        "fields": ",".join(fields),
    }
    r = requests.get(UNIPROT, params=params, timeout=60)
    r.raise_for_status()
    return list(csv.DictReader(io.StringIO(r.text), delimiter="\t"))


def go_terms(text):
    """Set of GO IDs in a UniProt GO column, e.g. "bone mineralization [GO:0030282]; ..." -> {"GO:0030282", ...}."""
    return set(re.findall(r"GO:\d{7}", text or ""))


def uniprot_features(text):
    """Parse a UniProt TSV feature column into (type, span, note) tuples, e.g. ("CHAIN", "52..100", "Osteocalcin")."""
    pattern = r'([A-Z_]+) (\d+(?:\.\.\d+)?)((?:; /[a-z_]+="[^"]*")*)'
    return [
        (kind, span, (re.search(r'/note="([^"]*)"', attrs) or [None, ""])[1])
        for kind, span, attrs in re.findall(pattern, text)
    ]


def ncbi_lineages(taxids):
    """NCBI Taxonomy lineage (root first, ending with the species) of each taxid, via E-utilities efetch."""
    r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                     params={"db": "taxonomy", "id": ",".join(map(str, taxids)), "retmode": "xml"}, timeout=60)
    r.raise_for_status()
    out = {}
    for taxon in ET.fromstring(r.text).findall("Taxon"):
        taxid = int(taxon.findtext("TaxId"))
        out[taxid] = taxon.findtext("Lineage").split("; ") + [taxon.findtext("ScientificName")]
    return out


def kegg(op, arg):
    """Raw text of a KEGG REST call, e.g. kegg("get", "hsa:632") or kegg("link", "pathway/hsa:632")."""
    r = requests.get(f"{KEGG}/{op}/{arg}", timeout=60)
    r.raise_for_status()
    return r.text


def kegg_pairs(op, arg):
    """KEGG REST output as a list of (left, right) tab-separated pairs."""
    return [tuple(line.split("\t")[:2]) for line in kegg(op, arg).splitlines() if line]


def kegg_id_from_entrez(entrez):
    """KEGG gene ID (e.g. hsa:632) of an Entrez Gene ID."""
    pairs = kegg_pairs("conv", f"hsa/ncbi-geneid:{entrez}")
    return pairs[0][1] if pairs else None
