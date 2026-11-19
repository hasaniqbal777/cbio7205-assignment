"""Q2 Part A (a) - Osteocalcin identifiers from UniProt and KEGG."""

from cmn.apis import kegg, kegg_id_from_entrez, uniprot_features, uniprot_human

SEED = "BGLAP"


def main():
    entry = uniprot_human(SEED, fields=(
        "accession", "id", "gene_primary", "protein_name", "length", "xref_geneid", "xref_kegg", "go_p",
        "ft_signal", "ft_propep", "ft_chain", "ft_mod_res",
    ))
    entrez = entry["GeneID"].strip(";")
    kegg_id = kegg_id_from_entrez(entrez)

    print("(a) Osteocalcin identifiers")
    print(f"  UniProt ID        {entry['Entry']}  ({entry['Entry Name']})")
    print(f"  Protein name      {entry['Protein names']}")
    print(f"  Entrez Gene ID    {entrez}")
    print(f"  Gene name         {entry['Gene Names (primary)']}")
    print(f"  KEGG ID           {kegg_id}  (UniProt cross-reference: {entry['KEGG'].strip(';')})")
    print(f"  Sequence length   {entry['Length']} aa")
    print("  Sequence features:")
    for field in ("Signal peptide", "Propeptide", "Chain", "Modified residue"):
        for kind, span, note in uniprot_features(entry[field]):
            print(f"    {kind:8s} {span:8s} {note}")
    print("  GO biological process:")
    for term in entry["Gene Ontology (biological process)"].split("; "):
        print(f"    {term}")

    print("\n  KEGG entry:")
    for line in kegg("get", kegg_id).splitlines():
        if line.split(" ")[0] in ("ENTRY", "SYMBOL", "NAME", "ORTHOLOGY", "POSITION", "AASEQ"):
            print(f"    {line}")



if __name__ == "__main__":
    main()
