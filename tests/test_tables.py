"""The map editor reads the annotation tables of Human-GEM and of yeast-GEM by column name, from the tables only."""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "maps"))
import mapedit  # noqa: E402

YAML = """\
- metaData:
    short_name: Test-GEM
    version: "1.0.0"
- metabolites:
  - id: s_1
    name: alpha
    compartment: c
  - id: s_2
    name: beta
    compartment: c
- reactions:
  - id: r_1
    name: alpha to beta
    metabolites:
      s_1: -1.0
      s_2: 1.0
    gene_reaction_rule: g1 or g2
    subsystem:
      - Test subsystem
- genes:
  - id: g1
    name: ORF1
  - id: g2
    name: ORF2
- compartments:
    c: cytosol
"""

# Human-GEM style: genes/mets/rxns as first column, columns like rxnKEGGID; quoted header as in Human-GEM 2.0.0
HUMAN = {
    "genes.tsv": '"genes"\t"geneSymbols"\n"g1"\t"SYM1"\ng2\t\n',
    "metabolites.tsv": "mets\tmetKEGGID\tmetHMR2ID\tmetRetired\ns_1\tC00001\tm001\tm001;old1\ns_2\t\t\t\n",
    "reactions.tsv": "rxns\trxnKEGGID\nr_1\tR00001\n",
}
# yeast-GEM style: id as first column, columns named after the annotation keys
YEAST = {
    "genes.tsv": "id\tuniprot\ng1\tP00001\ng2\tP00002\n",
    "metabolites.tsv": "id\tname\tchebi\tkegg.compound\ns_1\talpha\tCHEBI:1\tC00001\ns_2\tbeta\tCHEBI:2\t\n",
    "reactions.tsv": "id\tname\tec-code\tkegg.reaction\nr_1\talpha to beta\t1.1.1.1;1.1.1.2\tR00001\n",
}


def model(tables):
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "test.yml"), "w") as fh:
            fh.write(YAML)
        for name, content in tables.items():
            with open(os.path.join(d, name), "w") as fh:
                fh.write(content)
        return mapedit.Model(os.path.join(d, "test.yml"), d, [])


class TableStyles(unittest.TestCase):
    def test_human_gem_style(self):
        m = model(HUMAN)
        self.assertEqual(m.rxn_kegg, {"r_1": "R00001"})
        self.assertEqual(m.met_kegg, {"s_1": "C00001", "s_2": ""})
        self.assertEqual(m.symbol, {"g1": "SYM1", "g2": "ORF2"})  # symbol from genes.tsv, else the YAML name
        self.assertEqual(m.legacy["m001"], "s_1")  # retired ids resolve to the current id
        self.assertEqual(m.legacy["old1"], "s_1")

    def test_yeast_gem_style(self):
        m = model(YEAST)
        self.assertEqual(m.rxn_kegg, {"r_1": "R00001"})
        self.assertEqual(m.met_kegg, {"s_1": "C00001", "s_2": ""})
        self.assertEqual(m.symbol, {"g1": "ORF1", "g2": "ORF2"})  # gene names from the YAML
        self.assertEqual(m.legacy, {})  # yeast-GEM has no retired id columns

    def test_cross_references_come_from_the_tables_only(self):
        # the YAML carries a KEGG id the table does not have: it is not used
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "test.yml"), "w") as fh:
                fh.write(YAML.replace("    compartment: c\n- reactions:",
                                      "    compartment: c\n    annotation:\n      kegg.compound: C99999\n- reactions:"))
            for name, content in YEAST.items():
                with open(os.path.join(d, name), "w") as fh:
                    fh.write(content.replace("C00001", ""))
            m = mapedit.Model(os.path.join(d, "test.yml"), d, [])
        self.assertEqual(m.met_kegg["s_1"], "")

    def test_missing_or_empty_tables(self):
        m = model({"genes.tsv": "", "reactions.tsv": "id\n"})
        self.assertEqual(m.rxn_kegg, {})
        self.assertEqual(m.symbol, {"g1": "ORF1", "g2": "ORF2"})
        self.assertEqual(m.met_kegg, {})


if __name__ == "__main__":
    unittest.main()
