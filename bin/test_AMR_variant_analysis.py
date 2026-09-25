#!/usr/bin/env python3
"""Regression tests for https://github.com/CDCgov/neissflow/issues/37.

An in-frame duplication annotation (e.g. p.Gly70dup) in a monitored AMR gene
made get_FA19_calls() raise KeyError: 'dup', killing the entire AMR variant
analysis step for the isolate. In-frame deletions (p.Gly70del) and complex
indel annotations (e.g. p.Gly70delinsAA) hit the same KeyError through the
same AA-dictionary lookup.

These tests call get_FA19_calls() directly with a minimal fake tab-delimited
snippy VCF (EFFECT column at index 10, as the parser expects). stdlib only.

Run: python3 bin/test_AMR_variant_analysis.py
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import AMR_variant_analysis as amr

AA = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E',
      'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
      'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
      'Tyr': 'Y', 'Val': 'V'}

# Mirrors the real assets/AMR_defaults.tsv row "rplD aa70"
WG_DEFAULTS = {'rplD aa70': {'Gene': 'rplD', 'Locus': 'CP012026',
                             'Nucleotide Position': 'NA',
                             'AA Position': '70', 'Default': 'G'}}


def fake_vcf_row(effect):
    """One tab-delimited snippy row; EFFECT (protein HGVS) at column 10."""
    cols = ['CP012026', '1615178', 'ins', 'G', 'GGCC', 'GGCC:17', 'CDS', '+',
            '210/621', '70/206', effect, 'VT05_01662', 'rplD',
            '50S ribosomal protein L4']
    return '\t'.join(cols) + '\n'


class TestGetFA19CallsIndelAnnotations(unittest.TestCase):
    def run_calls(self, effect):
        with tempfile.NamedTemporaryFile('w', suffix='.tsv',
                                         delete=False) as fh:
            fh.write(fake_vcf_row(effect))
            path = fh.name
        try:
            return amr.get_FA19_calls(WG_DEFAULTS, path, AA)
        finally:
            os.unlink(path)

    def test_dup_annotation_does_not_raise(self):
        # Exact class of annotation from issue #37 (rplD p.Gly70dup)
        results = self.run_calls('conservative_inframe_insertion p.Gly70dup')
        self.assertEqual(results['rplD aa70'], 'Gdup')

    def test_del_annotation_does_not_raise(self):
        results = self.run_calls('conservative_inframe_deletion p.Gly70del')
        self.assertEqual(results['rplD aa70'], 'Gdel')

    def test_complex_indel_annotation_falls_back_to_default(self):
        results = self.run_calls('complex_substitution p.Gly70delinsAA')
        self.assertEqual(results['rplD aa70'], 'G')

    def test_standard_substitution_unchanged(self):
        # Guard against regressing the normal missense path
        results = self.run_calls('missense_variant p.Gly70Arg')
        self.assertEqual(results['rplD aa70'], 'R')

    def test_frameshift_unchanged(self):
        results = self.run_calls('frameshift_variant p.Gly70fs')
        self.assertEqual(results['rplD aa70'], 'G')


if __name__ == '__main__':
    unittest.main(verbosity=2)
