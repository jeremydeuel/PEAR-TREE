"""Scored poly-A needs >= 2 independent fragments (tools/rte/annotator.py _supported_tails).

The evidence consensus keeps only columns covered by >= min_independent_fragments pooled
fragments, so a tail read through by one fragment is missing there while the combined.txt.gz
consensus has it (PD37590: no clip at 230/414 ends). The annotator uses the richer junction
string for structure, but scores the tail only as far as >= 2 fragments reach."""
from types import SimpleNamespace

from tools.rte.annotator import RteAnnotator, _richer_junction


def seg(kind, qlen, target="", strand=1):
    return SimpleNamespace(kind=kind, qlen=qlen, target=target, strand=strand)


def lay(frag, segs, role="CLIP"):
    return SimpleNamespace(role=role, frag_key=frag, segments=segs)


def tails(layouts, ev_left="", ev_right="", k=2):
    ann = RteAnnotator.__new__(RteAnnotator)
    ann.cfg = {"rte_polya_min_fragments": k}
    return ann._supported_tails(SimpleNamespace(raw_layouts=layouts), ev_left, ev_right, None, None)


def test_richer_junction_prefers_longer_clip():
    ev = "AAATTTTCCTCA"                       # depth-trimmed: no clip left
    comb = "aaaaaaaaaaaaaaaAAATTTTCC"         # combined consensus still has the tail
    assert _richer_junction(ev, comb) == comb
    assert _richer_junction("aaaaaGGCC", "aaaaGGCC") == "aaaaaGGCC"   # tie / evidence longer
    assert _richer_junction("", comb) == comb


def test_single_fragment_tail_is_not_scored():
    one = [lay(("s1", "f1"), [seg("ELEMENT", 40), seg("POLYA", 20, "A"), seg("REF", 60)]),
           lay(("s1", "f1"), [seg("POLYA", 18, "A"), seg("REF", 80)])]          # same fragment
    assert tails(one) == (0.0, 0.0)


def test_two_fragments_score_the_shorter_tail():
    two = [lay(("s1", "f1"), [seg("POLYA", 25, "A"), seg("REF", 60)]),
           lay(("s2", "f9"), [seg("POLYA", 14, "A"), seg("REF", 60)]),
           lay(("s3", "f3"), [seg("REF", 60), seg("POLYA", 30, "T")])]          # one T-tail only
    assert tails(two) == (14.0, 0.0)


def test_junction_string_layout_is_not_a_fragment():
    js = [lay(("junction",), [seg("POLYA", 30, "A"), seg("REF", 60)], role="JUNCTION"),
          lay(("s1", "f1"), [seg("POLYA", 20, "A"), seg("REF", 60)])]
    assert tails(js) == (0.0, 0.0)


def test_evidence_consensus_tail_counts_as_supported():
    # the evidence consensus is built at >= 2 fragment depth: its tail is proven
    assert tails([], ev_left="ggcaaaaaaaaaaaaAAATTTT") == (12.0, 0.0)
    assert tails([], ev_right="CCTGAtttttttttttcag") == (0.0, 11.0)
