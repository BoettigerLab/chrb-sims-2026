"""cohesin_collision = 'share_no_bypass': arms may share a monomer, never pass.

Checks, on dense 1D systems (many extruders on a short chain):
  1. no bypass: while two extruders stay bound, the order of their arms
     never inverts (a sign change of any arm-arm difference that skips or
     crosses zero would be a pass);
  2. head-on arms do end up sharing a monomer;
  3. the sticky rule, by contrast, does let arms pass after a stall releases
     (documents why the new rule exists);
  4. an arm resumes as soon as its blocker unbinds;
  5. CTCF still stalls arms;
  6. check_run_dir accepts folders made before cohesin_collision existed.

Runs under pytest, or standalone: python tests/test_share_collision.py
"""

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
from numba import typed

SIM_DIR = Path(__file__).resolve().parents[1] / "simulation"
sys.path.insert(0, str(SIM_DIR))

from run_sim_LE import SMC, update_SMC_sim  # noqa: E402


def _smcs(n, n_beads=300, lifetime=1e6, crash_lifetime=1000.0, ctcf=()):
    out = typed.List()
    for i in range(n):
        out.append(SMC(
            N_beads=n_beads,
            bound_lifetime=lifetime,
            unbound_lifetime=1,
            SMC_crash_lifetime=crash_lifetime,
            SMC_crash_prob=1.0,
            extrusion_sided=2,
            extrusion_rate=1.0,
            extrusion_rate_sd=0.0,
            cohesin_loading_probability_list=np.full(n_beads - 3, 1.0 / (n_beads - 3)),
            ctcf_site_location_list=np.array(ctcf, dtype=np.int64),
            ctcf_site_direction_list=np.zeros(len(ctcf), dtype=np.int64),
            ctcf_site_stall_probability_list=np.ones(len(ctcf)),
            ctcf_site_stall_time_list=np.full(len(ctcf), 1e9),
            smc_id=i,
        ))
    return out


def _run(smcs, steps, share):
    return np.stack([update_SMC_sim(smcs, share).copy() for _ in range(steps)])


def _passes(pos):
    """Count arm-order inversions between extruders that stay bound and move
    continuously (|step| <= 1 per arm, i.e. no unbind/rebind in between)."""
    n = 0
    T, K, _ = pos.shape
    for t in range(1, T):
        a, b = pos[t - 1], pos[t]
        ok = (a >= 0).all(1) & (b >= 0).all(1) & (np.abs(b - a) <= 1).all(1)
        for i in range(K):
            for j in range(i + 1, K):
                if not (ok[i] and ok[j]):
                    continue
                for ai in range(2):
                    for aj in range(2):
                        d0 = a[i, ai] - a[j, aj]
                        d1 = b[i, ai] - b[j, aj]
                        if d0 * d1 < 0:
                            n += 1
    return n


def test_share_never_bypasses():
    np.random.seed(1)
    pos = _run(_smcs(12), 3000, share=True)
    assert _passes(pos) == 0


def test_share_heads_meet_on_one_monomer():
    np.random.seed(2)
    pos = _run(_smcs(12), 2000, share=True)
    shared = 0
    for t in range(len(pos)):
        r = pos[t, :, 1][pos[t, :, 1] >= 0]
        l = pos[t, :, 0][pos[t, :, 0] >= 0]
        shared += np.intersect1d(r, l).size
    assert shared > 0


def test_sticky_can_bypass_after_release():
    # Short stall time, so stalls release often: the original rule then lets
    # arms step past each other. This is the behaviour the new rule removes.
    np.random.seed(3)
    pos = _run(_smcs(12, crash_lifetime=5.0), 3000, share=False)
    assert _passes(pos) > 0


def test_share_resumes_when_blocker_unbinds():
    np.random.seed(4)
    smcs = _smcs(2, n_beads=200)
    # place them by hand, head-on and touching: A = [50, 100], B = [100, 150]
    for s, (l, r) in zip(smcs, [(50, 100), (100, 150)]):
        s.bound = True
        s.l_pos, s.r_pos = l, r
        s.current_lifetime = 1e6
    update_SMC_sim(smcs, True)
    assert smcs[0].r_pos == 100 and smcs[1].l_pos == 100   # stuck, sharing 100
    smcs[1].unbind()
    smcs[1].bound_lifetime = 0                               # never rebinds
    update_SMC_sim(smcs, True)
    assert smcs[0].r_pos == 101                              # moving again at once


def test_share_ctcf_still_stalls():
    np.random.seed(5)
    pos = _run(_smcs(3, ctcf=(150,)), 2000, share=True)
    arms = pos[pos >= 0]
    at_site = np.mean(arms == 150)
    assert at_site > 20 / 300   # far above the uniform 1/300


def test_check_run_dir_accepts_pre_field_folders():
    from config_loader import check_run_dir
    from sim_config import SimConfig

    cfg = SimConfig(attraction_energy=0)
    with tempfile.TemporaryDirectory() as tmp:
        run_dir = Path(tmp)
        old = json.loads(cfg.model_dump_json())
        old.pop("cohesin_collision")               # as written by older code
        (run_dir / "config_resolved.json").write_text(json.dumps(old))
        (run_dir / "0").mkdir()
        (run_dir / "0" / "blocks_0-0.h5").write_bytes(b"")
        check_run_dir(cfg, run_dir)                # same behaviour: accepted
        try:
            check_run_dir(cfg.model_copy(update={"cohesin_collision": "share_no_bypass"}), run_dir)
        except SystemExit:
            pass
        else:
            raise AssertionError("switching the collision rule should refuse the old folder")


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"  PASS {name}")
            except Exception as exc:  # noqa: BLE001
                failed += 1
                print(f"  FAIL {name}: {type(exc).__name__}: {exc}")
    print(f"\n{'FAILED' if failed else 'ALL PASSED'}")
    sys.exit(1 if failed else 0)
