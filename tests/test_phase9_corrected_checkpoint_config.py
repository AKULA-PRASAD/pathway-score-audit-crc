"""Checkpointing, recovery, configuration fail-fast, H4 and output-dir isolation.
All writes go to pytest's tmp_path; output_dir fixtures inside the repository are only
validated, never written."""
import json
from pathlib import Path

import pytest

from src.phase9_corrected.checkpoint import (CheckpointCorruption, CheckpointError,
                                             CheckpointStore, DuplicateReplicateError,
                                             LockHeldError)
from src.phase9_corrected.config import (ConfigurationError, CorrectedConfig, OWNER_FIELDS,
                                         assert_output_dir_allowed, template_config)
from src.phase9_corrected.data import LAM_GRID
from src.phase9_corrected.engine import compute_outer_replicate
from src.phase9_corrected.integrity import (FROZEN_INPUTS, HardFailureH4, check_hashes,
                                            h4_check, historical_preservation_check,
                                            write_manifest_once)
from src.phase9_corrected.orchestration import run_corrected
from src.phase9_corrected.rng import PlanStreams
from src.phase9_corrected.config import PLAN_BASE
from tests.p9c_support import FakePipeline, approved_test_config, synthetic_context

ROOT = Path(__file__).resolve().parents[1]


def store(tmp_path, run_id=0, cfg="c" * 64, man="m" * 64):
    # tmp_path is outside the repository: use it as its own "root" for isolation checks
    return CheckpointStore(tmp_path / "out", run_id, cfg, man, root=tmp_path)


def recs(n, run_id=0):
    ctx = synthetic_context()
    s = PlanStreams(PLAN_BASE, run_id, "independent", 7)
    return [compute_outer_replicate(ctx, s, FakePipeline(), run_id, k, 3, LAM_GRID)
            for k in range(n)]


# ------------------------------------------------------------------ checkpoint
def test_create_append_reload_order_and_hash(tmp_path):
    rr = recs(5)
    with store(tmp_path) as st:
        for r in reversed(rr):                       # completion order != k order
            st.append(r)
    st2 = store(tmp_path).open()
    assert [r["k"] for r in st2.ordered(5)] == [0, 1, 2, 3, 4]
    assert st2.ordered(5)[2] == json.loads(json.dumps(rr[2]))
    with pytest.raises(CheckpointError):
        st2.ordered(6)


def test_restart_reproducibility_same_bytes(tmp_path):
    rr = recs(4)
    with store(tmp_path / "a") as st:
        for r in rr:
            st.append(r)
    with store(tmp_path / "b") as st:               # interrupted after 2, then resumed
        for r in rr[:2]:
            st.append(r)
    with store(tmp_path / "b") as st:
        todo = [k for k in range(4) if k not in st.done]
        assert todo == [2, 3]
        for k in todo:
            st.append(rr[k])
    a = (tmp_path / "a/out/runs/run0.jsonl").read_bytes()
    b = (tmp_path / "b/out/runs/run0.jsonl").read_bytes()
    assert a == b


def test_duplicate_prevention(tmp_path):
    r = recs(1)[0]
    with store(tmp_path) as st:
        st.append(r)
        with pytest.raises(DuplicateReplicateError):
            st.append(r)
    # duplicate planted on disk is detected on load
    p = tmp_path / "out/runs/run0.jsonl"
    lines = p.read_bytes().split(b"\n")
    p.write_bytes(b"\n".join(lines[:2] + [lines[1]]) + b"\n")
    with pytest.raises(DuplicateReplicateError):
        store(tmp_path).open()


def test_lock_prevents_second_writer(tmp_path):
    with store(tmp_path):
        with pytest.raises(LockHeldError):
            store(tmp_path).acquire()
    with store(tmp_path):                           # released after exit
        pass


def test_append_requires_lock_and_matching_run(tmp_path):
    st = store(tmp_path).open()
    with pytest.raises(CheckpointError):
        st.append(recs(1)[0])
    with store(tmp_path, run_id=1) as st1:
        with pytest.raises(CheckpointError):
            st1.append(recs(1, run_id=0)[0])


def test_partial_tail_detected_and_explicit_recovery_preserves_original(tmp_path):
    rr = recs(3)
    with store(tmp_path) as st:
        for r in rr:
            st.append(r)
    p = tmp_path / "out/runs/run0.jsonl"
    full = p.read_bytes()
    p.write_bytes(full[:-40])                       # simulated crash mid-write
    with pytest.raises(CheckpointCorruption):
        store(tmp_path).open()
    st = store(tmp_path)
    info = st.recover_partial_tail(tmp_path / "quarantine")
    assert Path(info["quarantined_original"]).read_bytes() == full[:-40]
    st.open()
    assert sorted(st.done) == [0, 1]                # last record dropped, not repaired
    with st:
        st.append(rr[2])
    assert p.read_bytes() == full                   # resumed file equals uninterrupted file


def test_tampered_record_and_header_mismatch_detected(tmp_path):
    with store(tmp_path) as st:
        st.append(recs(1)[0])
    p = tmp_path / "out/runs/run0.jsonl"
    s = p.read_text().replace('"k":0,"outer"', '"k":0,"outer"', 1)
    obj = [json.loads(l) for l in s.splitlines()]
    obj[1]["record"]["dC_app"] = 0.5
    p.write_text("\n".join(json.dumps(o, sort_keys=True, separators=(",", ":")) for o in obj) + "\n")
    with pytest.raises(CheckpointCorruption, match="hash mismatch"):
        store(tmp_path).open()
    with pytest.raises(CheckpointCorruption, match="header mismatch"):
        store(tmp_path, cfg="d" * 64).open()


# ------------------------------------------------------------------ config
def test_template_is_not_executable_and_lists_every_owner_item():
    cfg = CorrectedConfig(template_config(), ROOT)
    with pytest.raises(ConfigurationError) as ex:
        cfg.validate_for_execution()
    msg = str(ex.value)
    for f in OWNER_FIELDS:
        assert f in msg
    assert "execution_authorization" in msg


def test_unapproved_value_is_not_accepted(tmp_path):
    c = approved_test_config(ROOT, "scratch_p9c_test_out")
    c["approvals"]["apparent_cv"]["approved"] = False
    with pytest.raises(ConfigurationError, match="apparent_cv: no explicit owner approval"):
        CorrectedConfig(c, ROOT).validate_for_execution()


def test_fully_specified_test_config_validates():
    CorrectedConfig(approved_test_config(ROOT, "scratch_p9c_test_out"), ROOT).validate_for_execution()


@pytest.mark.parametrize("mut,frag", [
    (lambda c: c["apparent_cv"].update(T=0), "apparent_cv.T"),
    (lambda c: c["apparent_cv"].update(T=None), "apparent_cv.T"),
    (lambda c: c.update(base=1), "base"),
    (lambda c: c.update(B_in_max=50), "B_in exceeds"),
    (lambda c: c.update(ep1_option="C"), "ep1_option"),
    (lambda c: c.update(w5={"action": "W5-A", "corrective_rule": None,
                             "definition": "adjudication_D5_as_proposed"}), "corrective_rule"),
    (lambda c: c.update(ladder={"rungs": [400, 750], "beyond_ladder_rule": "stop_and_label_imprecise"}), "ladder.rungs"),
    (lambda c: c.update(ladder={"rungs": [750], "beyond_ladder_rule": None}), "beyond_ladder_rule"),
    (lambda c: c.update(output_dir="reports/phase9/dc_interval"), "protected"),
    (lambda c: c.update(output_dir="reports/phase9/dc_interval/new"), "protected"),
    (lambda c: c.update(output_dir="reports/phase9"), "protected"),
    (lambda c: c.update(output_dir="/tmp/elsewhere"), "outside the repository"),
    (lambda c: c.update(run_ids=[0, 0]), "run_ids"),
])
def test_invalid_values_rejected(mut, frag):
    c = approved_test_config(ROOT, "scratch_p9c_test_out")
    mut(c)
    with pytest.raises(ConfigurationError, match=frag):
        CorrectedConfig(c, ROOT).validate_for_execution()


def test_run_refuses_unapproved_config_before_touching_disk(tmp_path):
    cfg = CorrectedConfig(template_config(), ROOT)
    with pytest.raises(ConfigurationError):
        run_corrected(cfg, 0, 500)


def test_run_refuses_topup_without_authorisation_and_non_ladder_size():
    cfg = CorrectedConfig(approved_test_config(ROOT, "scratch_p9c_test_out"), ROOT)
    with pytest.raises(ConfigurationError, match="authorisation"):
        run_corrected(cfg, 0, 750)
    with pytest.raises(ConfigurationError, match="ladder rung"):
        run_corrected(cfg, 0, 600)
    assert not (ROOT / "scratch_p9c_test_out").exists()


def test_output_dir_guard_direct():
    with pytest.raises(ConfigurationError):
        assert_output_dir_allowed(None, ROOT)
    assert_output_dir_allowed(ROOT / "reports" / "phase9_corrected_candidate", ROOT)


# ------------------------------------------------------------------ H4 / preservation
def test_h4_passes_on_frozen_inputs_and_fails_on_change(tmp_path):
    assert h4_check(ROOT, "test")["status"] == "PASS"
    # a copy of one input with one byte changed must fail
    rel = "data/processed/analysis/endpoint.csv"
    (tmp_path / Path(rel).parent).mkdir(parents=True)
    data = (ROOT / rel).read_bytes()
    (tmp_path / rel).write_bytes(data[:-1] + bytes([data[-1] ^ 1]))
    res = check_hashes({rel: FROZEN_INPUTS[rel]}, tmp_path)
    assert not res["all_ok"]
    with pytest.raises(HardFailureH4):
        h4_check(tmp_path, "test")                  # other inputs absent -> also FAIL


def test_historical_artifacts_preserved():
    assert historical_preservation_check(ROOT)["status"] == "PASS"


def test_manifest_write_once(tmp_path):
    p = tmp_path / "m.json"
    h1 = write_manifest_once(p, {"a": 1})
    assert write_manifest_once(p, {"a": 1}) == h1
    with pytest.raises(HardFailureH4):
        write_manifest_once(p, {"a": 2})
