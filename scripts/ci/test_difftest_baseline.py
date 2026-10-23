"""Fail-closed artifact tests: no network, native compiler or extra dependency."""
import copy
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile
from datetime import datetime, timezone
import difftest_baseline as b

BASE = "a" * 40
NOW = datetime(2026, 10, 23, tzinfo=timezone.utc)


def run_metadata():
    repo = {"id": b.REPOSITORY_ID, "full_name": b.REPOSITORY}
    return dict(id=101, run_attempt=1, head_sha=BASE, event="push", head_branch="main",
                path=b.WORKFLOW_PATH, workflow_id=b.WORKFLOW_ID, repository=repo,
                head_repository=repo.copy(), status="completed", conclusion="success",
                run_started_at="2026-10-22T00:00:00Z")


def artifact():
    return dict(id=201, name="difftest-reports", expired=False, size_in_bytes=100,
                expires_at="2027-01-01T00:00:00Z", created_at="2026-10-22T01:00:00Z",
                digest="sha256:" + "b" * 64, workflow_run=dict(id=101, head_sha=BASE,
                head_branch="main", repository_id=b.REPOSITORY_ID, head_repository_id=b.REPOSITORY_ID))


def report():
    return dict(unit="a.b", seeds=16, status="green", cases_total=1, cases_passed=1,
                cases_failed=0, cases_deferred=0, cases_vacuous=0,
                cases=[dict(name="integer", outcome="pass", pass_=True, values_compared=0)])


def valid_report():
    r = report()
    r["pass"] = True
    r["cases"][0]["pass"] = r["cases"][0].pop("pass_")
    return r


def archive(files):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for name, value in files:
            z.writestr(name, value)
    return buf.getvalue()


class BaselineTests(unittest.TestCase):
    def test_event_matrix(self):
        repo = dict(id=b.REPOSITORY_ID, full_name=b.REPOSITORY)
        for event, extra in [("pull_request", {"pull_request": {"base": {"sha": BASE}}}),
                             ("merge_group", {"merge_group": {"base_sha": BASE}}),
                             ("push", {"ref": "refs/heads/main", "before": BASE}),
                             ("workflow_dispatch", {"inputs": {"baseline_sha": BASE}})]:
            p = dict(repository=repo, **extra)
            self.assertEqual(b.resolve_base(event, BASE, p), BASE)
            for wrong in ["", "0" * 40, "b" * 40]:
                with self.subTest(event=event, wrong=wrong), self.assertRaises(b.BaselineError):
                    b.resolve_base(event, wrong, p)
        with self.assertRaises(b.BaselineError):
            b.resolve_base("workflow_dispatch", "", dict(repository=repo, inputs={}))
        self.assertIsNone(b.resolve_base("workflow_dispatch", "", dict(repository=repo, inputs={"diagnostic_only": "true"})))
        self.assertIsNone(b.resolve_base("push", "0" * 40, dict(repository=repo, ref="refs/heads/main", before="0" * 40)))
        with self.assertRaises(b.BaselineError):
            b.resolve_base("push", "", dict(repository=repo, ref="refs/heads/main", before=""))

    def test_run_trust_and_latest_failure(self):
        r = run_metadata()
        self.assertEqual(b.select_run([r], BASE, b.REPOSITORY, b.WORKFLOW_ID), r)
        for key, value in [("head_sha", "b" * 40), ("event", "pull_request"), ("head_branch", "fork"),
                           ("workflow_id", 1), ("path", "nightly.yml")]:
            with self.subTest(key=key), self.assertRaises(b.BaselineError):
                b.select_run([dict(r, **{key: value})], BASE, b.REPOSITORY, b.WORKFLOW_ID)
        for key in ["repository", "head_repository"]:
            wrong = copy.deepcopy(r); wrong[key]["id"] = 1
            with self.assertRaises(b.BaselineError): b.select_run([wrong], BASE, b.REPOSITORY, b.WORKFLOW_ID)
        for status, conclusion in [("in_progress", None), ("completed", "failure")]:
            with self.assertRaises(b.BaselineError):
                b.select_run([r, dict(r, id=102, status=status, conclusion=conclusion)], BASE, b.REPOSITORY, b.WORKFLOW_ID)

    def test_artifact_boundaries(self):
        r, a = run_metadata(), artifact()
        self.assertEqual(b.select_artifact([a], r, BASE, NOW), a)
        for items in [[], [a, a], [dict(a, expired=True)], [dict(a, digest="")],
                      [dict(a, expires_at="2020-01-01T00:00:00Z")], [dict(a, created_at="2020-01-01T00:00:00Z")]]:
            with self.assertRaises(b.BaselineError): b.select_artifact(items, r, BASE, NOW)
        for key in a["workflow_run"]:
            wrong = copy.deepcopy(a); wrong["workflow_run"][key] = "wrong"
            with self.assertRaises(b.BaselineError): b.select_artifact([wrong], r, BASE, NOW)

    def test_coverage_bridge_and_corrupt_reports(self):
        r = valid_report()
        self.assertEqual(b.validate_report(r, "a.b")["with_values"], 1)
        for field in ["cases_vacuous", "cases"]:
            wrong = copy.deepcopy(r); wrong.pop(field)
            with self.assertRaises(b.BaselineError): b.validate_report(wrong, "a.b")
        for mutation in [lambda x: x["cases"].append(x["cases"][0]), lambda x: x.update(cases_total=3),
                         lambda x: x.update(unit="wrong"), lambda x: x.update(status="fail"),
                         lambda x: x.update(**{"pass": False}), lambda x: x.update(harness_error="broken"),
                         lambda x: x.update(seeds=2**64), lambda x: x["cases"][0].update(values_compared=-1),
                         lambda x: x["cases"][0].update(outcome="known_mismatch")]:
            wrong = copy.deepcopy(r); mutation(wrong)
            with self.assertRaises(b.BaselineError): b.validate_report(wrong, "a.b")
        with self.assertRaises(b.BaselineError):
            json.loads('{"a":1,"a":2}', object_pairs_hook=b.unique_object)
        with self.assertRaises(b.BaselineError): b.validate_set({"a.b": r}, ["a.b", "missing"])

    def test_atomic_zip_and_digest(self):
        blob = archive([("a.b.json", json.dumps(valid_report()))])
        provenance = dict(digest="sha256:" + hashlib.sha256(blob).hexdigest())
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "baseline"
            b.install_zip(blob, out, ["a.b"], provenance)
            self.assertTrue((out / "reports/a.b.json").is_file())
            self.assertEqual(json.loads((out / "provenance.json").read_text())["units"], ["a.b"])
            with self.assertRaises(b.BaselineError): b.install_zip(blob, out, ["a.b"], provenance)
        for bad, expected, prov in [(b"truncated", ["a.b"], provenance),
                                   (blob, ["missing"], provenance),
                                   (blob, ["a.b"], dict(digest="sha256:" + "0" * 64))]:
            with tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / "baseline"
                with self.assertRaises((b.BaselineError, zipfile.BadZipFile)):
                    b.install_zip(bad, out, expected, prov)
                self.assertFalse(out.exists()); self.assertEqual(list(Path(tmp).iterdir()), [])
        for files in [[("../escape.json", "{}")], [("a.b.json", "{}"), ("a.b.json", "{}")], [("other.txt", "x")]]:
            bad = archive(files)
            with tempfile.TemporaryDirectory() as tmp, self.assertRaises(b.BaselineError):
                b.install_zip(bad, Path(tmp) / "out", ["a.b"], dict(digest="sha256:" + hashlib.sha256(bad).hexdigest()))

    def test_modern_counts_and_explicit_nightly_selection(self):
        r = valid_report(); r["schema_version"] = 1
        with self.assertRaises(b.BaselineError): b.validate_report(r, "a.b")
        r["cases"][0]["exact_compared"] = 7
        self.assertEqual(b.validate_report(r, "a.b")["with_values"], 1)
        self.assertEqual(r["cases"][0]["exact_compared"], 7)
        import importlib.util
        spec = importlib.util.spec_from_file_location("fuzz_units", Path(__file__).with_name("fuzz-units.py"))
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        plans = [dict(id="a.b", status="green"), dict(id="audit.recursion", status="green"), dict(id="other", status="todo")]
        self.assertEqual(module.select(plans, ["a.b"]), ["a.b"])
        with self.assertRaises(ValueError): module.select(plans, [])
        with self.assertRaises(ValueError): module.select([], ["a.b"])

    def test_main_installs_only_after_complete_trust_chain(self):
        import contextlib
        import os
        import sys
        from unittest.mock import patch
        blob = archive([("a.b.json", json.dumps(valid_report()))])
        art = artifact(); art.update(digest="sha256:" + hashlib.sha256(blob).hexdigest(),
                                    size_in_bytes=len(blob), expires_at="9999-01-01T00:00:00Z")
        class API(b.GitHub):
            def json(self, path):
                if path.endswith("/ci.yml"): return dict(id=b.WORKFLOW_ID, path=b.WORKFLOW_PATH)
                if "/runs?" in path: return {"workflow_runs": [run_metadata()]}
                if path.endswith("/artifacts?per_page=100&page=1"): return {"artifacts": [art]}
                if path.endswith("/artifacts/201"): return art
                if path.endswith("/runs/101"): return run_metadata()
                raise AssertionError(path)
            def raw(self, path):
                assert path.endswith("/artifacts/201/zip"), path
                return blob
        with tempfile.TemporaryDirectory() as tmp:
            event = Path(tmp) / "event.json"
            event.write_text(json.dumps(dict(repository=dict(id=b.REPOSITORY_ID, full_name=b.REPOSITORY),
                                              pull_request=dict(base=dict(sha=BASE)))))
            out = Path(tmp) / "baseline"
            stdout, stderr = io.StringIO(), io.StringIO()
            args = ["helper", "pull_request", BASE, str(out), "--wait-seconds", "0"]
            with patch.object(b, "GitHub", API), patch.object(b, "base_units", return_value=["a.b"]), \
                 patch.dict(os.environ, GITHUB_EVENT_PATH=str(event)), patch.object(sys, "argv", args), \
                 contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                self.assertEqual(b.main(), 0)
            self.assertEqual(stdout.getvalue().strip(), "--baseline " + str(out / "reports"))
            provenance = json.loads((out / "provenance.json").read_text())
            self.assertEqual((provenance["run_event"], provenance["run_branch"], provenance["base_sha"]), ("push", "main", BASE))
            self.assertEqual(provenance["source_event"], "pull_request")
            # A corrupted download must emit no usable argument and install nothing.
            art["digest"] = "sha256:" + "0" * 64
            args[3] = str(Path(tmp) / "corrupt")
            stdout = io.StringIO()
            with patch.object(b, "GitHub", API), patch.object(b, "base_units", return_value=["a.b"]), \
                 patch.dict(os.environ, GITHUB_EVENT_PATH=str(event)), patch.object(sys, "argv", args), \
                 contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(b.main(), 1)
            self.assertEqual(stdout.getvalue(), "")
            self.assertFalse(Path(args[3]).exists())

    def test_api_pagination_and_errors(self):
        class API(b.GitHub):
            def __init__(self): self.calls = []
            def json(self, path):
                self.calls.append(path)
                return {"items": [1] * 100 if path.endswith("&page=1") else [2]}
        api = API(); self.assertEqual(len(api.pages("path", "items")), 101)
        self.assertEqual(len(api.calls), 2)
        import subprocess
        from unittest.mock import patch
        for code in [403, 404, 410]:
            with self.subTest(status=code), patch("subprocess.check_output", side_effect=subprocess.CalledProcessError(1, ["gh"], stderr=f"HTTP {code}")), self.assertRaises(b.BaselineError):
                b.fetch_exact(b.GitHub(), BASE, 0)
        with patch("subprocess.check_output", side_effect=subprocess.TimeoutExpired(["gh"], 120)), self.assertRaises(b.BaselineError):
            b.GitHub().raw("repos/fixed/path")
        class Pending(b.GitHub):
            def json(self, path):
                if path.endswith("ci.yml"): return dict(id=b.WORKFLOW_ID, path=b.WORKFLOW_PATH)
                return {"workflow_runs": [dict(run_metadata(), status="in_progress", conclusion=None)]}
        with self.assertRaises(b.BaselineError): b.fetch_exact(Pending(), BASE, 0)
        with patch("time.monotonic", side_effect=[0, 0, 0, 2]), patch("time.sleep") as sleep, self.assertRaises(b.BaselineError):
            b.fetch_exact(Pending(), BASE, 1)
        sleep.assert_called_once_with(1)



if __name__ == "__main__": unittest.main()
