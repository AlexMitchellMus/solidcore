//! Smoke tests for the self-authored StepOutils corpus.
//!
//! These files sit beside the STEP Tools download but are ours and committed.
//! They exercise parser branches a real corpus reaches: a plain analytic solid,
//! a missing entity reference, a loopless face, a multi-winding loop and an open
//! sheet.

use std::path::Path;

use step_import::part21::{Limits, parse};

#[test]
fn every_stepoutils_case_parses_as_part21() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/stepoutils");
    let mut cases = 0;
    for entry in std::fs::read_dir(&dir).expect("StepOutils directory") {
        let path = entry.expect("entry").path();
        if path.extension().and_then(|ext| ext.to_str()) != Some("stp") {
            continue;
        }
        let bytes = std::fs::read(&path).expect("read case");
        parse(&bytes, &Limits::default())
            .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        cases += 1;
    }
    assert!(cases >= 5, "expected the StepOutils cases, found {cases}");
}
