use dame::perpcr::{FILTERED_READS_WARNING, USEARCH_NOTE};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::tempdir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../tests/fixtures/perpcr").join(name)
}

fn dame(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_dame"))
        .current_dir(dir)
        .arg("convert")
        .args(args)
        .output()
        .expect("run dame convert")
}

#[test]
fn per_pcr_writes_expected_outputs() {
    let dir = tempdir().unwrap();
    let comp = fixture("Comparisons_3PCRs.fasta");
    let ps = fixture("PSinfo.txt");
    let out = dame(dir.path(), &["-i", comp.to_str().unwrap(), "--per-pcr", "--ps-info", ps.to_str().unwrap()]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(out.stderr.is_empty());
    assert_eq!(
        fs::read_to_string(dir.path().join("FilteredReads.perpcr.fna")).unwrap(),
        fs::read_to_string(fixture("expected/FilteredReads.perpcr.fna")).unwrap()
    );
    assert_eq!(
        fs::read_to_string(dir.path().join("PCRinfo.txt")).unwrap(),
        fs::read_to_string(fixture("expected/PCRinfo.txt")).unwrap()
    );
}

#[test]
fn per_pcr_warns_for_filteredreads_and_notes_u() {
    let dir = tempdir().unwrap();
    let src = dir.path().join("FilteredReads.fna");
    fs::copy(fixture("FilteredReads.fna"), &src).unwrap();
    let out = dame(dir.path(), &["-i", "FilteredReads.fna", "--per-pcr", "-u"]);
    assert!(out.status.success());
    let stderr = String::from_utf8(out.stderr).unwrap();
    assert_eq!(stderr, format!("{}\n{}\n", FILTERED_READS_WARNING, USEARCH_NOTE));
}

#[test]
fn rejects_sample_fastas_with_per_pcr() {
    let dir = tempdir().unwrap();
    let out = dame(dir.path(), &["-i", "x.fasta", "--per-pcr", "-s"]);
    assert_eq!(out.status.code(), Some(1));
    assert_eq!(
        String::from_utf8(out.stderr).unwrap(),
        "Error: --sample-fastas cannot be used with --per-pcr\n"
    );
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
}

#[test]
fn rejects_ps_info_without_per_pcr() {
    let dir = tempdir().unwrap();
    let out = dame(dir.path(), &["-i", "x.fasta", "--ps-info", "PS.txt"]);
    assert_eq!(out.status.code(), Some(1));
    assert_eq!(String::from_utf8(out.stderr).unwrap(), "Error: --ps-info requires --per-pcr\n");
}

#[test]
fn per_pcr_error_message_and_exit_code() {
    let dir = tempdir().unwrap();
    fs::write(dir.path().join("in.fasta"), ">S1\tt1-t2.t3-t4_1\t1_2_3\nACGT\n").unwrap();
    let out = dame(dir.path(), &["-i", "in.fasta", "--per-pcr"]);
    assert_eq!(out.status.code(), Some(1));
    assert_eq!(
        String::from_utf8(out.stderr).unwrap(),
        "Error: record 1 for sample S1 has 3 counts but 2 tag pairs\n"
    );
}
