use dame::perpcr::{flag_error, run_per_pcr};
use std::fs;
use std::path::{Path, PathBuf};
use tempfile::tempdir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../tests/fixtures/perpcr").join(name)
}

fn read(p: &Path) -> String {
    fs::read_to_string(p).unwrap()
}

fn comp() -> String {
    fixture("Comparisons_3PCRs.fasta").to_str().unwrap().to_string()
}

#[test]
fn fixture_with_psinfo_matches_expected() {
    let out = tempdir().unwrap();
    let ps = fixture("PSinfo.txt");
    run_per_pcr(&comp(), Some(ps.to_str().unwrap()), 0, None, false, out.path()).unwrap();
    assert_eq!(
        read(&out.path().join("FilteredReads.perpcr.fna")),
        read(&fixture("expected/FilteredReads.perpcr.fna"))
    );
    assert_eq!(read(&out.path().join("PCRinfo.txt")), read(&fixture("expected/PCRinfo.txt")));
    let mut names: Vec<String> = fs::read_dir(out.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().into_string().unwrap())
        .collect();
    names.sort();
    assert_eq!(names, vec!["FilteredReads.perpcr.fna", "PCRinfo.txt"]);
}

#[test]
fn fixture_without_psinfo_matches_expected() {
    let out = tempdir().unwrap();
    run_per_pcr(&comp(), None, 0, None, false, out.path()).unwrap();
    assert_eq!(
        read(&out.path().join("FilteredReads.perpcr.fna")),
        read(&fixture("expected/FilteredReads.perpcr.fna"))
    );
    assert_eq!(
        read(&out.path().join("PCRinfo.txt")),
        read(&fixture("expected/PCRinfo_no_psinfo.txt"))
    );
}

#[test]
fn length_filters_drop_records_without_padding() {
    let out = tempdir().unwrap();
    let ps = fixture("PSinfo.txt");
    run_per_pcr(&comp(), Some(ps.to_str().unwrap()), 100, Some(200), false, out.path()).unwrap();
    let fasta = read(&out.path().join("FilteredReads.perpcr.fna"));
    let lines: Vec<&str> = fasta.lines().collect();
    assert_eq!(lines.len(), 28);
    assert!(lines.iter().skip(1).step_by(2).all(|s| s.len() == 112 || s.len() == 120));
    assert_eq!(lines[26], ">S3_PCR3.14;size=6;sample=S3_PCR3;");
    assert!(read(&out.path().join("PCRinfo.txt")).contains("S3_PCR3\tS3\t3\tt17-t18\t2\t6\n"));
}

#[test]
fn short_header_record_is_skipped() {
    let dir = tempdir().unwrap();
    let src = dir.path().join("in.fasta");
    fs::write(&src, ">S1\nACGT\n>S1\tt1-t2.t3-t4_1\t1_0\nACGT\n").unwrap();
    let out = dir.path().join("out");
    fs::create_dir(&out).unwrap();
    run_per_pcr(src.to_str().unwrap(), None, 0, None, false, &out).unwrap();
    assert_eq!(
        read(&out.join("FilteredReads.perpcr.fna")),
        ">S1_PCR1.1;size=1;sample=S1_PCR1;\nACGT\n"
    );
}

#[test]
fn flag_errors() {
    assert_eq!(
        flag_error(true, None, true),
        Some("--sample-fastas cannot be used with --per-pcr")
    );
    assert_eq!(flag_error(false, Some("PS.txt"), false), Some("--ps-info requires --per-pcr"));
    assert_eq!(flag_error(true, Some("PS.txt"), false), None);
    assert_eq!(flag_error(false, None, true), None);
}

/// (name, fasta, psinfo, message); {psinfo} and {input} are replaced with paths.
const ERROR_CASES: &[(&str, &str, Option<&str>, &str)] = &[
    ("counts_mismatch", ">S1\tt1-t2.t3-t4_1\t1_2_3\nACGT\n", None,
     "record 1 for sample S1 has 3 counts but 2 tag pairs"),
    ("non_integer", ">S1\tt1-t2.t3-t4_1\t1_x\nACGT\n", None,
     "record 1 for sample S1 has a non-integer count 'x'"),
    ("inconsistent_x", ">S1\tt1-t2.t3-t4_1\t1_0\nACGT\n>S2\tt5-t6.t7-t8.t9-t10_1\t1_0_0\nACGT\n", None,
     "record 2 for sample S2 has 3 PCRs but earlier records have 2"),
    ("inconsistent_tags", ">S1\tt1-t2.t3-t4_1\t1_1\nACGT\n>S1\tt1-t2.empty-empty_2\t1_0\nAACC\n", None,
     "sample S1 has tag pairs t1-t2.empty-empty in record 2 but t1-t2.t3-t4 in an earlier record"),
    ("sample_not_in_psinfo", ">S9\tt1-t2.t3-t4_1\t1_0\nACGT\n", Some("S1\tt1\tt2\t1\nS1\tt3\tt4\t1\n"),
     "sample S9 is in the input but not in {psinfo}"),
    ("tag_mismatch", ">S1\tt1-t2.t3-t4_1\t1_0\nACGT\n", Some("S1\tt1\tt2\t1\nS1\tt3\tt9\t1\n"),
     "sample S1 PCR 2: tag pair t3-t4 in the input but t3-t9 in PSinfo"),
    ("psinfo_rows", ">S1\tt1-t2.t3-t4_1\t1_0\nACGT\n", Some("S1\tt1\tt2\t1\n"),
     "sample S1 does not have exactly one PSinfo row for each PCR 1..2"),
    ("empty_input", "", None, "no records in {input}"),
];

#[test]
fn errors_leave_no_output() {
    for (name, fasta, psinfo, message) in ERROR_CASES {
        let dir = tempdir().unwrap();
        let src = dir.path().join("in.fasta");
        fs::write(&src, fasta).unwrap();
        let ps = psinfo.map(|text| {
            let p = dir.path().join("PSinfo.txt");
            fs::write(&p, text).unwrap();
            p
        });
        let out = dir.path().join("out");
        fs::create_dir(&out).unwrap();
        let err = run_per_pcr(
            src.to_str().unwrap(),
            ps.as_ref().map(|p| p.to_str().unwrap()),
            0,
            None,
            false,
            &out,
        )
        .expect_err(name);
        let expected = message
            .replace("{psinfo}", ps.as_ref().map(|p| p.to_str().unwrap()).unwrap_or(""))
            .replace("{input}", src.to_str().unwrap());
        assert_eq!(err.to_string(), expected, "case {}", name);
        assert_eq!(fs::read_dir(&out).unwrap().count(), 0, "case {} left output", name);
    }
}
