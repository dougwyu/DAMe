use ahash::HashMap;
use dame::filter::{
    get_seqs_sets_and_fr_counts, make_ps_num_files, make_sample_name_array, read_ps_num_files,
};
use std::io::Write;
use std::sync::Mutex;
use tempfile::tempdir;

// Mutex to serialize tests that change the current directory
static CWD_LOCK: Mutex<()> = Mutex::new(());

fn write_psinfo(dir: &std::path::Path, lines: &[&str]) -> std::path::PathBuf {
    let path = dir.join("PSinfo.txt");
    let mut f = std::fs::File::create(&path).unwrap();
    for line in lines {
        writeln!(f, "{}", line).unwrap();
    }
    path
}

#[test]
fn test_make_sample_name_array() {
    let dir = tempdir().unwrap();
    let psinfo = write_psinfo(
        dir.path(),
        &[
            "SampleA\tTag1\tTag2\t1",
            "SampleA\tTag3\tTag4\t1",
            "SampleB\tTag5\tTag6\t1",
            "SampleB\tTag7\tTag8\t1",
        ],
    );

    let names = make_sample_name_array(psinfo.to_str().unwrap()).unwrap();
    assert_eq!(names, vec!["SampleA", "SampleB"]);
}

#[test]
fn test_make_sample_name_array_deduplicates() {
    let dir = tempdir().unwrap();
    let psinfo = write_psinfo(
        dir.path(),
        &["S1\tTag1\tTag2\t1", "S1\tTag3\tTag4\t1"],
    );

    let names = make_sample_name_array(psinfo.to_str().unwrap()).unwrap();
    assert_eq!(names.len(), 1);
    assert_eq!(names[0], "S1");
}

#[test]
fn test_make_ps_num_files_creates_files() {
    let dir = tempdir().unwrap();
    let psinfo = write_psinfo(
        dir.path(),
        &["S1\tTag1\tTag2\t1", "S1\tTag3\tTag4\t1"],
    );

    let _lock = CWD_LOCK.lock().unwrap();
    let orig = std::env::current_dir().unwrap();
    std::env::set_current_dir(dir.path()).unwrap();

    make_ps_num_files(psinfo.to_str().unwrap(), 2, 1, false).unwrap();

    let ps1 = std::fs::read_to_string(dir.path().join("PS1_files.txt")).unwrap();
    let ps2 = std::fs::read_to_string(dir.path().join("PS2_files.txt")).unwrap();

    std::env::set_current_dir(orig).unwrap();

    assert!(
        ps1.contains("pool1/Tag1_Tag2.txt"),
        "PS1 should contain pool1/Tag1_Tag2.txt, got: {ps1}"
    );
    assert!(
        ps2.contains("pool1/Tag3_Tag4.txt"),
        "PS2 should contain pool1/Tag3_Tag4.txt, got: {ps2}"
    );
}

#[test]
fn test_read_ps_num_files_returns_lines() {
    let dir = tempdir().unwrap();

    let _lock = CWD_LOCK.lock().unwrap();
    let orig = std::env::current_dir().unwrap();
    std::env::set_current_dir(dir.path()).unwrap();

    // Create PS1_files.txt and PS2_files.txt
    std::fs::write(
        dir.path().join("PS1_files.txt"),
        "pool1/Tag1_Tag2.txt\npool1/Tag5_Tag6.txt\n",
    )
    .unwrap();
    std::fs::write(
        dir.path().join("PS2_files.txt"),
        "pool1/Tag3_Tag4.txt\npool1/Tag7_Tag8.txt\n",
    )
    .unwrap();

    let ps_ins_lines = read_ps_num_files(2).unwrap();

    std::env::set_current_dir(orig).unwrap();

    assert_eq!(ps_ins_lines[&0].len(), 2);
    assert_eq!(ps_ins_lines[&1].len(), 2);
    assert!(ps_ins_lines[&0][0].contains("Tag1_Tag2.txt"));
    assert!(ps_ins_lines[&1][0].contains("Tag3_Tag4.txt"));
}

#[test]
fn test_get_seqs_sets_and_fr_counts_empty() {
    let haps: HashMap<usize, Vec<Vec<String>>> = {
        let mut m = HashMap::default();
        m.insert(0, vec![]);
        m.insert(1, vec![]);
        m
    };

    let (seqs_all, f, r, counts, seqs) = get_seqs_sets_and_fr_counts(2, &haps);

    assert!(seqs_all.is_empty(), "seqs_all should be empty");
    assert!(f.is_empty(), "f should be empty");
    assert!(r.is_empty(), "r should be empty");
    assert!(counts.is_empty(), "counts should be empty");
    assert!(seqs.is_empty(), "seqs should be empty");
}

#[test]
fn test_get_seqs_sets_and_fr_counts_with_data() {
    let haps: HashMap<usize, Vec<Vec<String>>> = {
        let mut m = HashMap::default();
        m.insert(
            0,
            vec![
                vec![
                    "CO1".to_string(),
                    "Tag1".to_string(),
                    "Tag2".to_string(),
                    "3".to_string(),
                    "AAAA".to_string(),
                ],
                vec![
                    "CO1".to_string(),
                    "Tag1".to_string(),
                    "Tag2".to_string(),
                    "1".to_string(),
                    "CCCC".to_string(),
                ],
            ],
        );
        m.insert(
            1,
            vec![vec![
                "CO1".to_string(),
                "Tag1".to_string(),
                "Tag2".to_string(),
                "2".to_string(),
                "AAAA".to_string(),
            ]],
        );
        m
    };

    let (seqs_all, f, r, counts, _seqs) = get_seqs_sets_and_fr_counts(2, &haps);

    assert!(seqs_all.contains("AAAA"), "seqs_all should contain AAAA");
    assert!(seqs_all.contains("CCCC"), "seqs_all should contain CCCC");
    assert_eq!(f[&0], "Tag1");
    assert_eq!(r[&0], "Tag2");
    assert_eq!(counts[&0], vec!["3".to_string(), "1".to_string()]);
}
