use ahash::HashMap;
use dame::filter::{
    all_sequences, index_haps, make_ps_num_files, make_sample_name_array, ps_info_rows,
    read_ps_num_files,
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
fn test_index_haps_empty_replicates_have_no_index_or_sequences() {
    let haps: HashMap<usize, Vec<Vec<String>>> = {
        let mut m = HashMap::default();
        m.insert(0, vec![]);
        m.insert(1, vec![]);
        m
    };

    let index = index_haps(2, &haps);
    assert!(index.is_empty());
    assert!(all_sequences(&index).is_empty());
}

#[test]
fn test_index_haps_builds_union_and_preserves_count_behavior() {
    let haps: HashMap<usize, Vec<Vec<String>>> = {
        let mut m = HashMap::default();
        m.insert(
            0,
            vec![
                vec![
                    "CO1".to_string(),
                    "Tag1".to_string(),
                    "Tag2".to_string(),
                    "7".to_string(),
                    "AAAA".to_string(),
                ],
                vec![
                    "CO1".to_string(),
                    "Tag1".to_string(),
                    "Tag2".to_string(),
                    "9".to_string(),
                    "AAAA".to_string(),
                ],
                vec![
                    "CO1".to_string(),
                    "Tag1".to_string(),
                    "Tag2".to_string(),
                    "not-a-count".to_string(),
                    "CCCC".to_string(),
                ],
                vec![
                    "CO1".to_string(),
                    "Tag1".to_string(),
                    "Tag2".to_string(),
                    "-1".to_string(),
                    "GGGG".to_string(),
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
                "TTTT".to_string(),
            ]],
        );
        m
    };

    let index = index_haps(2, &haps);
    assert_eq!(index[&0].forward_tag, "Tag1");
    assert_eq!(index[&0].reverse_tag, "Tag2");
    assert_eq!(index[&0].counts_by_sequence["AAAA"], 7);
    assert_eq!(index[&0].counts_by_sequence["CCCC"], 0);
    assert_eq!(index[&0].counts_by_sequence["GGGG"], -1_i64);
    assert_eq!(
        index[&0]
            .counts_by_sequence
            .get("missing")
            .copied()
            .unwrap_or(0),
        0
    );
    assert_eq!(
        all_sequences(&index),
        ["AAAA", "CCCC", "GGGG", "TTTT"].into_iter().collect()
    );
}

#[test]
fn test_index_haps_short_first_row_uses_empty_tags() {
    let haps: HashMap<usize, Vec<Vec<String>>> = {
        let mut m = HashMap::default();
        m.insert(
            0,
            vec![
                vec!["short".to_string()],
                vec![
                    "CO1".to_string(),
                    "Tag1".to_string(),
                    "Tag2".to_string(),
                    "4".to_string(),
                    "AAAA".to_string(),
                ],
            ],
        );
        m
    };

    let index = index_haps(1, &haps);
    assert_eq!(index[&0].forward_tag, "");
    assert_eq!(index[&0].reverse_tag, "");
    assert_eq!(index[&0].counts_by_sequence["AAAA"], 4);
}

#[test]
fn test_ps_info_rows_assigns_pcr_numbers() {
    let dir = tempdir().unwrap();
    let psinfo = write_psinfo(
        dir.path(),
        &[
            "S1\tt1\tt2\t1",
            "S1\tt3\tt4\t1",
            "S1\tt5\tt6\t2",
            "S2\tt7\tt8\t1",
            "S2\tt9\tt10\t1",
            "S2\tt11\tt12\t2",
        ],
    );
    let rows = ps_info_rows(psinfo.to_str().unwrap(), 3).unwrap();
    let expected: Vec<(String, usize, String, String)> = vec![
        ("S1", 1, "t1-t2", "1"),
        ("S1", 2, "t3-t4", "1"),
        ("S1", 3, "t5-t6", "2"),
        ("S2", 1, "t7-t8", "1"),
        ("S2", 2, "t9-t10", "1"),
        ("S2", 3, "t11-t12", "2"),
    ]
    .into_iter()
    .map(|(s, k, t, p)| (s.to_string(), k, t.to_string(), p.to_string()))
    .collect();
    assert_eq!(rows, expected);
}

#[test]
fn test_ps_info_rows_blank_line_consumes_a_slot() {
    let dir = tempdir().unwrap();
    let psinfo = write_psinfo(dir.path(), &["S1\tt1\tt2\t1", "", "S1\tt3\tt4\t1"]);
    let rows = ps_info_rows(psinfo.to_str().unwrap(), 2).unwrap();
    assert_eq!(rows[0].1, 1);
    assert_eq!(rows[1].1, 1);
    assert_eq!(rows[1].2, "t3-t4");
}

#[test]
fn test_ps_info_rows_skips_short_lines() {
    let dir = tempdir().unwrap();
    let psinfo = write_psinfo(dir.path(), &["S1\tt1\tt2\t1", "S1\tt3", "S1\tt5\tt6\t1"]);
    let rows = ps_info_rows(psinfo.to_str().unwrap(), 3).unwrap();
    assert_eq!(rows.len(), 2);
    assert_eq!((rows[1].1, rows[1].2.as_str()), (3, "t5-t6"));
}
