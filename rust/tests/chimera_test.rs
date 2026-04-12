use dame::chimera_check::{make_fas_seq_one_line, make_tag_files, make_tag_files_with_pools};
use std::io::Write;
use std::sync::Mutex;
use tempfile::tempdir;

// Mutex to serialize tests that change the current directory
static CWD_LOCK: Mutex<()> = Mutex::new(());

#[test]
fn test_make_tag_files_creates_files() {
    let dir = tempdir().unwrap();
    let psinfo_path = dir.path().join("PSinfo.txt");

    {
        let mut f = std::fs::File::create(&psinfo_path).unwrap();
        writeln!(f, "S1\tTag1\tTag2\t1").unwrap();
        writeln!(f, "S1\tTag3\tTag4\t1").unwrap();
    }

    let _lock = CWD_LOCK.lock().unwrap();
    let orig = std::env::current_dir().unwrap();
    std::env::set_current_dir(dir.path()).unwrap();

    make_tag_files(psinfo_path.to_str().unwrap(), 2).unwrap();

    let ps1 = std::fs::read_to_string(dir.path().join("PS1.tags.txt")).unwrap();
    let ps2 = std::fs::read_to_string(dir.path().join("PS2.tags.txt")).unwrap();

    std::env::set_current_dir(orig).unwrap();

    assert!(ps1.contains("Tag1\tTag2"), "PS1 should contain Tag1/Tag2, got: {ps1}");
    assert!(ps2.contains("Tag3\tTag4"), "PS2 should contain Tag3/Tag4, got: {ps2}");
}

#[test]
fn test_make_tag_files_with_pools_includes_pool() {
    let dir = tempdir().unwrap();
    let psinfo_path = dir.path().join("PSinfo.txt");

    {
        let mut f = std::fs::File::create(&psinfo_path).unwrap();
        writeln!(f, "S1\tTag1\tTag2\t1").unwrap();
        writeln!(f, "S1\tTag3\tTag4\t2").unwrap();
    }

    let _lock = CWD_LOCK.lock().unwrap();
    let orig = std::env::current_dir().unwrap();
    std::env::set_current_dir(dir.path()).unwrap();

    make_tag_files_with_pools(psinfo_path.to_str().unwrap(), 2).unwrap();

    let ps1 = std::fs::read_to_string(dir.path().join("PS1.tags.txt")).unwrap();

    std::env::set_current_dir(orig).unwrap();

    assert!(
        ps1.contains("Tag1\tTag2\t1"),
        "PS1 should contain pool number, got: {ps1}"
    );
}

#[test]
fn test_make_fas_seq_one_line() {
    let dir = tempdir().unwrap();

    let fasta_content = ">header1\nACGT\nACGT\n>header2\nTTTT\n";
    std::fs::write(dir.path().join("Pool1.noChim.fasta"), fasta_content).unwrap();

    let _lock = CWD_LOCK.lock().unwrap();
    let orig = std::env::current_dir().unwrap();
    std::env::set_current_dir(dir.path()).unwrap();

    make_fas_seq_one_line(1).unwrap();

    let result = std::fs::read_to_string(dir.path().join("Pool1.noChim.oneLiner.fasta")).unwrap();

    std::env::set_current_dir(orig).unwrap();

    let lines: Vec<&str> = result.lines().collect();
    assert_eq!(lines[0], ">header1");
    assert_eq!(lines[1], "ACGTACGT");
    assert_eq!(lines[2], ">header2");
    assert_eq!(lines[3], "TTTT");
}
