use dame::decollapse::{run, DecollapseArgs};
use std::io::Write;
use tempfile::tempdir;

#[test]
fn test_decollapse_single_seq() {
    let dir = tempdir().unwrap();
    let input_path = dir.path().join("input.txt");
    let output_path = dir.path().join("out.fasta");

    let mut f = std::fs::File::create(&input_path).unwrap();
    writeln!(f, "CO1\tTag1\tTag2\t3\tACGT").unwrap();

    run(DecollapseArgs {
        input: input_path.to_str().unwrap().to_string(),
        out_fas: output_path.to_str().unwrap().to_string(),
    })
    .unwrap();

    let contents = std::fs::read_to_string(&output_path).unwrap();
    let lines: Vec<&str> = contents.lines().collect();
    let headers: Vec<&str> = lines.iter().copied().filter(|l| l.starts_with('>')).collect();
    let seqs: Vec<&str> = lines.iter().copied().filter(|l| !l.starts_with('>')).collect();

    assert_eq!(headers.len(), 3);
    assert!(seqs.iter().all(|s| *s == "ACGT"));
    assert!(headers[0].contains("Tag1.Tag2.3_1"));
    assert!(headers[1].contains("Tag1.Tag2.3_2"));
    assert!(headers[2].contains("Tag1.Tag2.3_3"));
}

#[test]
fn test_decollapse_multiple_seqs() {
    let dir = tempdir().unwrap();
    let input_path = dir.path().join("input.txt");
    let output_path = dir.path().join("out.fasta");

    let mut f = std::fs::File::create(&input_path).unwrap();
    writeln!(f, "CO1\tTag1\tTag2\t2\tAAAA").unwrap();
    writeln!(f, "CO1\tTag1\tTag2\t1\tCCCC").unwrap();

    run(DecollapseArgs {
        input: input_path.to_str().unwrap().to_string(),
        out_fas: output_path.to_str().unwrap().to_string(),
    })
    .unwrap();

    let contents = std::fs::read_to_string(&output_path).unwrap();
    let headers: Vec<&str> = contents.lines().filter(|l| l.starts_with('>')).collect();

    assert_eq!(headers.len(), 3); // 2 + 1
}

#[test]
fn test_decollapse_global_seq_id_increments() {
    let dir = tempdir().unwrap();
    let input_path = dir.path().join("input.txt");
    let output_path = dir.path().join("out.fasta");

    let mut f = std::fs::File::create(&input_path).unwrap();
    writeln!(f, "CO1\tTag1\tTag2\t2\tAAAA").unwrap();
    writeln!(f, "CO1\tTag3\tTag4\t2\tCCCC").unwrap();

    run(DecollapseArgs {
        input: input_path.to_str().unwrap().to_string(),
        out_fas: output_path.to_str().unwrap().to_string(),
    })
    .unwrap();

    let contents = std::fs::read_to_string(&output_path).unwrap();
    let headers: Vec<&str> = contents.lines().filter(|l| l.starts_with('>')).collect();

    assert_eq!(headers.len(), 4);
    assert!(headers[0].contains("Tag1.Tag2.2_1"));
    assert!(headers[1].contains("Tag1.Tag2.2_2"));
    assert!(headers[2].contains("Tag3.Tag4.2_3"));
    assert!(headers[3].contains("Tag3.Tag4.2_4"));
}
