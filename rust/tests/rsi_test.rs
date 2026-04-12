use dame::rsi::{compare, run, RsiArgs};
use ndarray::array;
use std::io::Write;
use tempfile::tempdir;

#[test]
fn test_compare_identical_replicates() {
    // matrix = [[10,10],[20,20],[30,30]]
    // col_a = [10/60, 20/60, 30/60] = col_b
    // min sums = 1.0, RSI = 1 - 1.0 = 0.0
    let matrix = array![[10i64, 10], [20, 20], [30, 30]];
    let result = compare(&matrix, "sample1", 1, 2);
    assert!(
        result.abs() < 1e-10,
        "Expected 0.0 for identical replicates, got {}",
        result
    );
}

#[test]
fn test_compare_completely_different() {
    // matrix = [[100,0],[0,100]]
    // col_a = [1.0, 0.0], col_b = [0.0, 1.0]
    // min sums = 0.0, RSI = 1 - 0.0 = 1.0
    let matrix = array![[100i64, 0], [0, 100]];
    let result = compare(&matrix, "sample1", 1, 2);
    assert!(
        (result - 1.0).abs() < 1e-10,
        "Expected 1.0 for completely different replicates, got {}",
        result
    );
}

#[test]
fn test_compare_partial_overlap() {
    // matrix = [[50,0],[50,100]]
    // col_a = [0.5, 0.5], col_b = [0.0, 1.0]
    // min(0.5,0.0) + min(0.5,1.0) = 0.0 + 0.5 = 0.5
    // RSI = 1 - 0.5 = 0.5
    let matrix = array![[50i64, 0], [50, 100]];
    let result = compare(&matrix, "sample1", 1, 2);
    assert!(
        (result - 0.5).abs() < 1e-10,
        "Expected 0.5 for partial overlap, got {}",
        result
    );
}

#[test]
fn test_compare_zero_replicate_handled() {
    // matrix = [[0,10],[0,20]] — first column all zeros
    // total[0] = 0 → replaced with 1 to avoid div by zero
    // col_a = [0/1, 0/1] = [0.0, 0.0], col_b = [10/30, 20/30]
    // min sums = 0.0, RSI = 1.0 (in [0,1] range)
    let matrix = array![[0i64, 10], [0, 20]];
    let result = compare(&matrix, "sample1", 1, 2);
    assert!(
        result >= 0.0 && result <= 1.0,
        "Expected result in [0, 1], got {}",
        result
    );
}

#[test]
fn test_run_rsi_produces_output() {
    let dir = tempdir().unwrap();
    let input_path = dir.path().join("input.txt");
    let output_path = dir.path().join("RSI_output.txt");

    // Format: sample\tF-R\tcount0\tF-R\tcount1\tseq
    // 2 rows for sample S1, 2 replicates → no_rep = (6-2)/2 = 2
    {
        let mut f = std::fs::File::create(&input_path).unwrap();
        writeln!(f, "S1\tA-B\t10\tA-B\t10\tACGT").unwrap();
        writeln!(f, "S1\tA-B\t20\tA-B\t20\tGGGG").unwrap();
    }

    let args = RsiArgs {
        input: input_path.to_str().unwrap().to_string(),
        explicit: false,
        output: Some(output_path.to_str().unwrap().to_string()),
    };
    run(args).expect("run should succeed");

    assert!(output_path.exists(), "Output file should exist");
    let contents = std::fs::read_to_string(&output_path).unwrap();
    assert!(
        contents.contains("Sample\tRSI"),
        "Output should contain header 'Sample\\tRSI', got: {}",
        contents
    );
    assert!(
        contents.contains("S1"),
        "Output should contain sample S1, got: {}",
        contents
    );
}
