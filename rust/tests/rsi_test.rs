use dame::rsi::{compare, fmt_rsi, run, RsiArgs};
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

// ── malformed input: too few columns ──────────────────────────────────────────
// `(num_cols - 2)` underflows usize when a line has fewer than two fields,
// producing a huge no_rep that defeated the `no_rep < 2` guard and looped
// effectively forever. Python's floor division yields a negative no_rep and
// exits cleanly, so these pin the Rust side to the same behaviour.

#[test]
fn test_single_column_input_returns_instead_of_hanging() {
    let dir = tempdir().unwrap();
    let input_path = dir.path().join("onecol.txt");
    let output_path = dir.path().join("out.txt");
    std::fs::write(&input_path, "Sample1\n").unwrap();

    let args = RsiArgs {
        input: input_path.to_str().unwrap().to_string(),
        explicit: false,
        output: Some(output_path.to_str().unwrap().to_string()),
    };
    run(args).expect("run should return cleanly on a single-column file");

    // Bails out before writing any output, as Python does via sys.exit(0).
    assert!(
        !output_path.exists(),
        "no output file should be written when there are no replicates"
    );
}

#[test]
fn test_too_few_columns_for_two_replicates() {
    // 4 columns -> no_rep = 1, which is below the two-replicate minimum.
    let dir = tempdir().unwrap();
    let input_path = dir.path().join("narrow.txt");
    let output_path = dir.path().join("out.txt");
    std::fs::write(&input_path, "S1\tt1-t2\t10\tACGT\n").unwrap();

    let args = RsiArgs {
        input: input_path.to_str().unwrap().to_string(),
        explicit: false,
        output: Some(output_path.to_str().unwrap().to_string()),
    };
    run(args).expect("run should return cleanly");
    assert!(!output_path.exists());
}

#[test]
fn test_exactly_two_replicates_still_computed() {
    // 6 columns -> no_rep = 2: the boundary case must still produce output.
    let dir = tempdir().unwrap();
    let input_path = dir.path().join("two.txt");
    let output_path = dir.path().join("out.txt");
    std::fs::write(&input_path, "S1\tt1-t2\t10\tt3-t4\t8\tACGT\n").unwrap();

    let args = RsiArgs {
        input: input_path.to_str().unwrap().to_string(),
        explicit: false,
        output: Some(output_path.to_str().unwrap().to_string()),
    };
    run(args).expect("run should succeed");
    let contents = std::fs::read_to_string(&output_path).unwrap();
    assert!(contents.contains("S1"), "got: {}", contents);
}


// ── fmt_rsi ───────────────────────────────────────────────────────────────────
// RSI values are written to a data file that dame-py also writes, so the text
// must match Python's float repr exactly, not merely parse to the same number.

#[test]
fn test_fmt_rsi_keeps_decimal_point_on_whole_numbers() {
    // Display would print these as "0" and "1"; Python writes "0.0" and "1.0".
    assert_eq!(fmt_rsi(0.0), "0.0");
    assert_eq!(fmt_rsi(1.0), "1.0");
}

#[test]
fn test_fmt_rsi_pads_single_digit_exponents() {
    // Debug would print "1e-5"; Python pads the exponent to two digits.
    assert_eq!(fmt_rsi(1e-5), "1e-05");
    assert_eq!(fmt_rsi(2e-5), "2e-05");
    assert_eq!(fmt_rsi(1e-7), "1e-07");
}

#[test]
fn test_fmt_rsi_leaves_two_digit_exponents_alone() {
    assert_eq!(fmt_rsi(1.1102230246251565e-16), "1.1102230246251565e-16");
    assert_eq!(fmt_rsi(1e-22), "1e-22");
}

#[test]
fn test_fmt_rsi_ordinary_values_unchanged() {
    assert_eq!(fmt_rsi(0.5), "0.5");
    assert_eq!(fmt_rsi(0.0001), "0.0001");
    assert_eq!(fmt_rsi(0.5870445344129555), "0.5870445344129555");
    assert_eq!(fmt_rsi(0.3333333333333333), "0.3333333333333333");
}
