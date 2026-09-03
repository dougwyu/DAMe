use std::process::Command;

#[test]
fn cli_reports_version_3_1_0() {
    let output = Command::new(env!("CARGO_BIN_EXE_dame"))
        .arg("--version")
        .output()
        .expect("run dame --version");

    assert!(output.status.success());
    assert_eq!(String::from_utf8(output.stdout).unwrap(), "dame 3.1.0\n");
    assert!(output.stderr.is_empty());
}
