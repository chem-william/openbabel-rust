//!
//! Test cases for OpenBabel Wrapper
//!
//! Only support single-thread
//! cargo test -- --test-threads=1

#[cfg(test)]
mod test {
    use crate::ob;
    use tempfile::NamedTempFile;

    const TEST_XYZ: &str = "2

Au      0.020000    0.000000    0.000000
Au      1.442498    2.498480    0.000000";

    #[macro_export]
    macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            assert!(($x - $y).abs() < $d, "abs({} - {}) > {}", $x, $y, $d);
        };
    }

    #[test]
    fn test_mol() {
        let test_data: Vec<(String, (u32, u32, u32, f64))> = vec![
            (String::from("c1ccccc1N"), (7, 7, 7, 93.126)),
            (String::from("CCC(COC(=O)[C@@H](NP(=O)(Oc1ccccc1)OC[C@H]1O[C@@]([C@@H]([C@@H]1O)O)(C#N)c1ccc2n1ncnc2N)C)CC"), (42, 45, 42, 602.576))
        ];
        for (s, r) in test_data.iter() {
            cxx::let_cxx_string!(smiles = s.as_str());
            let mol = ob::OBMol_from_smiles(&smiles);
            assert_eq!(ob::OBMol_num_atoms(&mol), r.0);
            assert_eq!(ob::OBMol_num_bonds(&mol), r.1);
            assert_eq!(ob::OBMol_num_hvy_atoms(&mol), r.2);
            assert_delta!(ob::OBMol_get_mol_wt(&mol), r.3, 1e-3);
        }
    }

    #[test]
    fn test_fingerprint() {
        cxx::let_cxx_string!(smiles = "c1ccccc1");
        let mol = ob::OBMol_from_smiles(&smiles);
        for fp_name in vec![
            "FP2", "FP3", "FP4",
            "ECFP0", "ECFP2", "ECFP4", "ECFP6", "ECFP8", "ECFP10",
            ] {
            cxx::let_cxx_string!(name = fp_name);
            let p_data = ob::OBFingerprint_get_fingerprint(&name, &mol, 4096);
            assert_eq!(p_data.as_ref().unwrap().len(), 128);
        }
    }

    #[test]
    fn test_smarts_pattern() {
        cxx::let_cxx_string!(smiles = "NCC(=O)NCC");
        let mol = ob::OBMol_from_smiles(&smiles);
        let test_data: Vec<(String, (u32, u32, u32, Vec<i32>))> = vec![
            (String::from("O=CN"), (3, 2, 1, vec![4, 3, 5])),
            (String::from("CN"), (2, 1, 3, vec![2, 1, 3, 5, 6, 5])),
        ];

        for (s, (num_atoms, num_bonds, num_match, match_indexes)) in test_data.iter() {
            cxx::let_cxx_string!(smarts = s);
            let sp = ob::OBSmartsPattern_from_smarts(&smarts);
            assert_eq!(ob::OBSmartsPattern_num_atoms(&sp), *num_atoms);
            assert_eq!(ob::OBSmartsPattern_num_bonds(&sp), *num_bonds);
            let match_cxx_vec = ob::OBSmartsPattern_match(&sp, &mol);
            assert_eq!(ob::OBSmartsPattern_num_matches(&sp), *num_match);
            assert_eq!(match_cxx_vec.as_slice(), match_indexes);
        }
    }

    #[test]
    fn test_supported_formats() {
        let supported_formats = ob::OBConversion_get_supported_input_format();
        assert_eq!(
            supported_formats.len(),
            22,
            "wrong amount of supported input formats"
        );

        let supported_formats = ob::OBConversion_get_supported_output_format();
        assert_eq!(
            supported_formats.len(),
            15,
            "wrong amount of supported output formats"
        );
    }

    #[test]
    fn test_read_from_str() {
        let obmol = ob::OBMol_new();
        let obconv = ob::OBConversion_new();
        cxx::let_cxx_string!(input = TEST_XYZ);
        cxx::let_cxx_string!(input_format = "xyz");
        let result = ob::OBConversion_set_in_format(&obconv, &input_format);
        assert_eq!(result, 0, "unable to set input format");

        let result = ob::OBConversion_read_string(&obconv, &obmol, &input);
        assert_eq!(result, 0, "unable to read molecule");

        assert_eq!(ob::OBMol_num_atoms(&obmol), 2, "wrong number of atoms");
    }

    #[test]
    fn test_write_to_str() {
        let obmol = ob::OBMol_new();
        let obconv = ob::OBConversion_new();
        cxx::let_cxx_string!(input = TEST_XYZ);

        // set input and output format separately
        cxx::let_cxx_string!(input_format = "xyz");
        let result = ob::OBConversion_set_in_format(&obconv, &input_format);
        assert_eq!(result, 0, "unable to set input format");

        cxx::let_cxx_string!(output_format = "smi");
        let result = ob::OBConversion_set_out_format(&obconv, &output_format);
        assert_eq!(result, 0, "unable to set output format");

        let result = ob::OBConversion_read_string(&obconv, &obmol, &input);
        assert_eq!(result, 0, "unable to read molecule");

        let result_str = ob::OBConversion_write_string(&obconv, &obmol);
        assert_eq!(
            "[Au][Au]\t\n", result_str,
            "output string did not match expected"
        );
    }

    #[test]
    fn test_read_write_from_file() {
        let file_handle = NamedTempFile::new().expect("unable to create temp file");
        cxx::let_cxx_string!(file_path_cxx = file_handle.path().to_str().unwrap());

        let write_obmol = ob::OBMol_new();
        let obconv = ob::OBConversion_new();
        cxx::let_cxx_string!(input = TEST_XYZ);

        cxx::let_cxx_string!(input_format_cxx = "xyz");
        cxx::let_cxx_string!(output_format_cxx = "xyz");
        let result =
            ob::OBConversion_set_in_and_out_formats(&obconv, &input_format_cxx, &output_format_cxx);
        assert_eq!(result, 0, "unable to set input and output format");

        let result = ob::OBConversion_read_string(&obconv, &write_obmol, &input);
        assert_eq!(result, 0, "unable to read molecule");

        let result = ob::OBConversion_write_file(&obconv, &write_obmol, &file_path_cxx);
        assert_eq!(result, 0, "unable to write to file");

        let read_obmol = ob::OBMol_new();
        let result = ob::OBConversion_read_file(&obconv, &read_obmol, &file_path_cxx);
        assert_eq!(result, 0, "unable to read file");
        assert_eq!(ob::OBMol_num_atoms(&read_obmol), 2, "wrong number of atoms");
    }
}
