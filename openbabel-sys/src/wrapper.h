#pragma once
#include <memory>
#include "rust/cxx.h"
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>

using namespace std;
namespace OpenBabel {
    class OBMol;
    class OBSmartsPattern;
    class OBConversion;
    class OBForceField;

    // Debug
    void print_global_instances();

    // OBConversion
    // unique_ptr<OBMol> OBConversion_smi_to_mol(const string &smiles);
    unique_ptr<OBConversion> OBConversion_new();
    unsigned int OBConversion_set_in_format(const unique_ptr<OBConversion> & pConv, const string &input_format);
    unsigned int OBConversion_set_out_format(const unique_ptr<OBConversion> & pConv, const string &output_format);
    unsigned int OBConversion_set_in_and_out_formats(const unique_ptr<OBConversion> & pConv, const string &input_format, const string &output_format);
    unsigned int OBConversion_read_string(const unique_ptr<OBConversion> & pConv, const unique_ptr<OBMol> & pMol, const string &input);
    rust::String OBConversion_write_string(const unique_ptr<OBConversion> & pConv, const unique_ptr<OBMol> & pMol);
    unsigned int OBConversion_read_file(const unique_ptr<OBConversion> & pConv, const unique_ptr<OBMol> & pMol, const string &input_path);
    rust::Vec<rust::String> OBConversion_get_supported_input_format();
    rust::Vec<rust::String> OBConversion_get_supported_output_format();


    // OBForceField
    OBForceField* OBForceField_find_forcefield(const string &ff_name);
    unsigned int OBForceField_setup(const unique_ptr<OBMol> & pMol, OBForceField* pFF);
    void OBForceField_conjugate_gradients(OBForceField* pFF, u_int32_t steps, double econv);
    void OBForceField_conjugate_gradients_initialize(OBForceField* pFF, u_int32_t steps, double econv);
    bool OBForceField_conjugate_gradients_take_n_steps(OBForceField* pFF, u_int32_t n);

    void OBForceField_steepest_descent(OBForceField* pFF, u_int32_t steps, double econv);
    void OBForceField_steepest_descent_initialize(OBForceField* pFF, u_int32_t steps, double econv);
    bool OBForceField_steepest_descent_take_n_steps(OBForceField* pFF, u_int32_t n);
    double OBForceField_energy(OBForceField* pFF);


    // OBMol
    unique_ptr<OBMol> OBMol_new();
    unique_ptr<OBMol> OBMol_from_smiles(const string &smiles);
    unsigned int OBMol_num_atoms(const unique_ptr<OBMol> & pMol);
    unsigned int OBMol_num_bonds(const unique_ptr<OBMol> & pMol);
    unsigned int OBMol_num_hvy_atoms(const unique_ptr<OBMol> & pMol);
    double OBMol_get_mol_wt(const unique_ptr<OBMol> & pMol);

    // OBFingerprint
    typedef vector<unsigned int> FPData;
    unique_ptr<FPData> OBFingerprint_get_fingerprint(const string &fp_name, const unique_ptr<OBMol> & pMol, u_int32_t nbits);
    // unique_ptr<FPData> OBFingerprint_get_fingerprint_in_batch(const string &fp_thread_name, const rust::Vec<rust::String> & smiles_vec, u_int32_t nbits);
    // deprecated: slow performance, root cause to be identified

    // OBSmartsPattern
    unique_ptr<OBSmartsPattern> OBSmartsPattern_from_smarts(const string &smarts);
    unsigned int OBSmartsPattern_num_atoms(const unique_ptr<OBSmartsPattern> & pSP);
    unsigned int OBSmartsPattern_num_bonds(const unique_ptr<OBSmartsPattern> & pSP);
    unsigned int OBSmartsPattern_num_matches(const unique_ptr<OBSmartsPattern> & pSP);
    unique_ptr<vector<int>> OBSmartsPattern_match(const unique_ptr<OBSmartsPattern> & pSP, const unique_ptr<OBMol> &pMol);
}
